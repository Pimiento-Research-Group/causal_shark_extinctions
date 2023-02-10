# load packages
library(chronosphere)
library(tidyverse)
library(here)
library(brms)


# plotting configurations
source(here("R", "config_file.R"))

# stage data --------------------------------------------------------------

# load 95 bin Phanerozoic time scale based on the stratigraphic stages 
# of Gradstein et al. 2020
load(here("data", 
          "stages.Rdata"))

# clean up for binning
dat_stages <- stages %>% 
  as_tibble() %>% 
  select(stg, bottom, mid, top)


# fossil data -------------------------------------------------------------


# load PyRate estimates
dat_genus <- read_delim(here("data",
                               "fossil_occurrences",
                               "all_genera_pg_10_Grj_se_est_genus_name.txt")) %>% 
  # estimate fad and lad for each species
  pivot_longer(cols = contains("ts"), 
               names_to = "origination", 
               values_to = "origination_age") %>% 
  pivot_longer(cols = contains("te"), 
               names_to = "extinction", 
               values_to = "extinction_age") %>% 
  group_by(genus) %>% 
  summarise(ori_age = mean(origination_age),
            ext_age = mean(extinction_age)) %>% 
  # bin fad and lad to stages
  mutate(bin_ori = 95 - cut(ori_age, breaks = dat_stages$bottom,
                            include.lowest = TRUE,
                            labels = FALSE), 
         bin_ext = 95 - cut(ext_age, breaks = dat_stages$bottom,
                            include.lowest = TRUE,
                            labels = FALSE)) %>%
  drop_na(bin_ori, bin_ext) %>% 
  # fill in duration bins
  mutate(bin_occ = map2(.x = bin_ori, 
                        .y = bin_ext,
                        .f = ~ seq(.x, .y, by = 1))) %>% 
  select(genus, bin_occ, bin_ext) %>% 
  unnest(bin_occ) %>% 
  # create extinction signal
  group_by(genus) %>% 
  mutate(ext_signal = if_else(bin_occ == bin_ext, 1, 0)) %>% 
  # clean up
  ungroup() %>% 
  select(-bin_ext, 
         bin = bin_occ)



# occurrence data ---------------------------------------------------------

# load occurrence database
dat_occurrences <- read_rds(here("data",
                                 "fossil_occurrences",
                                 "database_occurrences_10_Jan_2023.rds"))

# bin the occurrences to stages
dat_occ_binned <- dat_occurrences %>% 
  mutate(bin_min = 95 - cut(Min_Ma, breaks = dat_stages$bottom,
                            include.lowest = TRUE,
                            labels = FALSE), 
         bin_max = 95 - cut(Max_Ma, breaks = dat_stages$bottom,
                            include.lowest = TRUE,
                            labels = FALSE)) %>% 
  # select entries, where the early and late interval fields indicate 
  # the same stg or the late_intervar field is empty
  filter(bin_min == bin_max | is.na(bin_max)) %>% 
  # in these entries, use the stg indicated by the early_interval
  select(bin = bin_min, everything(), -bin_max) %>% 
  # clean up
  drop_na(bin) %>% 
  # get the age estimate for the corresponding bin
  left_join(dat_stages %>% 
              select(bin = stg, age = mid))



# reconstruct paleo-coordinates -------------------------------------------

# reformat data
dat_coords <- dat_occ_binned %>% 
  distinct(Longitude, Latitude, age) %>% 
  # prepare for reconstruct function
  mutate(age = round(age)) %>% 
  select(Longitude, Latitude, age) %>% 
  as.data.frame() 


# preallocate vector
coord_vector <- data.frame("paleolong" = rep(NA, nrow(dat_coords)), 
                           "paleolat" = rep(NA, nrow(dat_coords)))

# iterate over rows
for (i in 1:nrow(dat_coords)) {
  
  long_lat <- dat_coords[i, 1:2]
  age <- dat_coords[i, 3]
  suppressMessages(
    tryCatch({
      coord_vector[i,] <- reconstruct(long_lat, age = age, chunk = 1)
    }, error = function(e){})
  )
  cat(paste0(i, "  "))
}

# merge back to original data
dat_occ_binned_coord <- dat_occ_binned %>% 
  distinct(Longitude, Latitude, age) %>% 
  select(Longitude, Latitude, age) %>% 
  bind_cols(coord_vector) %>% 
  full_join(dat_occ_binned) %>% 
  drop_na(paleolong, paleolat)  

# # save data
# dat_occ_binned_coord %>% 
#   write_rds(here("data",
#                  "fossil_occurrences",
#                  "binned_rotated_occurrences.rds"), 
#             compress = "gz")



# geographic range --------------------------------------------------------


# calculate latitude range in absolute degrees
dat_range_lat <- dat_occ_binned_coord %>% 
  # select only genera
  filter(rank == "genus") %>% 
  group_by(accepted_name, bin) %>% 
  summarise(min_lat = min(paleolat), 
            max_lat = max(paleolat)) %>% 
  ungroup() %>% 
  mutate(range_lat = abs(min_lat - max_lat)) %>% 
  select(accepted_name, bin, range_lat)


# calculate great circle distance via geodesic distance between two 
# points specified by radian latitude/longitude using the
# Spherical Law of Cosines (slc)
dat_dist <- dat_occ_binned_coord %>% 
  # select only species
  filter(rank == "genus") %>%
  group_by(accepted_name, bin) %>% 
  # get maximum distance
  summarise(min_lat = min(paleolat), 
            max_lat = max(paleolat), 
            min_long = min(paleolong), 
            max_long = max(paleolong)) %>%
  ungroup() %>% 
  # convert to radians
  mutate(across(min_lat:max_long, ~ .x*pi/180)) %>% 
  # use earth radius and slc
  mutate(geo_dist = acos(sin(min_lat)*sin(max_lat) + cos(min_lat)*cos(max_lat) * cos(max_long-min_long)) * 6371) %>% 
  drop_na(geo_dist) %>% 
  select(accepted_name, 
         bin,
         geo_dist)



# latitude ----------------------------------------------------------------

# calculate latitudinal preference of species
dat_latitude <- dat_occ_binned_coord %>% 
  # select only species
  filter(rank == "genus") %>% 
  group_by(accepted_name, bin) %>% 
  summarise(latitude_pref = mean(paleolat)) %>% 
  ungroup()


# calculate abundance -----------------------------------------------------

# number of occurrences
dat_abund <- dat_occ_binned_coord %>% 
  # select only species
  filter(rank == "genus") %>% 
  count(bin, accepted_name) 

# number of species within genera
dat_abund_genus <- dat_occ_binned_coord %>% 
  # select only species
  filter(rank == "species") %>% 
  count(bin, genus) 


# sampling effort ---------------------------------------------------------

# estimate sampling effort by the number of collections per bin
dat_sampling <- dat_occ_binned_coord %>% 
  distinct(bin, collection_no) %>% 
  count(bin) %>% 
  # add bins with zero counts
  mutate(bin = factor(bin, levels = 95:69)) %>% 
  complete(bin, fill = list(n = 0))

# alternatively the number of collections per bin from the pbdb

# download pbdb collections from bin 72 (Hauterivian) to bin 95 (Holocene)
# via the api call
dat_pbdb <- read_csv("https://paleobiodb.org/data1.2/colls/list.csv?interval=Hauterivian,Holocene")

# bin the data
dat_pbdb_sampling <- dat_pbdb %>% 
  mutate(bin_min = 95 - cut(min_ma, breaks = dat_stages$bottom,
                            include.lowest = TRUE,
                            labels = FALSE), 
         bin_max = 95 - cut(max_ma, breaks = dat_stages$bottom,
                            include.lowest = TRUE,
                            labels = FALSE)) %>% 
  # select entries, where the early and late interval fields indicate 
  # the same stg or the late_intervar field is empty
  filter(bin_min == bin_max | is.na(bin_max)) %>% 
  # in these entries, use the stg indicated by the early_interval
  select(bin = bin_min, everything(), -bin_max) %>% 
  # clean up
  drop_na(bin) %>% 
  # count collections per bin 
  distinct(bin, collection_no) %>% 
  count(bin) %>% 
  # add bins with zero counts
  mutate(bin = factor(bin, levels = 95:69)) %>% 
  complete(bin, fill = list(n = 0))



# preservation potential --------------------------------------------------

# get preservation rate from the PyRate output
dat_preservation_raw <- read_delim(here("data",
                                        "raw",
                                        "fossil_occurrences",
                                        "combined_10_mcmc.log"))
# bring in right format
dat_preservation <- dat_preservation_raw %>% 
  select(mean_q, contains("TS")) %>% 
  pivot_longer(cols = contains("TS"), 
               names_to = "species", 
               values_to = "age") %>% 
  # bin to 10 myr
  mutate(bin = 95 - cut(age, breaks = dat_stages$bottom,
                        include.lowest = TRUE,
                        labels = FALSE), 
         .before = 1) %>% 
  # remove unbinned occurrences
  drop_na(bin) %>% 
  # get to genus level
  mutate(species = str_remove(species, "_TS"), 
         accepted_name = str_extract(species, "^[^_]+(?=_)")) %>% 
  group_by(bin, accepted_name) %>% 
  summarise(mean_q = mean(mean_q)) %>% 
  ungroup() %>% 
  select(bin, accepted_name, mean_q) %>% 
  arrange(accepted_name) %>% 
  # use only occurrences from species that could get binned to stages
  filter(accepted_name %in% dat_abund$accepted_name)




# merge and combine -------------------------------------------------------

# combine all datasets to one
dat_full <- dat_range_lat %>% 
  # geographic range
  full_join(dat_dist) %>% 
  full_join(dat_abund) %>% 
  # latitudinal preference
  full_join(dat_latitude) %>% 
  # abundance
  rename(abund = n) %>% 
  # preservation rate
  left_join(dat_preservation) %>% 
  left_join(dat_preservation %>%
              group_by(bin) %>%
              summarise(mean_q_av = mean(mean_q))) %>% 
  mutate(mean_q = if_else(is.na(mean_q), mean_q_av, mean_q)) %>%
  # sampling
  left_join(dat_pbdb_sampling %>% 
              mutate(bin = as.numeric(as.character(bin)))) %>% 
  rename(pbdb_collections = n) %>% 
  left_join(dat_sampling %>% 
              mutate(bin = as.numeric(as.character(bin))) %>% 
              rename(shark_collections = n)) 

# add species extinction signal
dat_full <- dat_full %>% 
  left_join(dat_genus %>%
              rename(accepted_name = genus) %>% 
              filter(accepted_name %in% dat_full$accepted_name)) %>%
  drop_na(ext_signal) %>% 
  # get taxonomy
  left_join(dat_occ_binned %>% 
              distinct(accepted_name, family, order)) %>% 
  # add species counts
  left_join(dat_abund_genus %>% rename(n_genus = n, 
                                       accepted_name = genus)) %>% 
  select(order, family, genus = accepted_name,
         everything())  %>% 
  # make missing values explicit
  replace_na(list(geo_dist = 0, range_lat = 0,  
                  latitude_pref = mean(dat_latitude$latitude_pref), 
                  n_genus = 0))

# save dataset
dat_full %>% 
  write_rds(here("data",
                 "fossil_occurrences",
                 "processed_fossil_data_genus.rds"), 
            compress = "gz")




# scale predictors --------------------------------------------------------


# environmental proxy data
dat_proxy <- read_rds(here("data",
                           "processed_proxy_data.rds"))


# bring predictors on meaningful scales to enable comparisons
dat_proxy <- dat_proxy %>% 
  mutate(across(c(d13C, sr_value, # scale all productivity parameters
                  n_units, outcrop_area), # same for outcrop parameters
                ~ scale(.x)[,1], 
                .names = "{col}_std"))

# same for fossil data
dat_fossil <- dat_full %>% 
  # add zero genus counts
  replace_na(list(n_genus = 0)) %>% 
  mutate(across(c(range_lat, geo_dist, # scale geographic range parameters
                  mean_q, # and preservation rate
                  pbdb_collections, shark_collections), # and sampling effort
                ~ scale(.x)[,1], 
                .names = "{col}_std"), 
         latitude_pref_abs = abs(latitude_pref)) # use absolute latitude 


# combine data sources

# merge
dat_merged <- left_join(dat_proxy,
                        dat_fossil,
                        by = "bin") %>% 
  drop_na(genus) %>% 
  replace_na(list(abund = 0)) %>% 
  select(bin, order, family, genus,
         latitude_pref_abs,
         n_genus,
         abund, 
         ext_signal, 
         sea_level, 
         cont_area, 
         contains("temp"), 
         contains("std")) 



# abundance ---------------------------------------------------------------


# use same adjustment set as in main analysis
# first adjustment set
mod1 <- brm_logistic("ext_signal ~ abund + geo_dist_std + order")
mod2 <- brm_logistic("ext_signal ~ abund + range_lat_std + order")

# number of genera
mod3 <- brm_logistic("ext_signal ~ n_genus + geo_dist_std + order")
mod4 <- brm_logistic("ext_signal ~ n_genus + range_lat_std + order")


# second adjustment subset
mod5 <- brm_logistic("ext_signal ~ n_genus + mean_q_std + cont_area + order +
                     n_units_std + temp_deep_st:temp_deep_lt3 + sr_value_std +
                     shark_collections_std + temp_gat_binned")

mod6 <- brm_logistic("ext_signal ~ n_genus + mean_q_std + cont_area + order +
                     outcrop_area_std + temp_deep_st:temp_deep_lt3 + d13C_std +
                     shark_collections_std + temp_gat_binned")

mod7 <- brm_logistic("ext_signal ~ abund + mean_q_std + cont_area + order +
                     n_units_std + temp_deep_st:temp_deep_st:temp_deep_lt2 + sr_value_std +
                     shark_collections_std + temp_deep_binned")

mod8 <- brm_logistic("ext_signal ~ abund + mean_q_std + cont_area + order +
                     outcrop_area_std + temp_deep_st:temp_deep_lt3 + sr_value_std +
                     shark_collections_std + temp_gat_binned")

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2,
                                   mod3, mod4,
                                   mod5, mod6,
                                   mod7, mod8,
                                   variable = c("b_abund", "b_n_genus"),
                                   seed = 1708,
                                   ndraws =  1e4, 
                                   missing = NA) %>% 
  pivot_longer(cols = everything(), 
               names_to = "coef_name", 
               values_to = "coef_val") %>% 
  drop_na(coef_val)

# save trend predictions
dat_pred_post %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_abund_genus.rds"), 
            compress = "gz")



# geographic range --------------------------------------------------------

# first adjustment set
mod1 <- brm_logistic("ext_signal ~ geo_dist_std + abund + order")
mod2 <- brm_logistic("ext_signal ~ range_lat_std + abund + order")

# number of genera
mod3 <- brm_logistic("ext_signal ~ geo_dist_std + n_genus + order")
mod4 <- brm_logistic("ext_signal ~ range_lat_std + n_genus + order")


# second adjustment subset
mod5 <- brm_logistic("ext_signal ~ range_lat_std + n_units_std  + temp_deep_st:temp_deep_lt3 +
                     mean_q_std + sr_value_std + shark_collections_std  + cont_area +
                     order + temp_gat_binned")

mod6 <- brm_logistic("ext_signal ~ range_lat_std + outcrop_area_std  + temp_deep_st:temp_deep_lt3 +
                     mean_q_std + d13C_std + shark_collections_std  + cont_area +
                     order + temp_gat_binned")

mod7 <- brm_logistic("ext_signal ~ geo_dist_std  + n_units_std  + temp_deep_st:temp_deep_lt2 +
                     mean_q_std + sr_value_std + shark_collections_std  + cont_area +
                     order + temp_deep_binned")

mod8 <- brm_logistic("ext_signal ~ geo_dist_std + outcrop_area_std  + temp_deep_st:temp_deep_lt3 +
                     mean_q_std + sr_value_std + shark_collections_std  + cont_area +
                     order + temp_gat_binned")

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2,
                                   mod3, mod4,
                                   mod5, mod6,
                                   mod7, mod8,
                                   variable = c("b_geo_dist_std", "b_range_lat_std"),
                                   seed = 1708,
                                   ndraws =  1e4, 
                                   missing = NA) %>% 
  pivot_longer(cols = everything(), 
               names_to = "coef_name", 
               values_to = "coef_val") %>% 
  drop_na(coef_val)

# save trend predictions
dat_pred_post %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_geo_range_genus.rds"), 
            compress = "gz")



# paleotemperature --------------------------------------------------------

# start with deep ocean temperature
# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_deep_st:temp_deep_lt1")
mod2 <- brm_logistic("ext_signal ~ temp_deep_st:temp_deep_lt2")
mod3 <- brm_logistic("ext_signal ~ temp_deep_st:temp_deep_lt3")
mod4 <- brm_logistic("ext_signal ~ temp_deep_st:temp_deep_lt4")

# same for global average temperature
mod5 <- brm_logistic("ext_signal ~ temp_gat_st:temp_gat_lt1")
mod6 <- brm_logistic("ext_signal ~ temp_gat_st:temp_gat_lt2")
mod7 <- brm_logistic("ext_signal ~ temp_gat_st:temp_gat_lt3")
mod8 <- brm_logistic("ext_signal ~ temp_gat_st:temp_gat_lt4")

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2,
                                   mod3, mod4,
                                   mod5, mod6,
                                   mod7, mod8,
                                   seed = 1708,
                                   ndraws =  1e4, 
                                   missing = NA) %>% 
  select(contains("temp")) %>% 
  pivot_longer(cols = everything(), 
               names_to = "coef_name", 
               values_to = "coef_val") %>% 
  drop_na(coef_val)


# save trend predictions
dat_pred_post %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_paleotemperature_genus.rds"), 
            compress = "gz")


# productivity ------------------------------------------------------------

## first adjustment set 
# for d13C
mod1 <- brm_logistic("ext_signal ~ d13C_std + 
                     latitude_pref_abs + cont_area +
                     shark_collections_std + temp_gat_binned")

mod2 <- brm_logistic("ext_signal ~ d13C_std + 
                     latitude_pref_abs + cont_area + 
                     pbdb_collections_std + temp_gat_binned")

mod3 <- brm_logistic("ext_signal ~ d13C_std + 
                     latitude_pref_abs + cont_area + 
                     pbdb_collections_std + temp_deep_binned")

mod4 <- brm_logistic("ext_signal ~ d13C_std + 
                     latitude_pref_abs + cont_area + 
                     shark_collections_std + temp_deep_binned")

# for Sr ratio 
mod5 <- brm_logistic("ext_signal ~ sr_value_std + 
                     latitude_pref_abs + cont_area + 
                     shark_collections_std + temp_gat_binned")

mod6 <- brm_logistic("ext_signal ~ sr_value_std + 
                     latitude_pref_abs + cont_area + 
                     pbdb_collections_std + temp_gat_binned")

mod7 <- brm_logistic("ext_signal ~ sr_value_std + 
                     latitude_pref_abs + cont_area + 
                     pbdb_collections_std + temp_deep_binned")

mod8 <- brm_logistic("ext_signal ~ sr_value_std + 
                     latitude_pref_abs + cont_area + 
                     shark_collections_std + temp_deep_binned")


## second adjustment set 
# for d13C
mod9 <- brm_logistic("ext_signal ~ d13C_std + latitude_pref_abs + 
                     sea_level + cont_area + 
                     shark_collections_std + temp_gat_binned")

mod10 <- brm_logistic("ext_signal ~ d13C_std + latitude_pref_abs +
                      sea_level + cont_area +
                      pbdb_collections_std + temp_gat_binned")

mod11 <- brm_logistic("ext_signal ~ d13C_std + latitude_pref_abs + 
                      sea_level + cont_area +
                      pbdb_collections_std + temp_deep_binned")

mod12 <- brm_logistic("ext_signal ~ d13C_std + latitude_pref_abs + 
                      sea_level + cont_area +
                      shark_collections_std + temp_deep_binned")

# for Sr ratio 
mod13 <- brm_logistic("ext_signal ~ sr_value_std + latitude_pref_abs + 
                     sea_level + cont_area + 
                     shark_collections_std + temp_gat_binned")

mod14 <- brm_logistic("ext_signal ~ sr_value_std + latitude_pref_abs +
                      sea_level + cont_area +
                      pbdb_collections_std + temp_gat_binned")

mod15 <- brm_logistic("ext_signal ~ sr_value_std + latitude_pref_abs + 
                      sea_level + cont_area +
                      pbdb_collections_std + temp_deep_binned")

mod16 <- brm_logistic("ext_signal ~ sr_value_std + latitude_pref_abs + 
                      sea_level + cont_area +
                      shark_collections_std + temp_deep_binned")

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2, mod3, mod4, 
                                   mod5, mod6, mod7, mod8,
                                   mod9, mod10, mod11, mod12,
                                   mod13, mod14, mod15, mod16,
                                   variable = c("b_d13C_std", "b_sr_value_std"),
                                   seed = 1708,
                                   ndraws =  1e4, 
                                   missing = NA) %>% 
  pivot_longer(cols = everything(), 
               names_to = "coef_name", 
               values_to = "coef_val") %>% 
  drop_na(coef_val)


# save trend predictions
dat_pred_post %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_productivity_genus.rds"))


# sea level ---------------------------------------------------------------

# average over temperature estimates 
mod1 <- brm_logistic("ext_signal ~ sea_level + temp_gat_binned")
mod2 <- brm_logistic("ext_signal ~ sea_level + temp_deep_binned")

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2,
                                   variable = "b_sea_level",
                                   seed = 1708,
                                   ndraws =  1e4, 
                                   missing = NA) %>% 
  pivot_longer(cols = everything(), 
               names_to = "coef_name", 
               values_to = "coef_val") %>% 
  drop_na(coef_val)


# save trend predictions
dat_pred_post %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_sea_level_genus.rds"), 
            compress = "gz")



# shelf area --------------------------------------------------------------

# one model with sea level 
mod1 <- brm_logistic("ext_signal ~ cont_area + sea_level")


# average posterior draws by model stacking
dat_pred_post <- as_draws_df(mod1, 
                             variable = "b_cont_area") %>% 
  as_tibble() %>% 
  select(coef_val = b_cont_area) %>% 
  slice_sample(n = 1e4)


# save trend predictions
dat_pred_post %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_shelf_area_genus.rds"), 
            compress = "gz")




# temperature -------------------------------------------------------------

# start with deep ocean temperature
# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt1")
mod2 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt2")
mod3 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt3")
mod4 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt4")

# same for global average temperature
mod5 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt1")
mod6 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt2")
mod7 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt3")
mod8 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt4")


# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2,
                                   mod3, mod4,
                                   mod5, mod6,
                                   mod7, mod8,
                                   variable = c("b_temp_gat_binned", "b_temp_deep_binned"),
                                   seed = 1708,
                                   ndraws =  1e4, 
                                   missing = NA) %>% 
  pivot_longer(cols = everything(), 
               names_to = "coef_name", 
               values_to = "coef_val") %>% 
  drop_na(coef_val)

# save predictions
dat_pred_post %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_temperature_genus.rds"), 
            compress = "gz")

