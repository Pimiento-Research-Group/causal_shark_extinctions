# load packages
library(tidyverse)
library(here)
library(readxl)
library(brms)


# plotting configurations
source(here("R", "config_file.R"))


# bin proxy data ----------------------------------------------------------

### sea level
dat_sealevel <- read_xlsx(here("data",
                                    "raw",
                                    "sea_level",
                                    "miller_et_al_2005.xlsx")) %>% 
  # clean up column names
  rename(sea_level = "Sea level best estimate", 
         age = "Age for best estimate (MA)")  %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(sea_level_mean = mean(sea_level)) 





### productivity
#d13C
dat_13C <- read_csv(here("data",
                              "raw",
                              "productivity",
                              "veizer_prokoph_2015.csv")) %>% 
  # clean up column names
  select(age = "gts2012", 
         d13C)  %>% 
  filter(age <= 180) %>% 
  drop_na(age, d13C) %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(d13C_mean = mean(d13C)) 

#Sr87/86
dat_SR <- read_csv(here("data",
                             "raw",
                             "productivity",
                             "mcarthur_howarth_shields_2012.csv"), 
                        col_names = FALSE) %>% 
  # clean up column names
  select(age = X1, 
         sr_value = X2) %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(sr_value_mean = mean(sr_value)) %>% 
  add_row(tibble(bin = 15, 
                 sr_value_mean = mean(read_csv(here("data",
                                                    "raw",
                                                    "productivity",
                                                    "mcarthur_howarth_shields_2012.csv"),
                                               col_names = FALSE)$X2)))



### outcrop area

# rock units
dat_marine_units <- read_csv(here("data",
                                       "raw",
                                       "outcrop_area",
                                       "macrostrat_24_11_2022.csv")) %>% 
  # select only marine environments
  filter(!str_detect(.$environ, "non-marine")) %>% 
  # bin to 10myr
  mutate(age = (t_age+b_age)/2, 
         bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  select(bin, unit_id, age, environ) %>% 
  count(bin) %>% 
  # add bins with zero counts
  rename(n_units = n)


# outcrop area
dat_outcrop <- read_xlsx(here("data",
                              "raw",
                              "outcrop_area",
                              "wall_ivany_wilkinson_2009.xlsx")) %>% 
  # clean up colnames
  rename(age = "age (ma)", 
         outcrop_area = "cumul_area (10^6 km^2)") %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  # use natural spline for interpolation
  { spline(.$bin, .$outcrop_area, 
           xout = 1:15) } %>% 
  as_tibble() %>% 
  rename(bin = x, outcrop_area = y) 



### shelf area
dat_cont_area <- read_csv(here("data",
                                    "raw",
                                    "shelf_area",
                                    "kocsis_scotese_2021.csv")) %>% 
  # clean up colnames
  rename(cont_area = "shelf-rgeos") %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(cont_area = mean(cont_area))


### temperature
dat_temp <- read_xlsx(here("data",
                                "raw",
                                "temperature",
                                "scotese_et_al_2021.xlsx")) %>% 
  # clean up colnames
  select(age = Age,
         temp_gat = GAT,
         temp_deep = "Deep Ocean") %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(250, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(temp_gat_binned = mean(temp_gat), 
            temp_deep_binned = mean(temp_deep))

### paleotemperature
# calculate lagged temperatures trends
dat_paleotemp <- dat_temp %>% 
  mutate(temp_gat_st = temp_gat_binned - lead(temp_gat_binned), 
         temp_gat_lt1 = lead(temp_gat_st), 
         temp_gat_lt2 = lead(temp_gat_st, n = 2), 
         temp_gat_lt3 = lead(temp_gat_st, n = 3), 
         temp_gat_lt4 = lead(temp_gat_st, n = 4), 
         # same for deep ocean temperature
         temp_deep_st = temp_deep_binned - lead(temp_deep_binned), 
         temp_deep_lt1 = lead(temp_deep_st), 
         temp_deep_lt2 = lead(temp_deep_st, n = 2), 
         temp_deep_lt3 = lead(temp_deep_st, n = 3), 
         temp_deep_lt4 = lead(temp_deep_st, n = 4)) %>% 
  # add missing bin
  fill(contains("temp"), .direction = "downup") %>% 
  filter(between(bin, 1, 15)) %>% 
  select(-contains("binned"))


# fossil data -------------------------------------------------------------


# get those species that survive until now
spec_modern <- read_delim(here("data",
                               "fossil_occurrences",
                               "combined_10_se_est_species_names.txt")) %>% 
  filter(te == 0) %>% 
  pull(species)

# load PyRate estimates
dat_ext <- read_delim(here("data",
                               "fossil_occurrences",
                               "combined_10_se_est_species_names.txt")) %>% 
  # estimate fad and lad for each species
  mutate(across(c(ts, te), abs)) %>% 
  group_by(species) %>% 
  summarise(ori_age = mean(ts),
            ext_age = mean(te)) %>% 
  # bin to 10myr
  mutate(bin_ori = cut(ori_age, breaks = seq(250, 0, by = -10),
                       include.lowest = TRUE,
                       labels = FALSE), 
         bin_ext = cut(ext_age, breaks = seq(250, 0, by = -10),
                       include.lowest = TRUE,
                       labels = FALSE)) %>% 
  drop_na(bin_ori, bin_ext) %>% 
  # fill in duration bins
  mutate(bin_occ = map2(.x = bin_ori, 
                        .y = bin_ext,
                        .f = ~ seq(.x, .y, by = -1))) %>% 
  select(species, bin_occ, bin_ext) %>% 
  unnest(bin_occ) %>% 
  # create extinction signal
  group_by(species) %>% 
  mutate(ext_signal = if_else(bin_occ == bin_ext, 1, 0), 
         # assign 0 to those species that survive until the modern
         ext_signal = if_else(species %in% spec_modern, 
                              0, ext_signal)) %>% 
  # clean up
  ungroup() %>% 
  select(-bin_ext, 
         bin = bin_occ) %>% 
  mutate(accepted_name = str_replace(species, "_", " ")) %>%
  select(-species)



# load occurrence database
dat_occurrences <- read_rds(here("data",
                                 "fossil_occurrences",
                                 "database_occurrences_15_Apr_2023.rds"))

# bin the occurrences to 10 myr
dat_occ_binned <- dat_occurrences %>% 
  mutate(bin_min = cut(Min_Ma, breaks = seq(250, 0, by = -10),
                       include.lowest = TRUE,
                       labels = FALSE), 
         bin_max = cut(Max_Ma, breaks = seq(250, 0, by = -10),
                       include.lowest = TRUE,
                       labels = FALSE)) %>% 
  # select entries, where the early and late interval fields indicate 
  # the same stg or the late_interval field is empty
  filter(bin_min == bin_max | is.na(bin_max)) %>% 
  # in these entries, use the stg indicated by the early_interval
  select(bin = bin_min, everything(), -bin_max) %>% 
  # clean up
  drop_na(bin) %>% 
  # get the age estimate for the corresponding bin
  left_join(tibble(bin = 25:1, age = seq(245, 5, by = -10)))




# number of occurrences
dat_abund <- dat_occ_binned %>% 
  # select only species
  filter(rank == "species") %>% 
  count(bin, accepted_name) 

# number of species within genera
dat_abund_genus <- dat_occ_binned %>% 
  # select only species
  filter(rank == "species") %>% 
  count(bin, genus) 


# estimate sampling effort by the number of collections per bin
dat_sampling <- dat_occ_binned %>% 
  distinct(bin, collection_no) %>% 
  count(bin, 
        name = "shark_collections") %>% 
  # add bins with zero counts
  complete(bin = 1:14, fill = list(shark_collections = 0))

# alternatively the number of collections per bin from the pbdb

# download pbdb collections from bin 72 (Hauterivian) to bin 95 (Holocene)
# via the api call
dat_pbdb <- read_csv("https://paleobiodb.org/data1.2/colls/list.csv?interval=Hauterivian,Holocene")

# bin the data
dat_pbdb_sampling <- dat_pbdb %>% 
  mutate(bin_min = cut(min_ma, breaks = seq(250, 0, by = -10),
                       include.lowest = TRUE,
                       labels = FALSE), 
         bin_max = cut(max_ma, breaks = seq(250, 0, by = -10),
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
  count(bin, 
        name = "pbdb_collections") %>% 
  # add bins with zero counts
  complete(bin = 1:14, fill = list(pbdb_collections = 0))



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
  mutate(bin = cut(age, breaks = seq(250, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)-1, 
         .before = 1) %>% 
  # summarise
  drop_na(bin) %>% 
  group_by(bin, species) %>% 
  summarise(mean_q = mean(mean_q)) %>% 
  ungroup() %>% 
  # clean up species names for joining
  mutate(species = str_remove(species, "_TS"), 
         accepted_name = str_replace(species, "_", " ")) %>% 
  select(bin, accepted_name, mean_q) %>% 
  arrange(accepted_name)


# combine all datasets to one
dat_full <- list(dat_ext,
                 dat_preservation,
                 dat_sealevel,
                 dat_13C,
                 dat_SR,
                 dat_marine_units,
                 dat_outcrop,
                 dat_cont_area,
                 dat_temp,
                 dat_paleotemp,
                 dat_sampling,
                 dat_pbdb_sampling) %>% 
  map(~ filter(.x, bin %in% 1:14)) %>% 
  reduce(full_join) %>%  
  # remove missing values
  left_join(dat_preservation %>%
              group_by(accepted_name) %>%
              summarise(mean_q_av = mean(mean_q))) %>% 
  mutate(mean_q = if_else(is.na(mean_q), mean_q_av, mean_q)) %>% 
  drop_na(everything()) %>% 
  # get taxonomy
  left_join(dat_occurrences %>% 
              distinct(accepted_name, genus, family, order)) %>% 
  select(order, family, genus, species = accepted_name,
         everything())  



# scale predictors --------------------------------------------------------


# bring predictors on meaningful scales to enable comparisons
dat_merged <- dat_full %>% 
  mutate(across(c(d13C_mean, sr_value_mean, # scale all productivity parameters
                  n_units, outcrop_area, # same for outcrop parameters
                  mean_q, # preservation rate
                  pbdb_collections, shark_collections), #  outcrop parameters
                ~ scale(.x)[,1], 
                .names = "{col}_std")) %>% 
  select(bin, order, family, genus, species,
         ext_signal, 
         sea_level_mean, 
         cont_area, 
         contains("temp"), 
         contains("std")) %>% 
  drop_na()

# save dataset
dat_merged %>% 
  write_rds(here("data",
                 "processed_fossil_data_10myr.rds"))






# productivity ------------------------------------------------------------

## first adjustment set 
# for d13C
mod1 <- brm_logistic("ext_signal ~ d13C_mean_std + 
                     pbdb_collections_std + sea_level_mean + 
                     cont_area + temp_gat_binned")

mod2 <- brm_logistic("ext_signal ~ d13C_mean_std +
                     outcrop_area_std + sea_level_mean + 
                     cont_area + temp_gat_binned")

mod3 <- brm_logistic("ext_signal ~ d13C_mean_std + 
                     pbdb_collections_std + sea_level_mean +
                     cont_area + temp_deep_binned")

mod4 <- brm_logistic("ext_signal ~ d13C_mean_std + 
                     outcrop_area_std + sea_level_mean + 
                     cont_area + temp_deep_binned")

# for Sr ratio 
mod5 <- brm_logistic("ext_signal ~ sr_value_mean_std + 
                     pbdb_collections_std + sea_level_mean + 
                     cont_area + temp_gat_binned")

mod6 <- brm_logistic("ext_signal ~ sr_value_mean_std + 
                     outcrop_area_std + sea_level_mean + 
                     cont_area + temp_gat_binned")

mod7 <- brm_logistic("ext_signal ~ sr_value_mean_std + 
                     pbdb_collections_std + sea_level_mean +
                     cont_area + temp_deep_binned")

mod8 <- brm_logistic("ext_signal ~ sr_value_mean_std + 
                     outcrop_area_std + sea_level_mean + 
                     cont_area + temp_deep_binned")


## second adjustment set 
# outcrop area, preservation potential, sampling effort, shelf area, 
# taxonomic identity, temperature
# for d13C
mod9 <- brm_logistic("ext_signal ~ d13C_mean_std + 
                     outcrop_area_std + mean_q_std + 
                     pbdb_collections_std + cont_area + 
                     order + temp_gat_binned")

mod10 <- brm_logistic("ext_signal ~ d13C_mean_std + 
                      outcrop_area_std + mean_q_std +
                      shark_collections_std + cont_area +
                      order + temp_gat_binned")

mod11 <- brm_logistic("ext_signal ~ d13C_mean_std + 
                      outcrop_area_std + mean_q_std +
                      pbdb_collections_std + cont_area +
                      order + temp_deep_binned")

mod12 <- brm_logistic("ext_signal ~ d13C_mean_std + 
                      outcrop_area_std + mean_q_std +
                      shark_collections_std + cont_area +
                      order + temp_deep_binned")

# for Sr ratio 
mod13 <- brm_logistic("ext_signal ~ sr_value_mean_std + 
                      outcrop_area_std + mean_q_std + 
                      pbdb_collections_std + cont_area + 
                      order + temp_gat_binned")

mod14 <- brm_logistic("ext_signal ~ sr_value_mean_std + 
                      outcrop_area_std + mean_q_std +
                      shark_collections_std + cont_area +
                      order + temp_gat_binned")

mod15 <- brm_logistic("ext_signal ~ sr_value_mean_std + 
                      outcrop_area_std + mean_q_std +
                      pbdb_collections_std + cont_area +
                      order + temp_deep_binned")

mod16 <- brm_logistic("ext_signal ~ sr_value_mean_std + 
                      outcrop_area_std + mean_q_std +
                      shark_collections_std + cont_area +
                      order + temp_deep_binned")




# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2, mod3, mod4, 
                                   mod5, mod6, mod7, mod8,
                                   mod9, mod10, mod11, mod12,
                                   mod13, mod14, mod15, mod16,
                                   variable = c("b_d13C_mean_std", 
                                                "b_sr_value_mean_std"),
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
                 "pred_trend_productivity_10myr.rds"))


# sea level ---------------------------------------------------------------

# average over temperature estimates 
mod1 <- brm_logistic("ext_signal ~ sea_level_mean + temp_gat_binned")
mod2 <- brm_logistic("ext_signal ~ sea_level_mean + temp_deep_binned")

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2,
                                   variable = "b_sea_level_mean",
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
                 "pred_trend_sea_level_10myr.rds"), 
            compress = "gz")



# shelf area --------------------------------------------------------------

# one model with sea level 
mod1 <- brm_logistic("ext_signal ~ cont_area + sea_level_mean")


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
                 "pred_trend_shelf_area_10myr.rds"), 
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
                 "pred_trend_temperature_10myr.rds"), 
            compress = "gz")


# estimate trend lines
# set up grid to average over
dat_new <- tibble(temp_deep_binned = 0:25) %>%
  expand_grid(temp_deep_st = -3:3, 
              temp_deep_lt = -3:3) %>% 
  mutate(temp_deep_lt1 = temp_deep_lt,
         temp_deep_lt2 = temp_deep_lt,
         temp_deep_lt3 = temp_deep_lt,
         temp_deep_lt4 = temp_deep_lt,
         temp_gat_binned = temp_deep_binned, 
         temp_gat_st = temp_deep_st,
         temp_gat_lt1 = temp_deep_lt,
         temp_gat_lt2 = temp_deep_lt,
         temp_gat_lt3 = temp_deep_lt,
         temp_gat_lt4 = temp_deep_lt)

# set up number of draws 
nr_draws <- 100

# average prediction by model stacking
dat_pred <- pp_average(mod1, mod2,
                       mod3, mod4, 
                       mod5, mod6,
                       mod7, mod8,
                       newdata = dat_new,
                       seed = 1708,
                       summary = FALSE, 
                       method = "posterior_epred", 
                       ndraws = nr_draws) %>% 
  as_tibble() %>% 
  mutate(nr_draw = rownames(.)) %>% 
  pivot_longer(cols = contains("V")) %>% 
  add_column(temperature = rep(dat_new$temp_deep_binned, nr_draws)) %>% 
  group_by(nr_draw, temperature) %>%
  ggdist::mean_qi(value) %>% 
  select(temperature, value, nr_draw)


# save predictions
dat_pred %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_temperature_tenmyr.rds"), 
            compress = "gz")
