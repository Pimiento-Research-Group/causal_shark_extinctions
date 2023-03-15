# load packages
library(tidyverse)
library(here)
library(readxl)
library(brms)
library(tidybayes)

# plotting configurations
source(here("R", "config_file.R"))



# sea level ---------------------------------------------------------------

### cenozoic sea level from miller et al 2020 ###
dat_sealevel <- read_xlsx(here("data",
                               "raw",
                               "sea_level",
                               "miller_et_al_2020.xlsx")) %>% 
  # clean up column names
  rename(sea_level = "Sea Level (m) Smoothed", 
         age = "Age (ka)")  %>% 
  # transform to million years
  mutate(age = age/1000) %>% 
  # bin to 1myr
  mutate(bin = cut(age, breaks = seq(66, 0, by = -1),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(sea_level_mean = mean(sea_level)) %>% 
  complete(bin = 66:1) %>%
  fill(sea_level_mean, .direction = "downup")


# productivity ------------------------------------------------------------


### Cenozoic delta13C values from Westerhold et al 2020 ###
dat_13C <- read_xlsx(here("data",
                          "raw",
                          "productivity",
                          "westerhold_et_al_2020.xlsx"),
                     sheet = "Table S34",
                     skip = 1) %>% 
  # remove redundant columns
  select(age = age_tuned, 
         d13C = ISOBENd13cLOESSsmooth) %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(66, 0, by = -1),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(d13C_mean = mean(d13C))



### Cenozoic diatom diversity from the Neptune database ###
# load data
dat_diatom <- read_tsv(here("data",
                            "raw",
                            "productivity",
                            "diatom_2022-11-27_09-58-09.csv")) %>% 
  # clean up column names
  select(age = "Age (Ma) Gradstein et al. 2012", 
         taxon_id = "Resolved Taxon ID") 

# rarefaction based on this paper: https://www.nature.com/articles/nature07435
# number of rarefaction iteration
nr_iter <- 1000
# number of samples
nr_samples <- 96
# pre-allocate list
dat_diatom_list <- vector(mode = "list", 
                          length = nr_iter)
# iterate
for (i in 1:nr_iter) {
  
  dat_diatom_list[[i]] <- dat_diatom_ceno %>% 
    # bin data
    mutate(bin = cut(age, breaks = seq(90, 0, by = -1))) %>% 
    group_by(bin) %>% 
    # rarefaction
    slice_sample(n = 96, replace = TRUE) %>% 
    # calculate diversity
    summarise(count = n_distinct(taxon_id)) %>% 
    # add bins with zero counts
    complete(bin, fill = list(count = 0)) %>% 
    ungroup() %>% 
    # add numeric age for plotting
    add_column(age = seq(0.5, 89.5), 
               rare_id = i)
  
}

# summarize
dat_diatom <- dat_diatom_list %>% 
  bind_rows() %>% 
  group_by(bin, age) %>% 
  summarise(div_mean = mean(count), 
            div_sd = sd(count)) %>% 
  mutate(div_low = div_mean - 1.96*div_sd, 
         div_high = div_mean + 1.96*div_sd) %>% 
  # bin to 1myr
  mutate(bin = cut(age, breaks = seq(66, 0, by = -1),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(diatom_rich = mean(div_mean)) 


# # save file
# dat_diatom %>% 
#   write_rds(here("data", 
#                  "raw", 
#                  "productivity", 
#                  "diatom_full.rds"))


# outcrop area ------------------------------------------------------------


### Phanerozoic rock units from the macrostrat database ###
dat_marine_units <- read_csv(here("data",
                                  "raw",
                                  "outcrop_area",
                                  "macrostrat_24_11_2022.csv")) %>% 
  # select only marine environments
  filter(!str_detect(.$environ, "non-marine")) %>% 
  # bin to 1myr
  mutate(age = (t_age+b_age)/2, 
         bin = cut(age, breaks = seq(66, 0, by = -1),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  select(bin, unit_id, age, environ) %>% 
  count(bin, 
        name = "n_units") %>% 
  arrange(bin) %>% 
  filter(bin <= 66)



# shelf area --------------------------------------------------------------

### Cenozoic flooded continental area (10^6 km^2) from Miller et al 2005  ###
dat_cont_area <- read_csv(here("data",
                               "raw",
                               "shelf_area",
                               "miller_et_al_2005.csv"),
                          col_names = FALSE) %>% 
  # clean up colnames
  rename(age = X2, 
         area = X1) %>% 
  # use natural spline for interpolation and binning
  { spline(.$age, .$area, 
           xout = 66:1) } %>% 
  as_tibble() %>% 
  rename(bin = x, cont_area = y) %>% 
  arrange(bin)


# temperature -------------------------------------------------------------

### Cenozoic deep-ocean temperature  used with equation 7afrom Cramer et al 2011 ###
dat_temp <- read_xlsx(here("data",
                           "raw",
                           "temperature",
                           "cramer_et_al_2011_7a.xlsx")) %>% 
  # clean up colnames
  select(age = Age,
         temp = Temperature,
         temp_low = "Temperature min",
         temp_high = "Temperature max") %>% 
  # bin to 1myr
  mutate(bin = cut(age, breaks = seq(66, 0, by = -1),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(temp_binned = mean(temp)) %>%
  complete(bin = 66:1) %>%
  fill(temp_binned, .direction = "downup")


# paleotemperature --------------------------------------------------------

# calculate lagged temperatures trends
dat_paleotemp <- dat_temp %>% 
  mutate(temp_st = temp_binned - lead(temp_binned), 
         temp_lt1 = lead(temp_st), 
         temp_lt2 = lead(temp_st, n = 2), 
         temp_lt3 = lead(temp_st, n = 3), 
         temp_lt4 = lead(temp_st, n = 4)) %>% 
  select(-contains("binned")) %>% 
  # add missing bin
  complete(bin = 66:1) %>% 
  fill(contains("temp"), .direction = "downup")


# extinction signal -------------------------------------------------------

# load PyRate estimates
dat_ext <- read_delim(here("data",
                           "fossil_occurrences",
                           "all_species_10_Grj_se_est_species_names.txt")) %>% 
  # estimate fad and lad for each species
  pivot_longer(cols = contains("ts"), 
               names_to = "origination", 
               values_to = "origination_age") %>% 
  pivot_longer(cols = contains("te"), 
               names_to = "extinction", 
               values_to = "extinction_age") %>% 
  group_by(species) %>% 
  summarise(ori_age = mean(origination_age),
            ext_age = mean(extinction_age)) %>% 
  # bin fad and lad to 1myr
  mutate(bin_ori = cut(ori_age, breaks = seq(66, 0, by = -1),
                       include.lowest = TRUE,
                       labels = FALSE), 
         bin_ext = cut(ext_age, breaks = seq(66, 0, by = -1),
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
  mutate(ext_signal = if_else(bin_occ == bin_ext, 1, 0))  %>% 
  # clean up
  ungroup() %>% 
  select(-bin_ext, 
         bin = bin_occ) %>% 
  mutate(accepted_name = str_replace(species, "_", " ")) %>%
  select(-species) %>% 
  filter(between(bin, 1, 66)) 




# sampling effort ---------------------------------------------------------

# load occurrence database
dat_occurrences <- read_rds(here("data",
                                 "fossil_occurrences",
                                 "database_occurrences_10_Jan_2023.rds"))

# bin the occurrences to 1myr
dat_occ_binned <- dat_occurrences %>% 
  mutate(bin_min = cut(Min_Ma, breaks = seq(66, 0, by = -1),
                       include.lowest = TRUE,
                       labels = FALSE), 
         bin_max = cut(Max_Ma, breaks = seq(66, 0, by = -1),
                       include.lowest = TRUE,
                       labels = FALSE)) %>% 
  # select entries, where the early and late interval fields indicate 
  # the same stg or the late_intervar field is empty
  filter(bin_min == bin_max | is.na(bin_max)) %>% 
  # in these entries, use the stg indicated by the early_interval
  select(bin = bin_min, everything(), -bin_max) %>% 
  # clean up
  drop_na(bin) 

# estimate sampling effort by the number of collections per bin
dat_sampling <- dat_occ_binned %>% 
  distinct(bin, collection_no) %>% 
  count(bin, 
        name = "shark_collections") %>% 
  # add bins with zero counts
  complete(bin = 66:1, fill = list(shark_collections = 0)) %>% 
  filter(between(bin, 1, 66))

# alternatively the number of collections per bin from the pbdb

# download pbdb collections from bin 72 (Hauterivian) to bin 95 (Holocene)
# via the api call
dat_pbdb <- read_csv("https://paleobiodb.org/data1.2/colls/list.csv?interval=Hauterivian,Holocene")

# bin the data
dat_pbdb_sampling <- dat_pbdb %>% 
  mutate(bin_min = cut(min_ma, breaks = seq(66, 0, by = -1),
                       include.lowest = TRUE,
                       labels = FALSE), 
         bin_max = cut(max_ma, breaks = seq(66, 0, by = -1),
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
  complete(bin = 66:1, fill = list(pbdb_collections = 0)) %>% 
  filter(between(bin, 1, 66))


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
  # bin to stages
  mutate(bin = cut(age, breaks = seq(66, 0, by = -1),
                   include.lowest = TRUE,
                   labels = FALSE)+1, 
         .before = 1) %>% 
  # summarise
  drop_na(bin) %>% 
  group_by(bin, species) %>% 
  summarise(mean_q = mean(mean_q)) %>% 
  ungroup() %>% 
  filter(between(bin, 1, 66)) %>% 
  # clean up species names for joining
  mutate(species = str_remove(species, "_TS"), 
         accepted_name = str_replace(species, "_", " ")) %>% 
  select(bin, accepted_name, mean_q) %>% 
  arrange(accepted_name) 


# merge and combine -------------------------------------------------------

# combine all datasets to one
dat_full <- list(dat_ext,
                 dat_preservation,
                 dat_sealevel,
                 dat_13C,
                 dat_diatom,
                 dat_marine_units,
                 dat_cont_area,
                 dat_temp,
                 dat_paleotemp,
                 dat_sampling,
                 dat_pbdb_sampling) %>% 
  reduce(full_join)  %>% 
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
         everything()) %>% 
  drop_na()


# scale predictors --------------------------------------------------------


# bring predictors on meaningful scales to enable comparisons
dat_merged <- dat_full %>% 
  rename(d13C = d13C_mean, 
         sea_level = sea_level_mean) %>% 
  mutate(across(c(d13C, diatom_rich ,# scale all productivity parameters
                  n_units, # same for outcrop parameters
                  mean_q, # preservation rate
                  pbdb_collections, shark_collections), # sampling effort
                ~ scale(.x)[,1], 
                .names = "{col}_std")) %>% 
  select(bin, order, family, genus, species,
         ext_signal, 
         sea_level, 
         diatom_rich,
         cont_area, 
         contains("temp"), 
         contains("std")) 

# save dataset
dat_merged %>% 
  write_rds(here("data",
                 "processed_fossil_data_cenozoic.rds"))


# mod temperature ---------------------------------------------------------

# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt1")
mod2 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt2")
mod3 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt3")
mod4 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt4")

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2,
                                   mod3, mod4,
                                   variable = "b_temp_binned",
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
                 "pred_trend_temperature_cenozoic.rds"), 
            compress = "gz")



# mod shelf area ----------------------------------------------------------

# one model with sea level 
mod1 <- brm_logistic("ext_signal ~ cont_area + sea_level")


# average posterior draws by model stacking
dat_pred_post <- as_draws_df(mod1,
                             variable = "b_cont_area") %>% 
  as_tibble() %>% 
  select(coef_val = b_cont_area) %>% 
  slice_sample(n = 1e4)


# save predictions
dat_pred_post %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_shelf_area_cenozoic.rds"), 
            compress = "gz")



# mod sea level -----------------------------------------------------------

# average over temperature estimates 
mod1 <- brm_logistic("ext_signal ~ sea_level + temp_binned")

# average posterior draws by model stacking
dat_pred_post <- as_draws_df(mod1,
                             variable = "b_sea_level") %>% 
  as_tibble() %>% 
  select(coef_val = b_sea_level) %>% 
  slice_sample(n = 1e4)


# save predictions
dat_pred_post %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_sea_level_cenozoic.rds"), 
            compress = "gz")



# mod productivity ------------------------------------------------------------

## first adjustment set 
# for d13C
mod1 <- brm_logistic("ext_signal ~ d13C_std + 
                     pbdb_collections_std + sea_level + 
                     cont_area + temp_binned")

mod2 <- brm_logistic("ext_signal ~ d13C_std +
                     n_units_std + sea_level + 
                     cont_area + temp_binned")


# for diatom richness 
mod3 <- brm_logistic("ext_signal ~ diatom_rich_std + 
                     pbdb_collections_std + sea_level + 
                     cont_area + temp_binned")

mod4 <- brm_logistic("ext_signal ~ diatom_rich_std + 
                     n_units_std + sea_level + 
                     cont_area + temp_binned")


## second adjustment set 
# outcrop area, preservation potential, sampling effort, shelf area, 
# taxonomic identity, temperature
# for d13C
mod5 <- brm_logistic("ext_signal ~ d13C_std + 
                     n_units_std + mean_q_std + 
                     pbdb_collections_std + cont_area + 
                     order + temp_binned")

mod6 <- brm_logistic("ext_signal ~ d13C_std + 
                      n_units_std + mean_q_std +
                      shark_collections_std + cont_area +
                      order + temp_binned")


# for diatom richness 
mod7 <- brm_logistic("ext_signal ~ diatom_rich_std + 
                      n_units_std + mean_q_std + 
                      pbdb_collections_std + cont_area + 
                      order + temp_binned")

mod8 <- brm_logistic("ext_signal ~ diatom_rich_std + 
                      n_units_std + mean_q_std +
                      shark_collections_std + cont_area +
                      order + temp_binned")

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2, 
                                   mod3, mod4, 
                                   mod5, mod6, 
                                   mod7, mod8,
                                   variable = c("b_d13C_mean_std", 
                                                "b_diatom_rich_std"),
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
                 "pred_trend_productivity_cenozoic.rds"))



