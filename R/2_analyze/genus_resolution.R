# load packages
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
  # bin to stages myr
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
  # use only occurrences from genera that could get binned to stages
  filter(accepted_name %in% dat_occ_binned$genus) %>% 
  rename(genus = accepted_name)




# merge and combine -------------------------------------------------------

dat_merged <- full_join(dat_genus, 
                        dat_preservation) %>% 
  # combine with those parameters that are the same on genus level
  full_join(read_rds(here("data",
                          "processed_merged_data.rds")) %>% 
              select(-c(ext_signal, mean_q_std, species))) %>% 
  # remove missing values
  left_join(dat_preservation %>%
              group_by(genus) %>%
              summarise(mean_q_av = mean(mean_q))) %>% 
  mutate(mean_q = if_else(is.na(mean_q), mean_q_av, mean_q), 
         mean_q_std = scale(mean_q)[,1]) %>% 
  drop_na()

# save dataset
dat_merged %>% 
  write_rds(here("data",
                 "processed_fossil_data_genus.rds"))



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
                     pbdb_collections_std + sea_level + 
                     cont_area + temp_gat_binned")

mod2 <- brm_logistic("ext_signal ~ d13C_std +
                     outcrop_area_std + sea_level + 
                     cont_area + temp_gat_binned")

mod3 <- brm_logistic("ext_signal ~ d13C_std + 
                     pbdb_collections_std + sea_level +
                     cont_area + temp_deep_binned")

mod4 <- brm_logistic("ext_signal ~ d13C_std + 
                     outcrop_area_std + sea_level + 
                     cont_area + temp_deep_binned")

# for Sr ratio 
mod5 <- brm_logistic("ext_signal ~ sr_mean_std + 
                     pbdb_collections_std + sea_level + 
                     cont_area + temp_gat_binned")

mod6 <- brm_logistic("ext_signal ~ sr_mean_std + 
                     outcrop_area_std + sea_level + 
                     cont_area + temp_gat_binned")

mod7 <- brm_logistic("ext_signal ~ sr_mean_std + 
                     pbdb_collections_std + sea_level +
                     cont_area + temp_deep_binned")

mod8 <- brm_logistic("ext_signal ~ sr_mean_std + 
                     outcrop_area_std + sea_level + 
                     cont_area + temp_deep_binned")


## second adjustment set 
# outcrop area, preservation potential, sampling effort, shelf area, 
# taxonomic identity, temperature
# for d13C
mod9 <- brm_logistic("ext_signal ~ d13C_std + 
                     outcrop_area_std + mean_q_std + 
                     pbdb_collections_std + cont_area + 
                     order + temp_gat_binned")

mod10 <- brm_logistic("ext_signal ~ d13C_std + 
                      outcrop_area_std + mean_q_std +
                      shark_collections_std + cont_area +
                      order + temp_gat_binned")

mod11 <- brm_logistic("ext_signal ~ d13C_std + 
                      outcrop_area_std + mean_q_std +
                      pbdb_collections_std + cont_area +
                      order + temp_deep_binned")

mod12 <- brm_logistic("ext_signal ~ d13C_std + 
                      outcrop_area_std + mean_q_std +
                      shark_collections_std + cont_area +
                      order + temp_deep_binned")

# for Sr ratio 
mod13 <- brm_logistic("ext_signal ~ sr_mean_std + 
                      outcrop_area_std + mean_q_std + 
                      pbdb_collections_std + cont_area + 
                      order + temp_gat_binned")

mod14 <- brm_logistic("ext_signal ~ sr_mean_std + 
                      outcrop_area_std + mean_q_std +
                      shark_collections_std + cont_area +
                      order + temp_gat_binned")

mod15 <- brm_logistic("ext_signal ~ sr_mean_std + 
                      outcrop_area_std + mean_q_std +
                      pbdb_collections_std + cont_area +
                      order + temp_deep_binned")

mod16 <- brm_logistic("ext_signal ~ sr_mean_std + 
                      outcrop_area_std + mean_q_std +
                      shark_collections_std + cont_area +
                      order + temp_deep_binned")




# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2, mod3, mod4, 
                                   mod5, mod6, mod7, mod8,
                                   mod9, mod10, mod11, mod12,
                                   mod13, mod14, mod15, mod16,
                                   variable = c("b_d13C_std", 
                                                "b_sr_mean_std"),
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
                 "pred_temperature_genus.rds"), 
            compress = "gz")

