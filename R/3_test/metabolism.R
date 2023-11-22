library(tidyverse)
library(here)
library(brms)
library(tidybayes)

# plotting configurations
source(here("R", "config_file.R"))

# set up number of draws 
nr_draws <- 100



# deep-time models --------------------------------------------------------------



# fossil data and environmental proxy data on species level
dat_merged <- read_rds(here("data",
                            "processed_merged_data.rds")) %>% 
  # remove families with less than 3 occurrences
  group_by(family) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 3)

  
  

# start with deep ocean temperature
# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt1 + (temp_deep_binned|family)")
mod2 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt2 + (temp_deep_binned|family)")
mod3 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt3 + (temp_deep_binned|family)")
mod4 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt4 + (temp_deep_binned|family)")

# same for global average temperature
mod5 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt1 + (temp_gat_binned|family)")
mod6 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt2 + (temp_gat_binned|family)")
mod7 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt3 + (temp_gat_binned|family)")
mod8 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt4 + (temp_gat_binned|family)")


# extract logit averaged over models

# average prediction by model stacking
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 mod5, mod6,
                                 mod7, mod8,
                                 method = "pseudobma")


# perform model averaging for every family
dat_pred_deep <- distinct(dat_merged, family) %>%
  drop_na(family) %>% 
  pull(family) %>% 
  map_df(.f = ~ pp_average(mod1, mod2,
                           mod3, mod4,
                           mod5, mod6,
                           mod7, mod8,
                           newdata = filter(dat_merged, 
                                            family == .x),
                           seed = 1708,
                           summary = FALSE,
                           method = "posterior_linpred",
                           ndraws = nr_draws, 
                           weights = mod_weights) %>% 
           as_tibble() %>% 
           pivot_longer(cols = contains("V")) %>% 
           add_column(family = .x) %>% 
           select(value, family))




# genus models --------------------------------------------------------

# read in genus resolution data
dat_merged <- read_rds(here("data",
                            "processed_fossil_data_genus.rds")) %>% 
  # remove families with less than 3 occurrences
  group_by(family) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 3)


# start with deep ocean temperature
# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt1 + (temp_deep_binned|family)")
mod2 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt2 + (temp_deep_binned|family)")
mod3 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt3 + (temp_deep_binned|family)")
mod4 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt4 + (temp_deep_binned|family)")

# same for global average temperature
mod5 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt1 + (temp_gat_binned|family)")
mod6 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt2 + (temp_gat_binned|family)")
mod7 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt3 + (temp_gat_binned|family)")
mod8 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt4 + (temp_gat_binned|family)")


# extract logit averaged over models

# average prediction by model stacking
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 mod5, mod6,
                                 mod7, mod8,
                                 method = "pseudobma",
                                 cores = parallelly::availableCores())


# perform model averaging
dat_pred_genus <- distinct(dat_merged, family) %>%
  drop_na(family) %>% 
  pull(family) %>% 
  map_df(.f = ~ pp_average(mod1, mod2,
                           mod3, mod4,
                           mod5, mod6,
                           mod7, mod8,
                           newdata = filter(dat_merged, 
                                            family == .x),
                           seed = 1708,
                           summary = FALSE,
                           method = "posterior_linpred",
                           ndraws = nr_draws, 
                           weights = mod_weights) %>% 
           as_tibble() %>% 
           pivot_longer(cols = contains("V")) %>% 
           add_column(family = .x) %>% 
           select(value, family))


# cenozoic models --------------------------------------------------------------


# read in cenozoic resolution data
dat_merged <- read_rds(here("data",
                            "processed_fossil_data_cenozoic.rds")) %>% 
  # remove families with less than 3 occurrences
  group_by(family) %>% 
  mutate(n = n()) %>% 
  ungroup() %>% 
  filter(n >= 3)


# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt1 + (temp_binned|family)")
mod2 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt2 + (temp_binned|family)")
mod3 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt3 + (temp_binned|family)")
mod4 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt4 + (temp_binned|family)")


# extract logit averaged over models
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 method = "pseudobma")


# perform model averaging
dat_pred_ceno <- distinct(dat_merged, family) %>%
  drop_na(family) %>% 
  pull(family) %>% 
  map_df(.f = ~ pp_average(mod1, mod2,
                           mod3, mod4,
                           newdata = filter(dat_merged, 
                                            family == .x),
                           seed = 1708,
                           summary = FALSE,
                           method = "posterior_linpred",
                           ndraws = nr_draws, 
                           weights = mod_weights) %>% 
           as_tibble() %>% 
           pivot_longer(cols = contains("V")) %>% 
           add_column(family = .x) %>% 
           select(value, family))




# merge together ----------------------------------------------------------


# merge together
dat_family <- dat_pred_genus %>% 
  add_column(scale = "genus") %>% 
  full_join(dat_pred_deep %>% 
              add_column(scale = "stages")) %>% 
  full_join(dat_pred_ceno %>% 
              add_column(scale = "ceno")) 

# save data
dat_family %>% 
  write_rds(here(here("data", 
                      "logits", 
                      "logit_family_full.rds")), 
            compress = "gz")





