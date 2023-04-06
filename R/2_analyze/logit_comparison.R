library(tidyverse)
library(here)
library(brms)
library(tidybayes)
library(patchwork)

# plotting configurations
source(here("R", "config_file.R"))


# load deep-time data ---------------------------------------------------------------


# fossil data and environmental proxy data on species level
dat_merged <- read_rds(here("data",
                            "processed_merged_data.rds")) %>% 
  mutate(bin = as.factor(bin))


# fit deep-time models --------------------------------------------------------------


# start with deep ocean temperature
# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt1 + (1|bin)")
mod2 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt2 + (1|bin)")
mod3 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt3 + (1|bin)")
mod4 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt4 + (1|bin)")

# same for global average temperature
mod5 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt1 + (1|bin)")
mod6 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt2 + (1|bin)")
mod7 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt3 + (1|bin)")
mod8 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt4 + (1|bin)")


# extract logit averaged over models
# set up number of draws 
nr_draws <- 100

# average prediction by model stacking
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 mod5, mod6,
                                 mod7, mod8,
                                 method = "pseudobma",
                                 cores = parallelly::availableCores())


# perform model averaging
dat_pred_deep <- distinct(dat_merged, bin) %>%
  mutate(bin = as.integer(as.character(bin))) %>% 
  pull(bin) %>% 
  map_df(.f = ~ pp_average(mod1, mod2,
                           mod3, mod4,
                           mod5, mod6,
                           mod7, mod8,
                           newdata = dat_merged %>% 
                             filter(bin == .x),
                           seed = 1708,
                           summary = FALSE,
                           method = "posterior_linpred",
                           ndraws = nr_draws, 
                           weights = mod_weights) %>% 
           as_tibble() %>% 
           pivot_longer(cols = contains("V")) %>% 
           mean_qi(value) %>% 
           add_column(bin = .x) %>% 
           select(value, .lower, .upper, bin))

# save predictions
dat_pred_deep %>% 
  write_rds(here("data", 
                 "logits", 
                 "logit_deep.rds"))


# fit cenozoic models --------------------------------------------------------------


# read in cenozoic resolution data
dat_merged <- read_rds(here("data",
                            "processed_fossil_data_cenozoic.rds")) %>% 
  mutate(bin = as.factor(bin))


# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt1 + (1|bin)")
mod2 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt2 + (1|bin)")
mod3 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt3 + (1|bin)")
mod4 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt4 + (1|bin)")


# extract logit averaged over models
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 method = "pseudobma")


# perform model averaging
dat_pred_ceno <- distinct(dat_merged, bin) %>%
  mutate(bin = as.integer(as.character(bin))) %>% 
  pull(bin) %>% 
  map_df(.f = ~ pp_average(mod1, mod2,
                           mod3, mod4,
                           newdata = dat_merged %>% 
                             filter(bin == .x),
                           seed = 1708,
                           summary = FALSE,
                           method = "posterior_linpred",
                           ndraws = nr_draws, 
                           weights = mod_weights) %>% 
           as_tibble() %>% 
           pivot_longer(cols = contains("V")) %>% 
           mean_qi(value) %>% 
           add_column(bin = .x) %>% 
           select(value, .lower, .upper, bin))


# save predictions
dat_pred_ceno %>% 
  write_rds(here("data", 
                 "logits", 
                 "logit_ceno.rds"))



# fit modern models --------------------------------------------------------------

dat_pred_ceno %>% 
  ggplot(aes(bin, value)) +
  geom_hline(yintercept = 0) +
  geom_pointrange(aes(ymin = .lower, 
                      ymax = .upper))
