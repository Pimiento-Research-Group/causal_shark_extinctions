library(tidyverse)
library(here)
library(brms)
library(tidybayes)
library(patchwork)

# plotting configurations
source(here("R", "config_file.R"))

# set up number of draws 
nr_draws <- 100




# fit deep-time models --------------------------------------------------------------



# fossil data and environmental proxy data on species level
dat_merged <- read_rds(here("data",
                            "processed_merged_data.rds")) %>% 
  mutate(bin = as.factor(bin))


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



# fit genus models --------------------------------------------------------

# read in cenozoic resolution data
dat_merged <- read_rds(here("data",
                            "processed_fossil_data_genus.rds")) %>% 
  mutate(bin = as.factor(bin))


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

# average prediction by model stacking
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 mod5, mod6,
                                 mod7, mod8,
                                 method = "pseudobma",
                                 cores = parallelly::availableCores())


# perform model averaging
dat_pred_genus <- distinct(dat_merged, bin) %>%
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
dat_pred_genus %>% 
  write_rds(here("data", 
                 "logits", 
                 "logit_genus.rds"))


# fit 10myr models --------------------------------------------------------

# read in cenozoic resolution data
dat_merged <- read_rds(here("data",
                            "processed_fossil_data_10myr.rds")) %>% 
  mutate(bin = as.factor(bin))


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

# average prediction by model stacking
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 mod5, mod6,
                                 mod7, mod8,
                                 method = "pseudobma",
                                 cores = parallelly::availableCores())


# perform model averaging
dat_pred_10myr <- distinct(dat_merged, bin) %>%
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
dat_pred_10myr %>% 
  write_rds(here("data", 
                 "logits", 
                 "logit_10myr.rds"))


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


# read iucn data
dat_iucn <- read_delim(here("data",
                            "iucn",
                            "CHONDRICHTHYES_iucn_history.txt"))

# reformat
dat_iucn <- dat_iucn %>% 
  pivot_longer(cols = -c(species), 
               names_to = "year", 
               values_to = "iucn_status") %>% 
  drop_na(iucn_status) %>% 
  mutate(ext_signal = if_else(iucn_status %in% c("LR/lc", "LC",
                                                 "LR / cd", "NT" ,
                                                 "LR / nt"),
                              1, 0))

# read temperature data
dat_ipcc <- tibble(filename = list.files(here("data",
                                  "ipcc")) %>%
         str_remove("tas_global_") %>%
         str_remove(".csv")) %>% 
  mutate(file_contents = map(list.files(here("data",
                                             "ipcc"),
                                        full.names = TRUE),  
                             ~ read_csv(.x))) %>% 
  filter(filename %in% c("Historical", "SSP1_1_9")) %>% 
  unnest(file_contents) %>% 
  select(year = Year, temp = Mean) %>% 
  filter(year %in% unique(dat_iucn$year))

# reformat
dat_ipcc <- dat_ipcc %>% 
  mutate(temp_st = temp - lead(temp), 
         temp_lt1 = lead(temp_st), 
         temp_lt2 = lead(temp_st, n = 2), 
         temp_lt3 = lead(temp_st, n = 3), 
         temp_lt4 = lead(temp_st, n = 4)) %>% 
  fill(contains("temp"), .direction = "downup") %>% 
  mutate(year = as.character(year))
  

# merge together
dat_merged <- left_join(dat_iucn, dat_ipcc)


# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp + temp_st:temp_lt1")
mod2 <- brm_logistic("ext_signal ~ temp + temp_st:temp_lt2")
mod3 <- brm_logistic("ext_signal ~ temp + temp_st:temp_lt3")
mod4 <- brm_logistic("ext_signal ~ temp + temp_st:temp_lt4")

# average posterior draws by model stacking
dat_pred_modern <- posterior_average(mod1, mod2,
                                     mod3, mod4,
                                     variable = c("b_temp"),
                                     seed = 1708,
                                     ndraws =  1e4,
                                     missing = NA) %>% 
  as_tibble() %>% 
  mean_qi(b_temp) %>% 
  select(value = b_temp, .lower, .upper)


# save predictions
dat_pred_modern %>% 
  write_rds(here("data", 
                 "logits", 
                 "logit_modern.rds"))





# visualise ---------------------------------------------------------------


# load 95 bin Phanerozoic time scale based on the stratigraphic stages 
# of Gradstein et al. 2020
load(here("data", 
          "stages.Rdata"))

# clean up for binning
dat_stages <- stages %>% 
  as_tibble() %>% 
  select(bin = stg, myr = mid, 
         xmin = bottom, 
         xmax = top)

# merge datasets
dat_pred_full <- dat_pred_deep %>% 
  # add stages
  left_join(dat_stages) %>%
  add_column(data_set = "deep_time") %>% 
  select(-bin) %>% 
  # same for genus resolution
  bind_rows(dat_pred_genus %>% 
              # add stages
              left_join(dat_stages) %>%
              add_column(data_set = "genus") %>% 
              select(-bin)) %>% 
  # add 10 myr resolution data
  bind_rows(dat_pred_10myr %>% 
              add_column(data_set = "10myr") %>% 
              arrange(bin) %>% 
              add_column(myr = seq(5, 135, by = 10)) %>% 
              mutate(xmin = myr - 5, 
                     xmax = myr + 5) %>% 
              select(-bin)) %>% 
  # add cenozoic resolution data
  bind_rows(dat_pred_ceno %>% 
              rename(myr = bin) %>% 
              add_column(data_set = "near_time") %>% 
              mutate(xmin = myr - 0.5, 
                     xmax = myr + 0.5)) %>% 
  # add modern data
  bind_rows(dat_pred_modern %>% 
              add_column(data_set = "modern", 
                         myr = 0, 
                         xmin = 0, 
                         xmax = 0)) 


# save predictions
dat_pred_full %>% 
  write_rds(here("data", 
                 "logits", 
                 "logit_full.rds"))


dat_pred_full %>% 
  ggplot(aes(myr, value, group = data_set)) +
  geom_hline(yintercept = 0, 
             colour = "grey50") +
  geom_linerange(aes(xmin = xmin,
                     xmax = xmax),
                 colour = "grey40") +
  geom_linerange(aes(ymin = .lower,
                     ymax = .upper),
                 colour = "grey40") +
  geom_point(aes(fill = data_set), 
             size = 3, 
             alpha = 0.7, 
             shape = 21, 
             colour = "grey20") +
  scale_x_continuous(breaks = seq(0, 150, 10), 
                     trans = "reverse") +
  theme(legend.position = "bottom")




                     