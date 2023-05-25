library(tidyverse)
library(here)
library(brms)
library(tidybayes)
library(patchwork)

# plotting configurations
source(here("R", "config_file.R"))

# prepare temperature data ------------------------------------------------

# load data
dat_ipcc <- tibble(filename = list.files(here("data",
                                              "ipcc")) %>%
                     str_remove("tas_global_") %>%
                     str_remove(".csv")) %>% 
  mutate(file_contents = map(list.files(here("data",
                                             "ipcc"),
                                        full.names = TRUE),  
                             ~ read_csv(.x))) %>% 
  unnest(file_contents)

# reformat
dat_ipcc <- dat_ipcc %>% 
  filter(filename %in% c("Historical", "SSP5_8_5")) %>% 
  select(year = Year, 
         temp_gat_binned = Mean) %>% 
  arrange(desc(year)) %>% 
  mutate(temp_gat_st = temp_gat_binned - lead(temp_gat_binned), 
         temp_gat_lt1 = lead(temp_gat_st), 
         temp_gat_lt2 = lead(temp_gat_st, n = 2), 
         temp_gat_lt3 = lead(temp_gat_st, n = 3), 
         temp_gat_lt4 = lead(temp_gat_st, n = 4))


# prepare iucn-sim data ---------------------------------------------------

# load data
dat_iucn <- read_rds(here("data",
                          "iucn_extinction_times.rds"))

# reformat
dat_iucn <- dat_iucn %>% 
  # add reference time
  mutate(ext_time = ext_time + 2023,   
         year = map(ext_time, 
                          ~ c(2023:.x))) %>% 
  unnest(year) %>% 
  # add extinction signal
  mutate(ext_signal = if_else(ext_time == year, 1, 0)) 



# merge -------------------------------------------------------------------

# combine datasets
dat_merged <- dat_iucn %>% 
  left_join(dat_ipcc) %>% 
  drop_na()



# fit models --------------------------------------------------------------

mod1 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt1")
mod2 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt2")
mod3 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt3")
mod4 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt4")





# average predictions -----------------------------------------------------


# for the future predictions
# set up number of draws 
nr_draws <- 1000

# extract logit averaged over models
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 method = "pseudobma")


# average prediction by model stacking
dat_pred_iucn <- pp_average(mod1, mod2,
                            mod3, mod4,
                            newdata = dat_ipcc %>%
                              filter(year >= 2023),
                            seed = 1708,
                            weights = mod_weights,
                            summary = FALSE,
                            method = "posterior_epred",
                            ndraws = nr_draws) %>%
  as_tibble() %>% 
  mutate(nr_draw = rownames(.)) %>% 
  pivot_longer(cols = contains("V")) %>% 
  add_column(year = rep(dat_ipcc %>% 
                          filter(year >= 2023) %>% 
                          pull(year), nr_draws)) %>% 
  # bin to 10 year
  # mutate(year_bin = cut(year, breaks = seq(2020, 2100, by = 10),
  #                       include.lowest = TRUE,
  #                       labels = FALSE)) %>% 
  group_by(year) %>%
  mean_qi(value) %>% 
  # add_column(year = seq(2020, 2090, by = 10)) %>% 
  select(year, value, .lower, .upper)


# same for fossil baseline
dat_pred_fossil <- pp_average(mod_fossil[[1]], mod_fossil[[2]],
                              mod_fossil[[3]], mod_fossil[[4]],
                              newdata = dat_ipcc %>%
                                filter(year >= 2023) %>% 
                                mutate(temp_gat_binned = temp_gat_binned + 14),
                              seed = 1708,
                              summary = FALSE,
                              method = "posterior_epred",
                              ndraws = nr_draws) %>%
  as_tibble() %>% 
  select(-name) %>% 
  mutate(nr_draw = rownames(.)) %>% 
  pivot_longer(cols = contains("V")) %>% 
  add_column(year = rep(dat_ipcc %>% 
                          filter(year >= 2023) %>% 
                          pull(year), nr_draws)) %>% 
  # bin to 10 year
  mutate(year_bin = cut(year, breaks = seq(2020, 2100, by = 10),
                        include.lowest = TRUE,
                        labels = FALSE)) %>% 
  group_by(year_bin) %>%
  mean_qi(value) %>% 
  add_column(year = seq(2020, 2090, by = 10)) %>% 
  select(year, value, .lower, .upper)


# combine data sources
dat_pred <- dat_pred_iucn %>% 
  # convert to percentage change
  mutate(value = (value - .$value[1]) / .$value[1], 
         .lower = (.lower - .$value[1]) / .$value[1],
         .upper = (.upper - .$value[1]) / .$value[1]) %>% 
  # reformat for joining
  pivot_longer(cols = c(value, .lower, .upper)) %>% 
  add_column(type = "iucn") %>% 
  full_join(dat_pred_fossil %>% # same for other dataframe
              mutate(value = (value - .$value[1]) / .$value[1], 
                     .lower = (.lower - .$value[1]) / .$value[1],
                     .upper = (.upper - .$value[1]) / .$value[1]) %>% 
              pivot_longer(cols = c(value, .lower, .upper)) %>% 
              add_column(type = "fossil")) %>% 
  pivot_wider(names_from = name, values_from = value)


# save predictions
dat_pred %>% 
  write_rds(here("data", 
                 "future_predictions.rds"))


dat_pred %>% 
  group_by(type) %>% 
  mean_qi(value)
