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


# read fossil models
mod_fossil <- read_rds(here("data", 
                             "fossil_temp_models.rds"))



# average predictions -----------------------------------------------------


# for the future predictions
# set up number of draws 
nr_draws <- 1000

# average prediction by model stacking
dat_pred_iucn <- pp_average(mod1, mod2,
                            mod3, mod4,
                            newdata = dat_ipcc %>%
                              filter(year >= 2023),
                            seed = 1708,
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
  mutate(year_bin = cut(year, breaks = seq(2020, 2100, by = 10),
                        include.lowest = TRUE,
                        labels = FALSE)) %>% 
  group_by(year_bin) %>%
  mean_qi(value) %>% 
  add_column(year = seq(2020, 2090, by = 10)) %>% 
  select(year, value, .lower, .upper)


# save predictions
dat_pred %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_temperature.rds"), 
            compress = "gz")

dat_ipcc %>% 
  filter(year >= 2023) %>% 
  mutate(temp_gat_binned = temp_gat_binned + 14) %>% 
  add_epred_draws(., mod_fossil[[1]]) %>% 
  group_by(year) %>% 
  median_qi(.epred) %>% 
  mutate(.epred = (.epred - .$.epred[1]) / .$.epred[1], 
         .lower = (.lower - .$.epred[1]) / .$.epred[1],
         .upper = (.upper - .$.epred[1]) / .$.epred[1]) %>% 
  ggplot(aes(year, .epred)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_line() +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              alpha = 0.2) +
  # scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
  #                    labels = c("0", "20", "40", "60", "80", "100")) +
  # coord_cartesian(ylim = c(0, 1)) +
  coord_cartesian(ylim = c(-0.5, 0.5)) +
  labs(y = "Extinction Risk [%]", 
       x = "Year")



dat_pred_iucn %>% 
  mutate(value = (value - .$value[1]) / .$value[1], 
         .lower = (.lower - .$value[1]) / .$value[1],
         .upper = (.upper - .$value[1]) / .$value[1]) %>% 
  ggplot(aes(year, value)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_ribbon(aes(ymin = .lower,
                  ymax = .upper),
              alpha = 0.2) +
  geom_line() +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(y = "Extinction Risk [%]", 
       x = "Year")
