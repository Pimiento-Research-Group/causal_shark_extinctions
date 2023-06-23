library(tidyverse)
library(here)
library(brms)
library(tidybayes)
library(deeptime)


# plotting configurations
source(here("R", "config_file.R"))

# set up number of draws 
nr_draws <- 100



# read data ---------------------------------------------------------------



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
                              0, 1))

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
dat_merged <- left_join(dat_iucn, dat_ipcc) %>% 
  mutate(metab = if_else(species %in% c(
    "Carcharodon carcharias", 
    "Lamna nasus", 
    "Lamna ditropis", 
    "Isurus oxyrinchus", 
    "Isurus paucus"
  ), 
  "mesotherm", 
  "ectotherm"))


# model -------------------------------------------------------------------


# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp:metab + temp_st:temp_lt1")
mod2 <- brm_logistic("ext_signal ~ temp:metab + temp_st:temp_lt2")
mod3 <- brm_logistic("ext_signal ~ temp:metab + temp_st:temp_lt3")
mod4 <- brm_logistic("ext_signal ~ temp:metab + temp_st:temp_lt4")

# extract model weights
mod_weights <- model_weights(mod1, 
                             mod2, 
                             mod3, 
                             mod4, 
                             weights = "pseudobma")

# create grid to average over
dat_new <- tibble(temp = seq(0.5, 1.3, by = 0.1),
                  temp_st = mean(dat_merged$temp_st),
                  temp_lt1 = mean(dat_merged$temp_lt1),
                  temp_lt2 = temp_lt1,
                  temp_lt3 = temp_lt1,
                  temp_lt4 = temp_lt1) %>% 
  expand_grid(metab = c("mesotherm", 
                        "ectotherm"))


# extract model averaged estimates
dat_av <- pp_average(mod1, mod2,
             mod3, mod4,
             newdata = dat_new,
             seed = 1708,
             summary = FALSE, 
             method = "posterior_epred", 
             weights = mod_weights,
             ndraws = 100)  %>% 
  as_tibble() %>% 
  mutate(nr_draw = rownames(.)) %>% 
  pivot_longer(cols = contains("V")) %>% 
  add_column(temperature = rep(dat_new$temp, 100), 
             metab = rep(dat_new$metab, 100)) %>% 
  group_by(temperature, metab) %>%
  mean_qi(value) %>% 
  select(temperature, metab, value, .lower, .upper)  

# save data
dat_av %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_metabolism.rds"))

# extract logits
posterior_average(mod1, mod2,
           mod3, mod4,
           seed = 1708,
           weights = mod_weights, 
           ndraws = 100) %>% 
  as_tibble() %>% 
  mean_qi() %>% 
  View()

# visualise ---------------------------------------------------------------


dat_av %>%  
  ggplot(aes(temperature, value)) +
  geom_ribbon(aes(ymin = .lower, 
                  ymax = .upper, 
                  fill = metab), 
              alpha = 0.3) +
  geom_line(aes(colour = metab))



