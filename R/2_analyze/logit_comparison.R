library(tidyverse)
library(here)
library(brms)
library(tidybayes)
library(deeptime)
library(patchwork)

# plotting configurations
source(here("R", "config_file.R"))

# set up number of draws 
nr_draws <- 100



# deep-time models --------------------------------------------------------------



# fossil data and environmental proxy data on species level
dat_merged <- read_rds(here("data",
                            "processed_merged_data.rds")) %>% 
  # merge bins
  mutate(bin_merg = cut(bin, 
                        breaks = 69:94))


# start with deep ocean temperature
# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt1 + (temp_deep_binned|bin_merg)")
mod2 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt2 + (temp_deep_binned|bin_merg)")
mod3 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt3 + (temp_deep_binned|bin_merg)")
mod4 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt4 + (temp_deep_binned|bin_merg)")

# same for global average temperature
mod5 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt1 + (temp_gat_binned|bin_merg)")
mod6 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt2 + (temp_gat_binned|bin_merg)")
mod7 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt3 + (temp_gat_binned|bin_merg)")
mod8 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt4 + (temp_gat_binned|bin_merg)")


# extract logit averaged over models

# average prediction by model stacking
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 mod5, mod6,
                                 mod7, mod8,
                                 method = "pseudobma")


# perform model averaging
dat_pred_deep <- distinct(dat_merged, bin_merg) %>%
  drop_na(bin_merg) %>% 
  pull(bin_merg) %>% 
  map_df(.f = ~ pp_average(mod1, mod2,
                           mod3, mod4,
                           mod5, mod6,
                           mod7, mod8,
                           newdata = filter(dat_merged, 
                                            bin_merg == .x),
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

# same for bayesian r-squared
dat_r_deep <- list(mod1, mod2,
                   mod3, mod4,
                   mod5, mod6,
                   mod7, mod8) %>%
  map_df(~ .x %>% 
           bayes_R2(summary = FALSE, 
                    ndraws = 1000) %>% 
           as_tibble()) %>% 
  median_qi(R2)


# genus models --------------------------------------------------------

# read in genus resolution data
dat_merged <- read_rds(here("data",
                            "processed_fossil_data_genus.rds")) %>% 
  # merge bins
  mutate(bin_merg = cut(bin, 
                        breaks = 69:94))


# start with deep ocean temperature
# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt1 + (temp_deep_binned|bin_merg)")
mod2 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt2 + (temp_deep_binned|bin_merg)")
mod3 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt3 + (temp_deep_binned|bin_merg)")
mod4 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt4 + (temp_deep_binned|bin_merg)")

# same for global average temperature
mod5 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt1 + (temp_gat_binned|bin_merg)")
mod6 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt2 + (temp_gat_binned|bin_merg)")
mod7 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt3 + (temp_gat_binned|bin_merg)")
mod8 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt4 + (temp_gat_binned|bin_merg)")


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
dat_pred_genus <- distinct(dat_merged, bin_merg) %>%
  drop_na(bin_merg) %>% 
  pull(bin_merg) %>% 
  map_df(.f = ~ pp_average(mod1, mod2,
                           mod3, mod4,
                           mod5, mod6,
                           mod7, mod8,
                           newdata = filter(dat_merged, 
                                            bin_merg == .x),
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


# same for bayesian r-squared
dat_r_genus <- list(mod1, mod2,
                    mod3, mod4,
                    mod5, mod6,
                    mod7, mod8) %>%
  map_df(~ .x %>% 
           bayes_R2(summary = FALSE, 
                    ndraws = 1000) %>% 
           as_tibble()) %>% 
  median_qi(R2)

# cenozoic models --------------------------------------------------------------


# read in cenozoic resolution data
dat_merged <- read_rds(here("data",
                            "processed_fossil_data_cenozoic.rds")) %>% 
  # merge bins
  mutate(bin_merg = cut(bin, 
                        breaks = 1:66))


# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt1 + (temp_binned|bin_merg)")
mod2 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt2 + (temp_binned|bin_merg)")
mod3 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt3 + (temp_binned|bin_merg)")
mod4 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt4 + (temp_binned|bin_merg)")


# extract logit averaged over models
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 method = "pseudobma")


# perform model averaging
dat_pred_ceno <- distinct(dat_merged, bin_merg) %>%
  drop_na(bin_merg) %>% 
  pull(bin_merg) %>% 
  map_df(.f = ~ pp_average(mod1, mod2,
                           mod3, mod4,
                           newdata = filter(dat_merged, 
                                            bin_merg == .x),
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

# same for bayesian r-squared
dat_r_ceno <- list(mod1, mod2,
                   mod3, mod4) %>%
  map_df(~ .x %>% 
           bayes_R2(summary = FALSE, 
                    ndraws = 1000) %>% 
           as_tibble()) %>% 
  median_qi(R2)

# modern models --------------------------------------------------------------


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
  arrange(desc(year)) %>% 
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

# same for bayesian r-squared
dat_r_modern <- list(mod1, mod2,
                     mod3, mod4) %>%
  map_df(~ .x %>% 
           bayes_R2(summary = FALSE, 
                    ndraws = 1000) %>% 
           as_tibble()) %>% 
  median_qi(R2)


# future models -------------------------------------------------------

# prepare temperature data 
# load data
dat_ipcc_fut <- tibble(filename = list.files(here("data",
                                              "ipcc")) %>%
                     str_remove("tas_global_") %>%
                     str_remove(".csv")) %>% 
  mutate(file_contents = map(list.files(here("data",
                                             "ipcc"),
                                        full.names = TRUE),  
                             ~ read_csv(.x))) %>% 
  unnest(file_contents) %>% # reformat
  filter(filename %in% c("Historical", "SSP5_8_5")) %>% 
  select(year = Year, 
         temp_gat_binned = Mean) %>% 
  arrange(desc(year)) %>% 
  mutate(temp_gat_st = temp_gat_binned - lead(temp_gat_binned), 
         temp_gat_lt1 = lead(temp_gat_st), 
         temp_gat_lt2 = lead(temp_gat_st, n = 2), 
         temp_gat_lt3 = lead(temp_gat_st, n = 3), 
         temp_gat_lt4 = lead(temp_gat_st, n = 4))

# prepare iucn-sim data 

# load data
dat_iucn_sim <- read_rds(here("data",
                          "iucn_extinction_times.rds")) %>% 
  # add reference time
  mutate(ext_time = ext_time + 2023,   
         year = map(ext_time, 
                    ~ c(2023:.x))) %>% 
  unnest(year) %>% 
  # add extinction signal
  mutate(ext_signal = if_else(ext_time == year, 1, 0)) 

# merge

# combine datasets
dat_merged <- dat_iucn_sim %>% 
  left_join(dat_ipcc_fut) %>% 
  drop_na()


# fit models

mod1 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt1")
mod2 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt2")
mod3 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt3")
mod4 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt4")


# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 method = "pseudobma")


# average posterior draws by model stacking
dat_pred_future <- posterior_average(mod1, mod2,
                                     mod3, mod4,
                                     variable = c("b_temp_gat_binned"),
                                     seed = 1708,
                                     ndraws =  1e4,
                                     missing = NA,
                                     weights = mod_weights) %>% 
  as_tibble() %>% 
  mean_qi(b_temp_gat_binned) %>% 
  select(value = b_temp_gat_binned, .lower, .upper)


# save predictions
dat_pred_future %>% 
  write_rds(here("data", 
                 "logits", 
                 "logit_future.rds"))


# same for bayesian r-squared
dat_r_future <- list(mod1, mod2,
                     mod3, mod4) %>%
  map_df(~ .x %>% 
           bayes_R2(summary = FALSE, 
                    ndraws = 1000) %>% 
           as_tibble()) %>% 
  median_qi(R2)


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
  mutate(bin = as.character(bin) %>% 
           str_sub(2, 3) %>% 
           as.numeric()) %>% 
  # add stages
  left_join(dat_stages) %>%
  add_column(data_set = "Stages") %>% 
  select(-bin) %>% 
  # same for genus resolution
  bind_rows(dat_pred_genus %>% 
              mutate(bin = as.character(bin) %>% 
                       str_sub(2, 3) %>% 
                       as.numeric()) %>%
              # add stages
              left_join(dat_stages) %>%
              add_column(data_set = "Genus") %>% 
              select(-bin)) %>% 
  # add cenozoic resolution data
  bind_rows(dat_pred_ceno %>% 
              mutate(bin = 65:1) %>%
              rename(myr = bin) %>% 
              add_column(data_set = "1myr") %>% 
              mutate(xmin = myr - 0.5, 
                     xmax = myr + 0.5)) %>% 
  # add modern data
  bind_rows(dat_pred_modern %>% 
              add_column(data_set = "Modern", 
                         myr = 0, 
                         xmin = 0, 
                         xmax = 0)) %>% 
  # add future data
  bind_rows(dat_pred_future %>% 
              add_column(data_set = "Future", 
                         myr = 0, 
                         xmin = -0.00002, 
                         xmax = 0.00002)) 


# save predictions
dat_pred_full %>% 
  write_rds(here("data", 
                 "logits", 
                 "logit_full.rds"))




# supplementary plot ------------------------------------------------------


# summarise per hypo- and hyperthermal
plot_log_hyp <- dat_pred_full %>%
  filter(data_set != "Future") %>% 
  # assign thermal status
  mutate(therm_stat = case_when(between(myr, 66, 72) | between(myr, 34, 38) ~ "Hypothermal",
                                between(myr, 94, 101) | between(myr, 56, 62) | between(myr, 12, 16) ~ "Hyperthermal",
                                .default = "Background")) %>% 
  group_by(therm_stat, data_set) %>% 
  median_qi(value) %>% 
  filter(therm_stat != "Background") %>%
  # select(therm_stat, value, data_set) %>%
  # pivot_wider(names_from = therm_stat,
  #             values_from = value) %>%
  # mutate(perc_inc = (Hypothermal - Hyperthermal)/ Hypothermal) %>% 
  # median_qi(perc_inc)
  ggplot(aes(therm_stat)) +
  geom_linerange(aes(ymin = .lower, 
                     ymax = .upper, 
                     colour = data_set), 
                 position = position_dodge(width = 0.4)) +
  geom_point(aes(y = value, 
                 fill = data_set), 
             shape = 21, 
             colour = "white", 
             size = 2, 
             stroke = 0.5, 
             position = position_dodge(width = 0.4)) +
  annotate("text",
           y = -2.3, 
           x = 1.5, 
           size = 8/.pt, 
           label = "+9%\n[-63%, +23%]", 
           colour = "grey40", 
           fontface = "bold") +
  annotate("rect", 
           xmin = 0.8, 
           xmax = 1.2, 
           ymin = -Inf, 
           ymax = Inf, 
           fill = "#C75E6B", 
           alpha = 0.2) + 
  annotate("rect", 
           xmin = 1.8, 
           xmax = 2.2,
           ymin = -Inf, 
           ymax = Inf, 
           fill = "#169199", 
           alpha = 0.2) + 
  labs(y = "Temperature dependancy\n[log-odds]", 
       x = NULL) +
  scale_fill_manual(values = c(colour_purple, 
                               colour_coral, 
                               "#4C634C")) +
  scale_colour_manual(values = c(colour_purple, 
                                 colour_coral, 
                                 "#4C634C")) +
  scale_y_continuous(breaks = c(0, -2, -4), 
                     limits = c(-5.4, 0)) +
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank())

# read in second plot for figure c
plot_thermal <- read_rds(here("data",
                              "predictions",
                              "thermal_plot.rds"))


# plot_full <- 
plot_full <- plot_log_hyp +
  plot_thermal +
  plot_annotation(tag_levels = "a")
 

 




