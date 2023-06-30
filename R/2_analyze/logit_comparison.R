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
                                 method = "pseudobma")


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
  # add stages
  left_join(dat_stages) %>%
  add_column(data_set = "Stages") %>% 
  select(-bin) %>% 
  # same for genus resolution
  bind_rows(dat_pred_genus %>% 
              # add stages
              left_join(dat_stages) %>%
              add_column(data_set = "Genus") %>% 
              select(-bin)) %>% 
  # add cenozoic resolution data
  bind_rows(dat_pred_ceno %>% 
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

# create plot
plot_logit <- dat_pred_full %>%
  mutate(data_set = factor(data_set, 
                           levels = c("Stages", 
                                      "Genus", 
                                      "1myr", 
                                      "Modern", 
                                      "Future"))) %>% 
  ggplot(aes(myr, value, group = data_set)) +
  geom_hline(yintercept = 0, 
             colour = "grey20") +
  geom_linerange(aes(xmin = xmin,
                     xmax = xmax),
                 colour = "grey75") +
  geom_linerange(aes(ymin = .lower,
                     ymax = .upper),
                 colour = "grey75") +
  geom_point(aes(fill = data_set), 
             size = 2, 
             alpha = 0.8, 
             shape = 21, 
             colour = "grey20") +
  scale_x_reverse() +
  scale_fill_manual(values = c("#4C634C", 
                               colour_coral, 
                               colour_purple, 
                               "#DF5B08", 
                               "#e71ed5"), 
                    name = NULL) +
  scale_y_continuous(limits = c(-7, 2), 
                     breaks = seq(-6, 2, by = 2)) +
  coord_geo(xlim = c(150, 0), 
            dat = list("epochs", "periods"),
            pos = list("b", "b"),
            alpha = 0.2, 
            height = unit(0.8, "line"), 
            size = list(5/.pt, 9/.pt),
            lab_color = "grey20", 
            color = "grey50", 
            abbrv = list(TRUE, FALSE), 
            fill = "white",
            expand = TRUE, 
            lwd = list(0.1, 0.2)) +
  labs(x = "Million years", 
       y = "Temperature dependancy\n(log-odds)") +
  theme(legend.position = "bottom") 
   



# same for r-squared
# combine data 
dat_rsq <- list(dat_r_deep, 
     dat_r_genus, 
     dat_r_ceno, 
     dat_r_modern, 
     dat_r_future) %>% 
  bind_rows() %>% 
  add_column(x_axis = factor(c(rep("Fossil", 3),
                               "Modern",
                               "Future"),
                             levels = c("Fossil",
                                        "Modern",
                                        "Future")), 
             data_set = factor(c("deep", 
                                 "genus", 
                                 "ceno", 
                                 "modern", 
                                 "future"), 
                               levels = c("deep", 
                                          "genus", 
                                          "ceno", 
                                          "modern", 
                                          "future")))
# and save
dat_rsq %>% 
  write_rds(here("data", 
                 "r_squared_full.rds"))

# visualise
plot_rsq <- dat_rsq %>%
  ggplot(aes(x_axis, R2)) +
  geom_linerange(aes(ymin = .lower,
                     ymax = .upper, 
                     group = data_set),
                 colour = "grey75", 
                 position = position_dodge(width = 0.2)) +
  geom_point(aes(fill = data_set),
             alpha = 0.8, 
             shape = 21, 
             size = 2, 
             colour = "grey20", 
             position = position_dodge(width = 0.2)) +
  scale_fill_manual(values = c("#4C634C", 
                               colour_coral, 
                               colour_purple, 
                               "#DF5B08", 
                               "#e71ed5"), 
                    name = NULL) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2)) +
  labs(x = NULL, 
       y = bquote("Bayesian" ~R^2)) +
  theme(legend.position = "none")



# patch together
plot_full <- plot_logit / plot_rsq +
  plot_annotation(tag_levels = "a") +
  plot_layout(heights = c(3, 1))


# save plot
ggsave(plot_full, filename = here("figures",
                                  "logit_r_squared.png"), 
       width = image_width, height = image_height*1.5,
       units = image_units, 
       bg = "white", device = ragg::agg_png)     
