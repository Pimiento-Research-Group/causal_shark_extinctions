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
                            "processed_merged_data.rds")) 


# start with deep ocean temperature
# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt1 + (temp_deep_binned|order)")
mod2 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt2 + (temp_deep_binned|order)")
mod3 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt3 + (temp_deep_binned|order)")
mod4 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt4 + (temp_deep_binned|order)")

# same for global average temperature
mod5 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt1 + (temp_gat_binned|order)")
mod6 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt2 + (temp_gat_binned|order)")
mod7 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt3 + (temp_gat_binned|order)")
mod8 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt4 + (temp_gat_binned|order)")


# extract logit averaged over models

# average prediction by model stacking
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 mod5, mod6,
                                 mod7, mod8,
                                 method = "pseudobma")


# perform model averaging
dat_pred_deep <- distinct(dat_merged, order) %>%
  drop_na(order) %>% 
  pull(order) %>% 
  map_df(.f = ~ pp_average(mod1, mod2,
                           mod3, mod4,
                           mod5, mod6,
                           mod7, mod8,
                           newdata = filter(dat_merged, 
                                            order == .x),
                           seed = 1708,
                           summary = FALSE,
                           method = "posterior_linpred",
                           ndraws = nr_draws, 
                           weights = mod_weights) %>% 
           as_tibble() %>% 
           pivot_longer(cols = contains("V"), 
                        values_to = "logit_val") %>% 
           add_column(order = .x) %>% 
           select(-name))

# save predictions
dat_pred_deep %>% 
  write_rds(here("data", 
                 "logits", 
                 "logit_order_deep.rds"), 
            compress = "gz")


# genus models --------------------------------------------------------

# read in genus resolution data
dat_merged <- read_rds(here("data",
                            "processed_fossil_data_genus.rds")) 


# start with deep ocean temperature
# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt1 + (temp_deep_binned|order)")
mod2 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt2 + (temp_deep_binned|order)")
mod3 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt3 + (temp_deep_binned|order)")
mod4 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt4 + (temp_deep_binned|order)")

# same for global average temperature
mod5 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt1 + (temp_gat_binned|order)")
mod6 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt2 + (temp_gat_binned|order)")
mod7 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt3 + (temp_gat_binned|order)")
mod8 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt4 + (temp_gat_binned|order)")


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
dat_pred_genus <- distinct(dat_merged, order) %>%
  drop_na(order) %>% 
  pull(order) %>% 
  map_df(.f = ~ pp_average(mod1, mod2,
                           mod3, mod4,
                           mod5, mod6,
                           mod7, mod8,
                           newdata = filter(dat_merged, 
                                            order == .x),
                           seed = 1708,
                           summary = FALSE,
                           method = "posterior_linpred",
                           ndraws = nr_draws, 
                           weights = mod_weights) %>% 
           as_tibble() %>% 
           pivot_longer(cols = contains("V"), 
                        values_to = "logit_val") %>% 
           add_column(order = .x) %>% 
           select(-name))

# save predictions
dat_pred_genus %>% 
  write_rds(here("data", 
                 "logits", 
                 "logit_order_genus.rds"), 
            compress = "gz")


# cenozoic models --------------------------------------------------------------


# read in cenozoic resolution data
dat_merged <- read_rds(here("data",
                            "processed_fossil_data_cenozoic.rds")) 


# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt1 + (temp_binned|order)")
mod2 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt2 + (temp_binned|order)")
mod3 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt3 + (temp_binned|order)")
mod4 <- brm_logistic("ext_signal ~ temp_binned + temp_st:temp_lt4 + (temp_binned|order)")


# extract logit averaged over models
# get model weights 
mod_weights <- loo_model_weights(mod1, mod2,
                                 mod3, mod4,
                                 method = "pseudobma")


# perform model averaging
dat_pred_ceno <- distinct(dat_merged, order) %>%
  drop_na(order) %>% 
  pull(order) %>% 
  map_df(.f = ~ pp_average(mod1, mod2,
                           mod3, mod4,
                           newdata = filter(dat_merged, 
                                            order == .x),
                           seed = 1708,
                           summary = FALSE,
                           method = "posterior_linpred",
                           ndraws = nr_draws, 
                           weights = mod_weights) %>% 
           as_tibble() %>% 
           pivot_longer(cols = contains("V"), 
                        values_to = "logit_val") %>% 
           add_column(order = .x) %>% 
           select(-name))


# save predictions
dat_pred_ceno %>% 
  write_rds(here("data", 
                 "logits", 
                 "logit_order_ceno.rds"), 
            compress = "gz")



# merge and summarise -----------------------------------------------------


# merge together
dat_order <- dat_pred_genus %>% 
  add_column(scale = "genus") %>% 
  full_join(dat_pred_deep %>% 
              add_column(scale = "stages")) %>% 
  full_join(dat_pred_ceno %>% 
              add_column(scale = "ceno")) %>% 
  group_by(order) %>% 
  median_qi(logit_val)

# save data
dat_order %>% 
  write_rds(here(here("data", 
                      "logits", 
                      "logit_order.rds")))


# visualise ---------------------------------------------------------------

plot_order <- dat_pred_genus %>%
  add_column(scale = "genus") %>% 
  full_join(dat_pred_deep %>% 
              add_column(scale = "stages")) %>% 
  full_join(dat_pred_ceno %>% 
              add_column(scale = "ceno")) %>% 
  group_by(order, scale) %>% 
  median_qi(logit_val) %>% 
  # reorder
  group_by(order) %>% 
  mutate(mean_val = mean(logit_val)) %>% 
  ungroup() %>% 
  filter(order != "incertae sedis") %>% 
  # abbreviate order for nicer plotting
  mutate(order = str_replace_all(order, "formes", "."), 
         order = fct_reorder(order, mean_val)) %>% 
  ggplot(aes(logit_val, order)) +
  geom_vline(xintercept = 0, 
             colour = "grey20", 
             linetype = "dashed") +
  geom_linerange(aes(xmin = .lower, 
                     xmax = .upper, 
                     group = scale), 
                 position = position_dodge(width = 1), 
                 colour = colour_grey) +
  geom_point(aes(fill = scale), 
             shape = 21,
             size = 2, 
             stroke = 0.6, 
             colour = "grey20",
             position = position_dodge(width = 1)) +
  scale_fill_manual(values = rev(c("#4C634C",
                               colour_coral,
                               colour_purple)),
                    labels = rev(c("Species - Stages",
                               "Genera - Stages",
                               "Species - Cenozoic subset")),
                    name = NULL) +
  labs(y = NULL, 
       x = "Temperature dependancy\n[log-odds]") +
  theme(legend.position = "top")

# save plot
ggsave(plot_order, filename = here("figures",
                                   "supplement",
                                   "logit_per_order.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png) 

  