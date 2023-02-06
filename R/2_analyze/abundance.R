# load packages
library(tidyverse)
library(here)
library(dagitty)
library(brms)
library(tidybayes)
library(patchwork)

# plotting configurations
source(here("R", "config_file.R"))



# load data  ------------------------------------------------------------

# fossil data and environmental proxy data on species level
dat_merged <- read_rds(here("data",
                            "processed_merged_data.rds")) 

# select taxa where we have an extinction event
dat_merged <- dat_merged %>% 
  filter(order %in% (dat_merged %>%
                       count(order,
                             ext_signal) %>%
                       filter(ext_signal == 1) %>%
                       pull(order)))

# total effect adjustment sets --------------------------------------------


# load the graph 
dag <- downloadGraph("dagitty.net/m_UM7hV")


# get adjustments sets for the total effect
adjustmentSets(dag, 
               exposure = "abundance",
               outcome = "extinction risk", 
               effect = "total")

# two sets of adjustments 
# { geographic range, taxonomic identity }
# { outcrop area, paleotemperature, preservation potential, 
# productivity, sampling effort, shelf area, taxonomic identity, temperature }

# fit models --------------------------------------------------------------


# first adjustment set
mod1 <- brm_logistic("ext_signal ~ abund + geo_dist_std + order")
mod2 <- brm_logistic("ext_signal ~ abund + range_lat_std + order")

# number of genera
mod3 <- brm_logistic("ext_signal ~ n_genus + geo_dist_std + order")
mod4 <- brm_logistic("ext_signal ~ n_genus + range_lat_std + order")

# the second adjustment set is too big to average over potential parameters 
# (would need 256 models), select 4 of these 256 models randomly instead
set.seed(1708)
model_subset <- list(abundance = c("abund", "n_genus"),
                     preservation = c("mean_q_std"),
                     shelf_area = c("cont_area"),
                     taxonomy = c("order"),
                     outcrop_area = c("outcrop_area_std", "n_units_std"),
                     paleotemperature = c("temp_deep_st:temp_deep_lt1",
                                          "temp_deep_st:temp_deep_lt2",
                                          "temp_deep_st:temp_deep_lt3",
                                          "temp_deep_st:temp_deep_lt4",
                                          "temp_gat_st:temp_gat_lt1",
                                          "temp_gat_st:temp_gat_lt2",
                                          "temp_gat_st:temp_gat_lt3",
                                          "temp_gat_st:temp_gat_lt4"),
                     productivity = c("d13C_std", "sr_value_std"),
                     sampling = c("pbdb_collections_std", "shark_collections_std"),
                     temperature = c("temp_gat_binned", "temp_deep_binned")) %>% 
  expand.grid() %>% 
  slice_sample(n = 4)

# second adjustment subset
model_subset[1, ]
mod5 <- brm_logistic("ext_signal ~ n_genus + mean_q_std + cont_area + order +
                     n_units_std + temp_deep_st:temp_deep_lt3 + sr_value_std +
                     shark_collections_std + temp_gat_binned")

model_subset[2, ]
mod6 <- brm_logistic("ext_signal ~ n_genus + mean_q_std + cont_area + order +
                     outcrop_area_std + temp_deep_st:temp_deep_lt3 + d13C_std +
                     shark_collections_std + temp_gat_binned")

model_subset[3, ]
mod7 <- brm_logistic("ext_signal ~ abund + mean_q_std + cont_area + order +
                     n_units_std + temp_deep_st:temp_deep_st:temp_deep_lt2 + sr_value_std +
                     shark_collections_std + temp_deep_binned")

model_subset[4, ]
mod8 <- brm_logistic("ext_signal ~ abund + mean_q_std + cont_area + order +
                     outcrop_area_std + temp_deep_st:temp_deep_lt3 + sr_value_std +
                     shark_collections_std + temp_gat_binned")


# average models ----------------------------------------------------------

# set up grid to average over
dat_new <- tibble(abund = 0:60, 
                  n_genus = abund) %>%
  # average over taxonomy
  expand_grid(order = unique(dat_merged$order)) %>% 
  # keep adjustment at average
  add_column(range_lat_std = mean(dat_merged$range_lat_std), 
             geo_dist_std = mean(dat_merged$geo_dist_std), 
             mean_q_std = mean(dat_merged$mean_q_std), 
             cont_area = mean(dat_merged$cont_area),
             n_units_std = mean(dat_merged$n_units_std), 
             outcrop_area_std = mean(dat_merged$outcrop_area_std), 
             temp_deep_st = mean(dat_merged$temp_deep_st), 
             temp_deep_lt3 = mean(dat_merged$temp_deep_lt3), 
             temp_deep_lt2 = mean(dat_merged$temp_deep_lt2), 
             d13C_std = mean(dat_merged$d13C_std), 
             sr_value_std = mean(dat_merged$sr_value_std),
             shark_collections_std = mean(dat_merged$shark_collections_std), 
             temp_gat_binned = mean(dat_merged$temp_gat_binned), 
             temp_deep_binned = mean(dat_merged$temp_deep_binned))

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
  add_column(abund = rep(dat_new$abund, nr_draws)) %>% 
  group_by(nr_draw, abund) %>%
  mean_qi(value) %>% 
  select(abund, nr_draw, value)

# save predictions
dat_pred %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_abundance.rds"), 
            compress = "gz")

# average over posterior draws
dat_pred_av <- dat_pred %>% 
  group_by(abund) %>% 
  summarise(value = mean(value))

# visualise
plot_abund <- dat_pred %>%
  ggplot(aes(abund, value)) +
  stat_dots(aes(y = as.numeric(ext_signal == 1), 
                side = ifelse(ext_signal == 1, "bottom", "top")),
            slab_colour = colour_grey,
            slab_size = 0.5,
            slab_alpha = 0.3,
            scale = 0.2,
            data = dat_merged %>% 
              pivot_longer(cols = c(abund, n_genus), 
                           values_to = "abund", 
                           names_to = "abund_name") %>% 
              # spread out a bit
              mutate(abund = abund + rnorm(nrow(.), 0, 0.7))) +
  geom_line(aes(group = nr_draw), 
            alpha = 0.2, color = "#BD6E4E") +
  geom_line(colour = "white",
            linewidth = 1.3, 
            data = dat_pred_av) +
  geom_line(colour = "#BD6E4E",
            linewidth = 1, 
            data = dat_pred_av) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  coord_cartesian(xlim = c(0, 60)) +
  labs(y = "Extinction Risk [%]", 
       x = "Abundance [Number of Occurrences]")

# estimate trend ----------------------------------------------------------

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2,
                                   mod3, mod4,
                                   mod5, mod6,
                                   mod7, mod8,
                                   variable = c("b_abund", "b_n_genus"),
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
                 "pred_trend_abund.rds"), 
            compress = "gz")

# visualise
plot_abund_beta <- dat_pred_post %>%
  ggplot(aes(coef_val)) +
  geom_vline(xintercept = 0, 
             linewidth = 0.5,
             colour = "grey40") +
  stat_slab(shape = 21, 
            slab_colour = NA, 
            slab_fill = "#BD6E4E", 
            slab_alpha = 0.7,
            point_size = 3, 
            point_fill = "white",
            point_colour = "#BD6E4E") +
  annotate("text", 
           label = "\u03B2", 
           x = -0.4, 
           y = 0.85, 
           size = 10/.pt, 
           colour = "grey40") +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = 0, limits = c(-0.45, 0.04)) +
  labs(y = NULL, 
       x = NULL) +
  theme(plot.background = elementalist::element_rect_round(radius = unit(0.85, "cm"), 
                                                           color = "#BD6E4E"), 
        axis.ticks = element_blank())




# patch together and save -------------------------------------------------


# patch together
plot_final <- plot_abund +
  inset_element(plot_abund_beta, 
                left = 0.75, 
                bottom = 0.6, 
                right = 0.9,
                top = 0.8) 

# save plot
ggsave(plot_final, filename = here("figures",
                                   "effect_abundance.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)
