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





# total effect adjustment sets --------------------------------------------


# load the graph 
dag <- downloadGraph("dagitty.net/m_UM7hV")


# get adjustments sets for temperature, for the total effect
adjustmentSets(dag, 
               exposure = "productivity",
               outcome = "extinction risk", 
               effect = "total")

# two sets of adjustments
# latitude, outcrop area, sampling effort, shelf area, temperature
# latitude, sampling effort, sea level, shelf area, temperature


# fit models --------------------------------------------------------------


## first adjustment set 
# for d13C
mod1 <- brm_logistic("ext_signal ~ d13C_std + 
                     latitude_pref_abs + cont_area +
                     shark_collections_std + temp_gat_binned")

mod2 <- brm_logistic("ext_signal ~ d13C_std + 
                     latitude_pref_abs + cont_area + 
                     pbdb_collections_std + temp_gat_binned")

mod3 <- brm_logistic("ext_signal ~ d13C_std + 
                     latitude_pref_abs + cont_area + 
                     pbdb_collections_std + temp_deep_binned")

mod4 <- brm_logistic("ext_signal ~ d13C_std + 
                     latitude_pref_abs + cont_area + 
                     shark_collections_std + temp_deep_binned")

# for Sr ratio 
mod5 <- brm_logistic("ext_signal ~ sr_value_std + 
                     latitude_pref_abs + cont_area + 
                     shark_collections_std + temp_gat_binned")

mod6 <- brm_logistic("ext_signal ~ sr_value_std + 
                     latitude_pref_abs + cont_area + 
                     pbdb_collections_std + temp_gat_binned")

mod7 <- brm_logistic("ext_signal ~ sr_value_std + 
                     latitude_pref_abs + cont_area + 
                     pbdb_collections_std + temp_deep_binned")

mod8 <- brm_logistic("ext_signal ~ sr_value_std + 
                     latitude_pref_abs + cont_area + 
                     shark_collections_std + temp_deep_binned")


## second adjustment set 
# for d13C
mod9 <- brm_logistic("ext_signal ~ d13C_std + latitude_pref_abs + 
                     sea_level + cont_area + 
                     shark_collections_std + temp_gat_binned")

mod10 <- brm_logistic("ext_signal ~ d13C_std + latitude_pref_abs +
                      sea_level + cont_area +
                      pbdb_collections_std + temp_gat_binned")

mod11 <- brm_logistic("ext_signal ~ d13C_std + latitude_pref_abs + 
                      sea_level + cont_area +
                      pbdb_collections_std + temp_deep_binned")

mod12 <- brm_logistic("ext_signal ~ d13C_std + latitude_pref_abs + 
                      sea_level + cont_area +
                      shark_collections_std + temp_deep_binned")

# for Sr ratio 
mod13 <- brm_logistic("ext_signal ~ sr_value_std + latitude_pref_abs + 
                     sea_level + cont_area + 
                     shark_collections_std + temp_gat_binned")

mod14 <- brm_logistic("ext_signal ~ sr_value_std + latitude_pref_abs +
                      sea_level + cont_area +
                      pbdb_collections_std + temp_gat_binned")

mod15 <- brm_logistic("ext_signal ~ sr_value_std + latitude_pref_abs + 
                      sea_level + cont_area +
                      pbdb_collections_std + temp_deep_binned")

mod16 <- brm_logistic("ext_signal ~ sr_value_std + latitude_pref_abs + 
                      sea_level + cont_area +
                      shark_collections_std + temp_deep_binned")



# average models ----------------------------------------------------------

# set up grid to average over
dat_new <- tibble(d13C_std = seq(-2, 2, by = 0.2), 
                  sr_value_std = seq(-2, 2, by = 0.2), 
                  # keep adjustments fixed at average
                  latitude_pref_abs = mean(dat_merged$latitude_pref_abs), 
                  sea_level = mean(dat_merged$sea_level), 
                  cont_area = mean(dat_merged$cont_area), 
                  shark_collections_std = mean(dat_merged$shark_collections_std), 
                  pbdb_collections_std = mean(dat_merged$pbdb_collections_std), 
                  temp_deep_binned = mean(dat_merged$temp_deep_binned), 
                  temp_gat_binned = mean(dat_merged$temp_gat_binned))  

# set up number of draws 
nr_draws <- 100

# average prediction by model stacking
dat_pred <- pp_average(mod1, mod2, mod3, mod4, 
                       mod5, mod6, mod7, mod8,
                       mod9, mod10, mod11, mod12,
                       mod13, mod14, mod15, mod16,
                       newdata = dat_new,
                       seed = 1708,
                       summary = FALSE, 
                       method = "posterior_epred", 
                       ndraws = nr_draws) %>% 
  as_tibble() %>% 
  mutate(nr_draw = rownames(.)) %>% 
  pivot_longer(cols = contains("V")) %>% 
  add_column(productivity = rep(dat_new$d13C_std, nr_draws)) %>% 
  group_by(nr_draw, productivity) %>%
  mean_qi(value) 

# average over posterior draws
dat_pred_av <- dat_pred %>% 
  group_by(productivity) %>% 
  summarise(value = mean(value))

# visualise
plot_prod <- dat_pred %>% 
  ggplot(aes(productivity, value)) +
  stat_dots(aes(y = as.numeric(ext_signal == 1), 
                side = ifelse(ext_signal == 1, "bottom", "top")),
            slab_colour = colour_grey,
            slab_size = 0.5,
            slab_alpha = 0.3,
            scale = 0.2,
            data = dat_merged %>% 
              pivot_longer(cols = c(d13C_std, sr_value_std), 
                           values_to = "productivity", 
                           names_to = "prod_name") %>% 
              # spread out a bit
              mutate(productivity = productivity + rnorm(nrow(.), 0, 0.1)) %>% 
              filter(between(productivity, -2, 2))) +
  geom_line(aes(group = nr_draw), 
            alpha = 0.2, color = colour_yellow) +
  geom_line(colour = "white",
            linewidth = 1.3, 
            data = dat_pred_av) +
  geom_line(colour = colour_yellow,
            linewidth = 1, 
            data = dat_pred_av) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  coord_cartesian(ylim = c(0, 1), 
                  xlim = c(-2, 2)) +
  labs(y = "Extinction Risk [%]", 
       x = "Global Primary Productivity [std]")


# estimate trend ----------------------------------------------------------

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2,
                                   mod3, mod4,
                                   mod5, mod6,
                                   mod7, mod8,
                                   variable = c("b_d13C_std", "b_sr_value_std"),
                                   seed = 1708,
                                   ndraws =  1e4, 
                                   missing = NA) %>% 
  pivot_longer(cols = everything(), 
               names_to = "coef_name", 
               values_to = "coef_val") %>% 
  drop_na(coef_val)



# visualise
plot_prod_beta <- dat_pred_post %>%
  ggplot(aes(coef_val)) +
  geom_vline(xintercept = 0, 
             linewidth = 0.5,
             colour = "grey40") +
  stat_slab(shape = 21, 
            slab_colour = NA, 
            slab_fill = colour_yellow, 
            slab_alpha = 0.7,
            point_size = 3, 
            point_fill = "white",
            point_colour = colour_yellow) +
  annotate("text", 
           label = "\u03B2", 
           x = -0.6, 
           y = 0.85, 
           size = 10/.pt, 
           colour = "grey40") +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = 0, limits = c(-1, 3.2)) +
  labs(y = NULL, 
       x = NULL) +
  theme(plot.background = elementalist::element_rect_round(radius = unit(0.85, "cm"), 
                                                           color = "#FFEFE1"), 
        axis.ticks = element_blank())




# patch together and save -------------------------------------------------


# patch together
plot_final <- plot_prod +
  inset_element(plot_prod_beta, 
                left = 0.1, 
                bottom = 0.6, 
                right = 0.25,
                top = 0.8) 



# save plot
ggsave(plot_final, filename = here("figures",
                                   "effect_productivity.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)
