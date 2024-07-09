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
dag <- downloadGraph("dagitty.net/mjiV5Qf")


# get adjustments sets for the total effect
adjustmentSets(dag, 
               exposure = "temperature",
               outcome = "extinction risk", 
               effect = "total")

# need to account for paleotemperature



# fit models --------------------------------------------------------------


# start with deep ocean temperature
# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt1")
mod2 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt2")
mod3 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt3")
mod4 <- brm_logistic("ext_signal ~ temp_deep_binned + temp_deep_st:temp_deep_lt4")

# same for global average temperature
mod5 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt1")
mod6 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt2")
mod7 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt3")
mod8 <- brm_logistic("ext_signal ~ temp_gat_binned + temp_gat_st:temp_gat_lt4")


# # save global average temperature models for reuse later
# list(mod1, mod2, 
#      mod3, mod4,
#      mod5, mod6, 
#      mod7, mod8) %>% 
#   write_rds(here("data", 
#                  "fossil_temp_models.rds"), 
#             compress = "gz")



# average models ----------------------------------------------------------

# set up grid to average over
dat_new <- tibble(temp_deep_binned = 0:25) %>%
  expand_grid(temp_deep_st = -3:3, 
              temp_deep_lt = -3:3) %>% 
  mutate(temp_deep_lt1 = temp_deep_lt,
         temp_deep_lt2 = temp_deep_lt,
         temp_deep_lt3 = temp_deep_lt,
         temp_deep_lt4 = temp_deep_lt,
         temp_gat_binned = temp_deep_binned, 
         temp_gat_st = temp_deep_st,
         temp_gat_lt1 = temp_deep_lt,
         temp_gat_lt2 = temp_deep_lt,
         temp_gat_lt3 = temp_deep_lt,
         temp_gat_lt4 = temp_deep_lt)

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
  add_column(temperature = rep(dat_new$temp_deep_binned, nr_draws)) %>% 
  group_by(nr_draw, temperature) %>%
  mean_qi(value) %>% 
  select(temperature, value, nr_draw)


# save predictions
dat_pred %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_temperature_stage.rds"), 
            compress = "gz")

# calculate marginal effect by model stacking
dat_deep <- posterior_average(mod1, mod2, mod3, mod4, 
                              seed = 1708,
                              ndraws = nr_draws)

dat_gat <- posterior_average(mod5, mod6, mod7, mod8, 
                             seed = 1708,
                             ndraws = nr_draws)
# average
dat_change <- dat_deep %>% 
  as_tibble() %>% 
  rename(b_temp = b_temp_deep_binned) %>% 
  bind_rows(dat_gat %>% 
              as_tibble() %>% 
              rename(b_temp = b_temp_gat_binned)) %>% 
  mutate(prob_change = plogis(b_temp  + 1) - plogis(b_temp)) %>% 
  mean_qi(prob_change)



dat_gat_change <- dat_gat %>% 
  as_tibble() %>% 
  mutate(prob_change = plogis(b_temp_gat_binned  + 1) - plogis(b_temp_gat_binned)) %>% 
  summarise(avg_change = mean(prob_change)) %>% 
  pull(avg_change)

(dat_deep_change + dat_gat_change)/ 2

# average over posterior draws
dat_pred_av <- dat_pred %>% 
  group_by(temperature) %>% 
  summarise(value = mean(value))

# visualise
plot_temp <- dat_pred %>%
  ggplot(aes(temperature, value)) +
  stat_dots(aes(y = as.numeric(ext_signal == 1), 
                side = ifelse(ext_signal == 1, "bottom", "top")),
    slab_colour = colour_grey,
    slab_size = 0.5,
    slab_alpha = 0.3,
    scale = 0.2,
    data = dat_merged %>% 
      pivot_longer(cols = c(temp_deep_binned, temp_gat_binned), 
                   values_to = "temperature", 
                   names_to = "temp_name") %>% 
      # spread out a bit
      mutate(temperature = temperature + rnorm(nrow(.), 0, 0.1))) +
  geom_line(aes(group = nr_draw), 
            alpha = 0.2, color = colour_coral) +
  geom_line(colour = "white",
            linewidth = 1.3, 
            data = dat_pred_av) +
  geom_line(colour = colour_coral,
            linewidth = 1, 
            data = dat_pred_av) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  coord_cartesian(ylim = c(0, 1), 
                  xlim = c(0, 25)) +
  labs(y = "Extinction Risk [%]", 
       x = "Global Temperature [°C]")


# estimate trend ----------------------------------------------------------

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2,
                                   mod3, mod4,
                                   mod5, mod6,
                                   mod7, mod8,
                                   variable = c("b_temp_gat_binned", "b_temp_deep_binned"),
                                   seed = 1708,
                                   ndraws =  1e4, 
                                   missing = NA) %>% 
  pivot_longer(cols = everything(), 
               names_to = "coef_name", 
               values_to = "coef_val") %>% 
  drop_na(coef_val)

# save predictions
dat_pred_post %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_temperature.rds"), 
            compress = "gz")

# visualise
plot_temp_beta <- dat_pred_post %>%
  ggplot(aes(coef_val)) +
  geom_vline(xintercept = 0, 
             linewidth = 0.5,
             colour = "grey40") +
  stat_slab(shape = 21, 
            slab_colour = NA, 
            slab_fill = colour_coral, 
            slab_alpha = 0.7,
            point_size = 3, 
            point_fill = "white",
            point_colour = colour_coral) +
  annotate("text", 
           label = "\u03B2", 
           x = -0.175, 
           y = 0.85, 
           size = 10/.pt, 
           colour = "grey40") +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = 0, limits = c(-0.22, 0.1)) +
  labs(y = NULL, 
       x = NULL) +
  theme(plot.background = elementalist::element_rect_round(radius = unit(0.85, "cm"), 
                                                           color = colour_coral), 
        axis.ticks = element_blank())




# patch together and save -------------------------------------------------


# patch together
plot_final <- plot_temp +
  inset_element(plot_temp_beta, 
                left = 0.75, 
                bottom = 0.6, 
                right = 0.9,
                top = 0.8) 


# # save plot
# ggsave(plot_final, filename = here("figures",
#                                    "effect_temperature.png"), 
#        width = image_width, height = image_height,
#        units = image_units, 
#        bg = "white", device = ragg::agg_png)
