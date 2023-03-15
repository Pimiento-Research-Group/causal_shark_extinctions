# load packages
library(tidyverse)
library(here)
library(dagitty)
library(brms)
library(tidybayes)

# plotting configurations
source(here("R", "config_file.R"))



# load data  ------------------------------------------------------------

# fossil data and environmental proxy data on species level
dat_merged <- read_rds(here("data",
                            "processed_merged_data.rds"))




# total effect adjustment sets --------------------------------------------


# load the graph 
dag <- downloadGraph("dagitty.net/m_UM7hV")


# get adjustments sets for the total effect
adjustmentSets(dag, 
               exposure = "paleotemperature",
               outcome = "extinction risk", 
               effect = "total")

# no need to adjust anything



# fit models --------------------------------------------------------------


# start with deep ocean temperature
# average over potential long-term trends
mod1 <- brm_logistic("ext_signal ~ temp_deep_st:temp_deep_lt1")
mod2 <- brm_logistic("ext_signal ~ temp_deep_st:temp_deep_lt2")
mod3 <- brm_logistic("ext_signal ~ temp_deep_st:temp_deep_lt3")
mod4 <- brm_logistic("ext_signal ~ temp_deep_st:temp_deep_lt4")

# same for global average temperature
mod5 <- brm_logistic("ext_signal ~ temp_gat_st:temp_gat_lt1")
mod6 <- brm_logistic("ext_signal ~ temp_gat_st:temp_gat_lt2")
mod7 <- brm_logistic("ext_signal ~ temp_gat_st:temp_gat_lt3")
mod8 <- brm_logistic("ext_signal ~ temp_gat_st:temp_gat_lt4")



# average models ----------------------------------------------------------

# set up grid to average over
dat_new <- tibble(temp_deep_st = -3:3) %>%
  expand_grid(temp_deep_lt = -3:3) %>% 
  mutate(temp_deep_lt1 = temp_deep_lt,
         temp_deep_lt2 = temp_deep_lt,
         temp_deep_lt3 = temp_deep_lt,
         temp_deep_lt4 = temp_deep_lt,
         temp_gat_st = temp_deep_st,
         temp_gat_lt1 = temp_deep_lt,
         temp_gat_lt2 = temp_deep_lt,
         temp_gat_lt3 = temp_deep_lt,
         temp_gat_lt4 = temp_deep_lt)

# set up number of draws 
nr_draws <- 100

# average prediction by model stacking
dat_pred <- pp_average(mod1, mod2, mod3, mod4, 
                       mod5, mod6, mod7, mod8,
                       newdata = dat_new,
                       seed = 1708,
                       summary = FALSE, 
                       method = "posterior_epred", 
                       ndraws = nr_draws) %>% 
  as_tibble() %>% 
  mutate(nr_draw = rownames(.)) %>% 
  pivot_longer(cols = contains("V")) %>% 
  add_column(temperature_st = rep(dat_new$temp_deep_st, nr_draws), 
             temperature_lt = rep(dat_new$temp_deep_lt, nr_draws)) %>% 
  mutate(trend = if_else(temperature_lt <= 0, "Warming", "Cooling")) %>% 
  group_by(nr_draw, temperature_st, trend) %>%
  mean_qi(value) %>% 
  select(temperature_st, trend, value, nr_draw)

# save predictions
dat_pred %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_paleotemperature.rds"), 
            compress = "gz")

# average over posterior draws
dat_pred_av <- dat_pred %>% 
  group_by(temperature_st, trend) %>% 
  summarise(value = mean(value)) %>% 
  ungroup()

# visualise
plot_paleotemp <- dat_pred %>%
  ggplot(aes(temperature_st, value)) +
  stat_dots(aes(y = as.numeric(ext_signal == 1),
                side = ifelse(ext_signal == 1, "bottom", "top")),
            slab_colour = colour_grey,
            slab_size = 0.5,
            slab_alpha = 0.3,
            scale = 0.2,
            data = dat_merged %>%
              pivot_longer(cols = c(temp_deep_st, temp_gat_st),
                           values_to = "temperature_st",
                           names_to = "temp_name") %>%
              # spread out a bit
              mutate(temperature_st = temperature_st + rnorm(nrow(.), 0, 0.1)) %>%
              filter(between(temperature_st, -3, 3))) +
  geom_line(aes(group = interaction(nr_draw, trend), 
                colour = trend), 
            alpha = 0.05) +
  geom_line(aes(group = trend), 
            colour = "white",
            linewidth = 1, 
            data = dat_pred_av) +
  geom_line(aes(colour = trend),
            linewidth = 0.7, 
            data = dat_pred_av) +
  annotate("text", 
           y = 0.363, x = -1.8, 
           label = "Long-term Cooling", 
           colour = "#9CBABF", 
           size = 10/.pt, 
           alpha = 0.8) +
  annotate(geom = "curve",
           x = -2.45, xend = -2.8,
           y = 0.36, yend = 0.275,
           curvature = 0.4,
           colour = "#9CBABF",
           arrow = arrow(length = unit(.2,"cm")), 
           alpha = 0.5) +
  annotate("text", 
           y = 0.323, x = 1.1, 
           label = "Long-term Warming", 
           colour = "#A76861", 
           size = 10/.pt, 
           alpha = 0.8) +
  annotate(geom = "curve",
           x = 1.75, xend = 2.3,
           y = 0.32, yend = 0.28,
           curvature = -0.2,
           colour = "#A76861",
           arrow = arrow(length = unit(.2,"cm")), 
           alpha = 0.5) +
  scale_color_manual(values = c("#A76861", "#9CBABF")) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  labs(y = "Extinction Risk [%]", 
       x = "Temperature Change [°C]") +
  theme(legend.position = "none")





# save plot ---------------------------------------------------------------



# save plot
ggsave(plot_final, filename = here("figures",
                                   "effect_paleotemperature.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)
