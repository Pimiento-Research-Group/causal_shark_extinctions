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


# get adjustments sets for the total effect
adjustmentSets(dag, 
               exposure = "sea level",
               outcome = "extinction risk", 
               effect = "total")

# only need to adjust for temperature


# fit models --------------------------------------------------------------


# average over temperature estimates 
mod1 <- brm_logistic("ext_signal ~ sea_level + temp_gat_binned")
mod2 <- brm_logistic("ext_signal ~ sea_level + temp_deep_binned")



# average models ----------------------------------------------------------

# set up grid to average over
dat_new <- tibble(sea_level = seq(-70, 70, by = 2)) %>%
  # keep adjustment at average
  add_column(temp_gat_binned = mean(dat_merged$temp_gat_binned), 
             temp_deep_binned = mean(dat_merged$temp_deep_binned))

# set up number of draws 
nr_draws <- 100

# average prediction by model stacking
dat_pred <- pp_average(mod1, mod2,
                       newdata = dat_new,
                       seed = 1708,
                       summary = FALSE, 
                       method = "posterior_epred", 
                       ndraws = nr_draws) %>% 
  as_tibble() %>% 
  mutate(nr_draw = rownames(.)) %>% 
  pivot_longer(cols = contains("V")) %>% 
  add_column(sea_level = rep(dat_new$sea_level, nr_draws)) %>% 
  group_by(nr_draw, sea_level) %>%
  mean_qi(value) %>% 
  select(sea_level, value, nr_draw)

# save predictions
dat_pred %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_sea_level.rds"))

# average over posterior draws
dat_pred_av <- dat_pred %>% 
  group_by(sea_level) %>% 
  summarise(value = mean(value))

# visualise
plot_sea <- dat_pred %>%
  ggplot(aes(sea_level, value)) +
  stat_dots(aes(y = as.numeric(ext_signal == 1), 
                side = ifelse(ext_signal == 1, "bottom", "top")),
            slab_colour = colour_grey,
            slab_size = 0.5,
            slab_alpha = 0.3,
            scale = 0.2,
            data = dat_merged %>% 
              pivot_longer(cols = sea_level, 
                           values_to = "sea_level", 
                           names_to = "sea_level_name") %>% 
              # spread out a bit
              mutate(sea_level = sea_level + rnorm(nrow(.), 0, 1))) +
  geom_line(aes(group = nr_draw), 
            alpha = 0.2, color = colour_mint) +
  geom_line(colour = "white",
            linewidth = 1.3, 
            data = dat_pred_av) +
  geom_line(colour = colour_mint,
            linewidth = 1, 
            data = dat_pred_av) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  coord_cartesian(xlim = c(-70, 70)) +
  labs(y = "Extinction Risk [%]", 
       x = "Global Mean Sea Level [m]")


# estimate trend ----------------------------------------------------------

# average posterior draws by model stacking
dat_pred_post <- posterior_average(mod1, mod2,
                                   variable = "b_sea_level",
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
                 "pred_trend_sea_level.rds"))

# visualise
plot_sea_beta <- dat_pred_post %>%
  ggplot(aes(coef_val)) +
  geom_vline(xintercept = 0, 
             linewidth = 0.5,
             colour = "grey40") +
  stat_slab(shape = 21, 
            slab_colour = NA, 
            slab_fill = colour_mint, 
            slab_alpha = 0.7,
            point_size = 3, 
            point_fill = "white",
            point_colour = colour_mint) +
  annotate("text", 
           label = "\u03B2", 
           x = -0.03, 
           y = 0.85, 
           size = 10/.pt, 
           colour = "grey40") +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = 0) +
  labs(y = NULL, 
       x = NULL) +
  theme(plot.background = elementalist::element_rect_round(radius = unit(0.85, "cm"), 
                                                           color = colour_mint), 
        axis.ticks = element_blank())




# patch together and save -------------------------------------------------


# patch together
plot_final <- plot_sea +
  inset_element(plot_sea_beta, 
                left = 0.75, 
                bottom = 0.6, 
                right = 0.9,
                top = 0.8) 

# save plot
ggsave(plot_final, filename = here("figures",
                                   "effect_sea_level.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)
