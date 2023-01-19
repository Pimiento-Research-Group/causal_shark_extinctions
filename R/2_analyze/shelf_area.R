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
               exposure = "shelf area",
               outcome = "extinction risk", 
               effect = "total")

# only need to adjust for sea level


# fit models --------------------------------------------------------------


# one model with sea level 
mod1 <- brm_logistic("ext_signal ~ cont_area + sea_level")



# extract predictions ----------------------------------------------------------

# set up grid to average over
dat_pred <- tibble(cont_area = seq(0, 0.3, by = 0.02), 
                   # keep adjustments at average
                   sea_level = mean(dat_merged$sea_level)) %>% 
  # add posterior draws
  add_epred_draws(mod1, 
                  ndraws = 100) %>% 
  ungroup() %>% 
  select(cont_area, .epred, .draw)

# save predictions
dat_pred %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_shelf_area.rds"))

dat_pred_av <- dat_pred %>% 
  group_by(cont_area) %>% 
  summarise(.epred = mean(.epred))

# visualise
plot_shelf <- dat_pred %>% 
  ggplot(aes(cont_area, .epred)) +
  stat_dots(aes(y = as.numeric(ext_signal == 1), 
                side = ifelse(ext_signal == 1, "bottom", "top")),
            slab_colour = colour_grey,
            slab_size = 0.5,
            slab_alpha = 0.3,
            scale = 0.2,
            data = dat_merged %>% 
              # spread out a bit
              mutate(cont_area = cont_area + rnorm(nrow(.), 0, 0.01))) +
  geom_line(aes(group = .draw), 
            alpha = 0.2, color = "#BD8D9E") +
  geom_line(colour = "white",
            linewidth = 1.3, 
            data = dat_pred_av) +
  geom_line(colour = "#BD8D9E",
            linewidth = 1, 
            data = dat_pred_av) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  coord_cartesian(ylim = c(0, 1), 
                  xlim = c(0, 0.25)) +
  labs(y = "Extinction Risk [%]", 
       x = "Continental Flooding [Proportion of Earth's Surface]")


# estimate trend ----------------------------------------------------------

# average posterior draws by model stacking
dat_pred_post <- as_draws_df(mod1, 
                             variable = "b_cont_area") %>% 
  as_tibble() %>% 
  select(coef_val = b_cont_area) %>% 
  slice_sample(n = 1e4)


# save trend predictions
dat_pred_post %>% 
  write_rds(here("data", 
                 "predictions", 
                 "pred_trend_shelf_area.rds"))

# visualise
plot_shelf_beta <- dat_pred_post %>%
  ggplot(aes(coef_val)) +
  geom_vline(xintercept = 0, 
             linewidth = 0.5,
             colour = "grey40") +
  stat_slab(shape = 21, 
            slab_colour = NA, 
            slab_fill = "#BD8D9E", 
            slab_alpha = 0.7,
            point_size = 3, 
            point_fill = "white",
            point_colour = "#BD8D9E") +
  annotate("text", 
           label = "\u03B2", 
           x = -17, 
           y = 0.85, 
           size = 10/.pt, 
           colour = "grey40") +
  scale_y_continuous(breaks = NULL) +
  scale_x_continuous(breaks = 0) +
  labs(y = NULL, 
       x = NULL) +
  theme(plot.background = elementalist::element_rect_round(radius = unit(0.85, "cm"), 
                                                           color = "#BD8D9E"), 
        axis.ticks = element_blank())




# patch together and save -------------------------------------------------


# patch together
plot_final <- plot_shelf +
  inset_element(plot_shelf_beta, 
                left = 0.1, 
                bottom = 0.6, 
                right = 0.25,
                top = 0.8) 



# save plot
ggsave(plot_final, filename = here("figures",
                                   "effect_shelf_area.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)
