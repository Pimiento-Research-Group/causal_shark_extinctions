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


# get adjustments sets for temperature, for the total effect
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
                  ndraws = 100) 

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


# save plot
ggsave(plot_shelf, filename = here("figures",
                                 "effect_shelf_area.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)
