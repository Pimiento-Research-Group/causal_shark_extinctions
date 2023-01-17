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


# get adjustments sets for temperature, for the total effect
adjustmentSets(dag, 
               exposure = "taxonomic identity",
               outcome = "extinction risk", 
               effect = "total")

# no need to account for confounders


# fit models --------------------------------------------------------------


# fit final model
mod1 <- brm_logistic("ext_signal ~ order")

conditional_effects(mod1)



# extract predictions ----------------------------------------------------------

# set up grid to average over
dat_pred <- dat_merged %>% 
  distinct(order) %>% 
  # add posterior draws
  add_epred_draws(mod1, 
                  ndraws = 10000) %>% 
  # add number of samples
  full_join(count(dat_merged, order)) %>% 
  # arrange
  group_by(order) %>% 
  mutate(median_epred = median(.epred), 
         max_epred = max(.epred)) %>% 
  ungroup() %>% 
  mutate(order = fct_reorder(order, median_epred))  
  
# visualise
plot_tax <- dat_pred %>% 
  ggplot(aes(.epred, order)) +
  stat_histinterval(shape = 21, 
                    interval_alpha = 0,
                    slab_colour = colour_purple, 
                    slab_fill = colour_purple, 
                    fill = "white") +
  geom_label(aes(label = n_lab, x = max_epred + 0.03), 
             data = dat_pred %>% 
               distinct(order, max_epred, n) %>% 
               mutate(n_lab = paste0("n = ", n)), 
             label.size = 0, 
             size = 10/.pt, 
             colour = colour_grey) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(y = NULL, 
       x = "Extinction Risk [%]")


# save plot
ggsave(plot_tax, filename = here("figures",
                                  "effect_taxonomy.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)
