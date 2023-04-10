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
                            "processed_merged_data.rds")) %>% 
  drop_na(order) %>% 
  mutate(superorder = fct_collapse(order,
                              Sharks = c("Carcharhiniformes",
                                         "Orectolobiformes",
                                         "Lamniformes",
                                         "Squaliformes",
                                         "Hexanchiformes",
                                         "Echinorhiniformes",
                                         "Heterodontiformes",
                                         "Pristiophoriformes",
                                         "Squatiniformes",
                                         "Synechodontiformes"),
                              Rays = c("Myliobatiformes",
                                       "Rajiformes",
                                       "Rhinopristiformes",
                                       "Torpediniformes"))) 


# total effect adjustment sets --------------------------------------------


# load the graph 
dag <- downloadGraph("dagitty.net/mjiV5Qf")


# get adjustments sets for the total effect
adjustmentSets(dag, 
               exposure = "taxonomic identity",
               outcome = "extinction risk", 
               effect = "total")

# no need to account for confounders


# fit models --------------------------------------------------------------


# fit final model
mod1 <- brm_logistic("ext_signal ~ order")

# same for superorders
mod2 <- brm_logistic("ext_signal ~ superorder")


# extract predictions ----------------------------------------------------------

# set up grid to average over
dat_pred_order <- dat_merged %>% 
  distinct(order, superorder) %>% 
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
  drop_na(order) %>% 
  mutate(order = fct_reorder(order, median_epred))  
  
  
# visualise
plot_order <- dat_pred_order %>%
  ggplot(aes(.epred, order)) +
  stat_histinterval(aes(fill = superorder), 
                    shape = 21, 
                    interval_alpha = 0,
                    slab_colour = "#9C8899",
                    slab_fill = "#9C8899",
                    point_size = 1.5, 
                    point_colour = "grey30") +
  geom_label(aes(label = n_lab, x = max_epred + 0.05),
             data = dat_pred_order %>%
               distinct(order, max_epred, n) %>%
               mutate(n_lab = paste0("n = ", n)),
             label.size = 0,
             size = 10/.pt,
             colour = colour_grey) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  scale_fill_manual(name = NULL, 
                    values = c("#EA8778", "#FFBE62"), 
                    limits = c("Sharks", "Rays")) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(y = NULL, 
       x = "Extinction Risk [%]") +
  theme(legend.position = c(0.8, 0.3), 
        legend.text = element_text(colour = "grey50", size = 10)) +
  guides(fill = guide_legend(override.aes = list(slab_fill = "white", 
                                                 slab_colour = "white")))

# ggsave(plot_order, filename = here("figures",
#                                    "presentations", 
#                                    "taxonomy_4.png"), 
#        width = image_width, height = image_height,
#        units = image_units, 
#        bg = "white", device = ragg::agg_png)

# same for superorders
# set up grid to average over
dat_pred_superorder <- dat_merged %>% 
  distinct(superorder) %>% 
  # add posterior draws
  add_epred_draws(mod2, 
                  ndraws = 10000) %>% 
  # arrange
  group_by(superorder) %>% 
  mutate(median_epred = median(.epred), 
         max_epred = max(.epred)) %>% 
  ungroup() %>% 
  drop_na(superorder) %>% 
  mutate(superorder = fct_reorder(superorder, median_epred))  

# visualise
plot_super <- dat_pred_superorder %>%
  filter(superorder != "incertae sedis") %>% 
  ggplot(aes(.epred, superorder)) +
  stat_histinterval(shape = 21, 
                    interval_alpha = 0,
                    slab_colour = "#705976", 
                    slab_fill = "#705976", 
                    fill = "white", 
                    point_size = 1.5, 
                    point_colour = "grey20") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100"), 
                     expand = c(0,0)) +
  scale_y_discrete(expand = c(0.3, 0), 
                   breaks = NULL) +
  coord_cartesian(xlim = c(0, 0.4)) +
  labs(y = NULL, 
       x = NULL) +
  theme(panel.border = element_rect(colour = "grey50",
                                    fill = NA, size = 0.5))


# save plot
ggsave(plot_order, filename = here("figures",
                                  "effect_taxonomy.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)
