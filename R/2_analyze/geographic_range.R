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
               exposure = "geographic range",
               outcome = "extinction risk", 
               effect = "total")

# two sets of adjustments 
# { abundance, taxonomic identity }
# { outcrop area, paleotemperature, preservation potential, productivity, 
# sampling effort, shelf area, taxonomic identity, temperature }


# fit models --------------------------------------------------------------


# first adjustment set
mod1 <- brm_logistic("ext_signal ~ geo_dist_std + abund + order")
mod2 <- brm_logistic("ext_signal ~ range_lat_std + abund + order")

# number of genera
mod3 <- brm_logistic("ext_signal ~ geo_dist_std + n_genus + order")
mod4 <- brm_logistic("ext_signal ~ range_lat_std + n_genus + order")


# average models ----------------------------------------------------------

# set up grid to average over
dat_new <- tibble(geo_dist_std = -2:2, 
                  n_genus = geo_dist_std) %>%
  # average over taxonomy
  expand_grid(order = unique(dat_merged$order)) %>% 
  mutate(range_lat_std = geo_dist_std) %>% 
  # keep adjustment at average
  add_column(abund = mean(dat_merged$abund))

# set up number of draws 
nr_draws <- 100

# average prediction by model stacking
dat_pred <- pp_average(mod1, mod2,
                       mod3, mod4,
                       newdata = dat_new,
                       seed = 1708,
                       summary = FALSE, 
                       method = "posterior_epred", 
                       ndraws = nr_draws) %>% 
  as_tibble() %>% 
  mutate(nr_draw = rownames(.)) %>% 
  pivot_longer(cols = contains("V")) %>% 
  add_column(geo_range = rep(dat_new$geo_dist_std, nr_draws)) %>% 
  group_by(nr_draw, geo_range) %>%
  mean_qi(value) 

# average over posterior draws
dat_pred_av <- dat_pred %>% 
  group_by(geo_range) %>% 
  summarise(value = mean(value))

# visualise
plot_range <- dat_pred %>%
  ggplot(aes(geo_range, value)) +
  stat_dots(aes(y = as.numeric(ext_signal == 1), 
                side = ifelse(ext_signal == 1, "bottom", "top")),
            slab_colour = colour_grey,
            slab_size = 0.5,
            slab_alpha = 0.3,
            scale = 0.2,
            data = dat_merged %>% 
              pivot_longer(cols = c(geo_dist_std, range_lat_std), 
                           values_to = "geo_range", 
                           names_to = "geo_range_name") %>% 
              # spread out a bit
              mutate(geo_range = geo_range + rnorm(nrow(.), 0, 0.1))) +
  geom_line(aes(group = nr_draw), 
            alpha = 0.2, color = colour_purple) +
  geom_line(colour = "white",
            linewidth = 1.3, 
            data = dat_pred_av) +
  geom_line(colour = colour_purple,
            linewidth = 1, 
            data = dat_pred_av) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  coord_cartesian(xlim = c(-2, 2)) +
  labs(y = "Extinction Risk [%]", 
       x = "Geographic Range [srd]")

# save plot
ggsave(plot_range, filename = here("figures",
                                   "effect_geographic_range.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)
