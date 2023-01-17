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
  mean_qi(value) 

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

# save plot
ggsave(plot_sea, filename = here("figures",
                                 "effect_sea_level.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)
