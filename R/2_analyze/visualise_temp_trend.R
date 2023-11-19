# load packages
library(tidyverse)
library(here)
library(ggdist)
library(patchwork)
library(gganimate)

# plotting configurations
source(here("R", "config_file.R"))



# load data  ------------------------------------------------------------

# get predictions from different scales/ data sets
dat_pred <- read_rds(here("data", 
                          "temp_effect.rds"))


# add counter for gif
dat_pred_cl <- dat_pred %>%
  mutate(dataset = case_when(
    dataset == "ceno" ~ "1myr", 
    dataset == "genus" ~ "Genus", 
    dataset == "stage" ~ "Stages"
  ), 
  dataset = factor(dataset, 
                   levels = c("Stages", 
                              "Genus", 
                              "1myr"))) %>% 
  mutate(nr_draw = as.integer(nr_draw), 
         counter = interaction(nr_draw, dataset), 
         counter = as.integer(counter))


# visualise ---------------------------------------------------------------

# empty plot
plot_temp_empty <- dat_pred_cl %>%
  ggplot(aes(temperature,
             value,
             group = counter,
             colour = dataset)) +
  geom_line(alpha = 0, 
            key_glyph = "rect") +
  scale_colour_manual(
    values = c("#4C634C",
               colour_coral, 
               colour_purple),
    labels = c("Species - Stages", 
               "Genera - Stages", 
               "Species - Cenozoic subset"),
    name = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.9, 
                                                   linewidth = 0.6))) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  coord_cartesian(ylim = c(0, 1), 
                  xlim = c(0, 25)) +
  labs(y = "Extinction Risk [%]", 
       x = "Global Temperature [°C]") +
  theme(legend.position = c(0.75, 0.75), 
        legend.key.width = unit(10, "mm"), 
        legend.key.height = unit(8, "mm"))

# save plot
ggsave(plot_temp_empty,
       filename = here("figures",
                       "temperature_empty.png"), 
       width = image_width, 
       height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)

# gif ---------------------------------------------------------------------

# raw plot
plot_temp <- dat_pred_cl %>%
  ggplot(aes(temperature,
             value,
             group = counter,
             colour = dataset)) +
  geom_line(alpha = 0.3, 
            key_glyph = "rect") +
  scale_colour_manual(
    values = c("#4C634C",
               colour_coral, 
               colour_purple),
    labels = c("Species - Stages", 
               "Genera - Stages", 
               "Species - Cenozoic subset"),
    name = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.9, 
                                                   linewidth = 0.6))) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  coord_cartesian(ylim = c(0, 1), 
                  xlim = c(0, 25)) +
  labs(y = "Extinction Risk [%]", 
       x = "Global Temperature [°C]") +
  theme(legend.position = c(0.75, 0.75), 
        legend.key.width = unit(10, "mm"), 
        legend.key.height = unit(8, "mm"))


# add a transition
gif_temp <- plot_temp + 
  transition_manual(counter, 
                    cumulative = TRUE)

# save animation
anim_save(animation = gif_temp, 
          filename = "gif_temperature.gif",
          path = here("figures"),
          nframes = 250, 
          res = 200, 
          height = image_height, 
          width = image_width, 
          height = image_height,
          units = image_units, 
          renderer = gifski_renderer(loop = FALSE))
