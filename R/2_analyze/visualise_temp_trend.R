# load packages
library(tidyverse)
library(here)
library(ggdist)
library(patchwork)


# plotting configurations
source(here("R", "config_file.R"))



# load data  ------------------------------------------------------------

# get predictions from different scales/ data sets
set.seed(123)
dat_pred <- paste0("pred_temperature_", 
       c("ceno", 
         "genus", 
         "stage"), 
         ".rds") %>% 
  map_df(~ read_rds(here("data",
                         "predictions",
                         .x)) %>% 
           group_by(temperature) %>% 
           slice_sample(n = 100) %>% 
           ungroup() %>% 
           add_column(dataset = .x) %>% 
           mutate(dataset = str_remove_all(dataset, 
                                           "pred_temperature_"), 
                  dataset = str_remove_all(dataset, 
                                           ".rds")))


# visualise ---------------------------------------------------------------

plot_temp <- dat_pred %>%
  mutate(dataset = case_when(
    dataset == "ceno" ~ "1myr", 
    dataset == "genus" ~ "Genus", 
    dataset == "stage" ~ "Stages"
  ), 
  dataset = factor(dataset, 
                   levels = c("Stages", 
                              "Genus", 
                              "1myr"))) %>% 
  ggplot(aes(temperature,
             value,
             group = interaction(nr_draw, dataset),
             colour = dataset)) +
  geom_line(alpha = 0.2) +
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
  theme(legend.position = c(0.75, 0.75))



# logit per order ---------------------------------------------------------

# read in data
dat_order <- read_rds(here(here("data",
                                "logits",
                                "logit_order.rds"))) 
  

# reformat
dat_order <- dat_order %>%
  filter(order != "incertae sedis") %>% 
  group_by(scale) %>%
  mutate(rank_val = rank(value)) %>%
  group_by(order) %>%
  mean_qi(rank_val) %>% 
  # add superorder
  left_join(read_rds(here("data",
                          "fossil_occurrences",
                          "database_occurrences_15_Apr_2023.rds")) %>% 
              count(order, superorder)) %>% 
  filter(superorder %in% c("Batoidea",
                           "Galeomorphii")) %>%
  # abbreviate order for nicer plotting
  mutate(order = str_replace_all(order, "formes", "."), 
         order = fct_reorder(order, rank_val))

# visualise
plot_order <- dat_order %>%
  ggplot(aes(y = order, 
             x = rank_val)) +
  geom_linerange(aes(xmin = .lower, 
                     xmax = .upper), 
                 linewidth = 0.4, 
                 colour = "grey70") +
  geom_point(aes(fill = superorder, 
                 size = n), 
             shape = 21, 
             colour = "grey20") +
  geom_text(aes(x = .lower,
                y = order,
                label = order),
            position = position_nudge(x = -0.9),
            size = 10/.pt,
            colour = "grey20") +
  annotate("text",
           y = -1, 
           x = 2, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "Low susceptibility") +
  annotate("text",
           y = -1, 
           x = 12, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "High susceptibility") +
  annotate("curve",
           y = -1, 
           yend = -1, 
           x = 4, 
           xend = 10,
           colour = "grey50", 
           curvature = 0,
           arrow = arrow(length = unit(.2,"cm"), 
                         ends = "both")) +
  labs(y = NULL, 
       x = "Temperature dependancy\n[ranked]") +
  scale_fill_manual(name = NULL, 
                    values = c("#EA8778", "#FFBE62"), 
                    limits = c("Galeomorphii", "Batoidea")) +
  scale_x_continuous(breaks = c(1, 5, 10, 15), 
                     expand = expansion(mult = c(0.1, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0.25, 0.1))) +
  guides(size = "none", 
         fill = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = c(0.88, 0.3), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank(), 
        legend.text = element_text(colour = "grey20", size = 9))




# patch together ----------------------------------------------------------

plot_full <- plot_temp / 
  plot_order +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "a")

# save plot
ggsave(plot_full, filename = here("figures",
                                  "effect_temperature.png"),
       width = image_width, height = image_height*1.7,
       units = image_units,
       bg = "white", device = ragg::agg_png)
