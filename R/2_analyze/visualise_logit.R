library(tidyverse)
library(here)

# plotting configurations
source(here("R", "config_file.R"))


# read data ---------------------------------------------------------------

# logit per order
dat_order <- read_rds(here(here("data", 
                      "logits", 
                      "logit_order.rds")))


# logit per order ---------------------------------------------------------


# reformat
dat_order <- dat_order %>%
  filter(order != "incertae sedis") %>% 
  # add superorder
  left_join(read_rds(here("data",
                          "fossil_occurrences",
                          "database_occurrences_01_Nov_2023.rds")) %>% 
              drop_na(order, superorder) %>% 
              filter(order != "incertae sedis") %>% 
              count(superorder, order)) %>% 
  # abbreviate order for nicer plotting
  mutate(order = str_replace_all(order, "formes", "."), 
         order = fct_reorder(order, logit_val))

# visualise
# plot_order <- 
dat_order %>%
  ggplot(aes(y = order, 
             x = logit_val)) +
  geom_vline(xintercept = 0, 
             linetype = "dotted", 
             colour = colour_grey) +
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
           x = -1.5, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "Low susceptibility") +
  annotate("text",
           y = -1, 
           x = -6.5, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "High susceptibility") +
  annotate("curve",
           y = -1, 
           yend = -1, 
           x = -5, 
           xend = -3,
           colour = "grey50", 
           curvature = 0,
           arrow = arrow(length = unit(.2,"cm"), 
                         ends = "both")) +
  labs(y = NULL, 
       x = "Temperature dependancy [log-odds]") +
  scale_fill_manual(name = NULL,
                    values = c("#FFBE62", "#C8A5D9", "#BD8D9E")) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0)), 
                     breaks = c(-5, -2.5, 0),
                     limits = c(-7, 3)) +
  scale_y_discrete(expand = expansion(mult = c(0.25, 0.1))) +
  guides(size = "none", 
         fill = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = c(0.88, 0.6), 
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
