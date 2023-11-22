library(tidyverse)
library(here)
library(ggsvg)
library(patchwork)

# plotting configurations
source(here("R", "config_file.R"))


# read data ---------------------------------------------------------------

# logit per order
dat_order <- read_rds(here(here("data", 
                      "logits", 
                      "logit_order.rds")))

# metabolism test on family level
dat_family <- read_rds(here(here("data",
                                 "metabolism",
                                 "metab_family.rds"))) %>% 
  # summarise
  group_by(family) %>% 
  ggdist::median_qi(value) %>% 
  # negate for visualisation
  mutate(across(2:4, ~ -.x))


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


# visualise order --------------------------------------------------------

# batoidea
svg_bato <- "https://images.phylopic.org/images/465ee9b8-3a2d-4c3d-bf43-6301f50d1f2e/vector.svg" %>%
  readLines() %>% 
  str_replace("#000000", "#FFBE62") %>% 
  paste(collapse = "\n")  

# galeomorphii
svg_galeo <- "https://images.phylopic.org/images/42135d61-3549-45d2-841c-4147548b0fad/vector.svg" %>%
  readLines() %>% 
  str_replace("#000000", "#169199") %>% 
  paste(collapse = "\n") 

# squalomorphii
svg_squalo <- "https://images.phylopic.org/images/60e7f957-1137-48db-8e1d-04b1c41c0c18/vector.svg" %>%
  readLines() %>% 
  str_replace("#000000", "#BD7CD5") %>% 
  paste(collapse = "\n") 


plot_order <- dat_order %>%
  ggplot(aes(y = order, 
             x = logit_val)) +
  geom_vline(xintercept = 0, 
             linetype = "dotted", 
             colour = colour_grey) +
  geom_linerange(aes(xmin = .lower,
                     xmax = .upper),
                 linewidth = 0.4,
                 colour = "white") +
  geom_point_svg(aes(size = n),
                 svg = svg_bato,
                 data = filter(dat_order, 
                               superorder == "Batoidea")) +
  geom_point_svg(aes(size = n),
                 svg = svg_galeo,
                 data = filter(dat_order, 
                               superorder == "Galeomorphii")) +
  geom_point_svg(aes(size = n),
                 svg = svg_squalo,
                 data = filter(dat_order, 
                               superorder == "Squalomorphii")) +
  geom_text(aes(x = logit_val,
                y = order,
                label = order),
            position = position_nudge(x = -0.35),
            colour = "grey20",
            size = 10/.pt,
            hjust = "right") +
  annotate("label",
           y = 15.5, 
           x = c(-1.1, -2.1, -3.1), 
           colour = c("#FFBE62", 
                      "#169199", 
                      "#BD7CD5"), 
           label.size = 0, 
           size = 11/.pt, 
           label = c("Batoidea", 
                     "Galeomorphii", 
                     "Squalomorphii"), 
           hjust = 0) +
  scale_size_continuous(range = c(5, 17)) +
  scale_x_continuous("Temperature dependancy [log-odds]", 
                     expand = expansion(mult = c(0.1, 0)), 
                     breaks = c(-4, -2, 0),
                     limits = c(-4.3, 1)) +
  scale_y_discrete(NULL, 
                   expand = expansion(mult = c(0.1, 0.2))) +
  guides(size = "none", 
         fill = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = "none", 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())


# save plot
ggsave(plot_order, filename = here("figures",
                                   "logit_per_group.png"),
       width = image_width, height = image_height,
       units = image_units,
       bg = "white", device = ragg::agg_png)
