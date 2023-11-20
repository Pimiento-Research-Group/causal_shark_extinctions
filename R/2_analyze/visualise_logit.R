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
  str_replace("#000000", "#C8A5D9") %>% 
  paste(collapse = "\n") 

# squalomorphii
svg_squalo <- "https://images.phylopic.org/images/60e7f957-1137-48db-8e1d-04b1c41c0c18/vector.svg" %>%
  readLines() %>% 
  str_replace("#000000", "#BD7CD5") %>% 
  paste(collapse = "\n") 


# visualise
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
            position = position_nudge(x = -0.7),
            size = 10/.pt,
            colour = "grey20", 
            hjust = "right") +
  annotate("text",
           y = c(13.5, 12, 10.5), 
           x = 0.5, 
           colour = c("#FFBE62", 
                      "#C8A5D9", 
                      "#BD7CD5"), 
           fontface = "bold",
           size = 11/.pt, 
           label = c("Batoidea", 
                     "Galeomorphii", 
                     "Squalomorphii"), 
           hjust = 0) +
  scale_size_continuous(range = c(5, 16)) +
  scale_x_continuous("Temperature dependancy [log-odds]", 
                     expand = expansion(mult = c(0.1, 0)), 
                     breaks = c(-4, -2, 0),
                     limits = c(-4.95, 3)) +
  scale_y_discrete(NULL, 
                   expand = expansion(mult = c(0.1, 0.1))) +
  guides(size = "none", 
         fill = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = c(0.88, 0.6), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())



# percentage change -------------------------------------------------------

# load in svgs

# alopiidae
svg_alop <- "https://images.phylopic.org/images/b65312ae-91b5-45b1-9553-c192f1000aba/vector.svg" %>%
  readLines() %>% 
  str_replace("#000000", "#282828") %>% 
  paste(collapse = "\n") 

# lamnidae
svg_lamn <- "https://images.phylopic.org/images/545d45f0-0dd1-4cfd-aad6-2b835223ea0d/vector.svg" %>%
  readLines() %>% 
  str_replace("#000000", "#282828") %>% 
  paste(collapse = "\n") 

# otodontidae
svg_oto <- "https://images.phylopic.org/images/4c45a2bc-81a2-4b4a-a91a-3decdec14499/vector.svg" %>%
  readLines() %>% 
  str_replace("#000000", "#282828") %>% 
  paste(collapse = "\n") 


# visualise
plot_family <- dat_family %>%
  ggplot(aes(x = value, y = family)) +
  geom_vline(xintercept = 0, 
             linetype = "dotted", 
             colour = colour_grey) +
  geom_point_svg(aes(x = -0.21), 
                 size = 36,
                 svg = svg_alop,
                 data = filter(dat_family,
                               family == "alo")) +
  geom_point_svg(aes(x = -0.178), 
                 size = 42,
                 svg = svg_lamn,
                 data = filter(dat_family,
                               family == "lamni")) +
  geom_point_svg(aes(x = 0.01), 
                 size = 69,
                 svg = svg_oto,
                 data = filter(dat_family,
                               family == "oto")) +
  geom_text(aes(label = round(value*100, 0) %>% paste0("%")), 
            position = position_nudge(y = 0.03),
            colour = "white", 
            size = 12/.pt) +
  annotate("text",
           y = c(2.5, 2, 1), 
           x = c(0.04, -0.08, -0.13), 
           colour = "#282828",
           size = 12/.pt, 
           label = c("Otodontidae", 
                     "Lamnidae", 
                     "Alopiidae")) +
  scale_x_continuous(name = "Difference in temperature dependancy [%]", 
                     limits = c(-0.3, 0.15)) +
  scale_y_discrete(name = NULL) +
  expand_limits(y = c(1.5, 4)) +
  theme(axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

# patch together ----------------------------------------------------------

plot_full <- plot_order / 
  plot_family +
  plot_annotation(tag_levels = "a")

# save plot
ggsave(plot_full, filename = here("figures",
                                  "logit_per_group.png"),
       width = image_width, height = image_height*1.9,
       units = image_units,
       bg = "white", device = ragg::agg_png)
