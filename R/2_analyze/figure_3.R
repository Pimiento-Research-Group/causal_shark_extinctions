library(tidyverse)
library(here)
library(ggsvg)
library(patchwork)
library(ggdist)
library(ggforce)

# plottiggforce# plotting configurations
source(here("R", "config_file.R"))

# read data ---------------------------------------------------------------

# family level logit values (temperature dependancy)
dat_family <- read_rds(here(here("data",
                                 "logits", 
                                 "logit_family_full.rds"))) %>% 
  group_by(family) %>% 
  median_qi(value) %>% 
  # add superorder
  left_join(read_rds(here("data",
                          "fossil_occurrences",
                          "database_occurrences_01_Nov_2023.rds")) %>% 
              drop_na(family, superorder) %>% 
              filter(order != "incertae sedis") %>% 
              count(superorder, family))
  



# visualise family --------------------------------------------------------

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


plot_fam <- dat_family %>%
  ggplot(aes(y = family, 
             x = value)) +
  geom_vline(xintercept = 0, 
             linetype = "dotted", 
             colour = colour_grey) +
  geom_point_svg(aes(size = 20), 
                 size = 10, 
                 svg = svg_bato,
                 data = filter(dat_family, 
                               superorder == "Batoidea")) +
  geom_point_svg(aes(size = 20),
                 size = 10,
                 svg = svg_galeo,
                 data = filter(dat_family, 
                               superorder == "Galeomorphii")) +
  geom_point_svg(aes(size = 20), 
                 size = 10,
                 svg = svg_squalo,
                 data = filter(dat_family, 
                               superorder == "Squalomorphii")) +
  geom_mark_ellipse(aes(label = family), 
                    label.fontsize = 10, 
                    label.fill = alpha("white", 0.5),  
                    label.colour = "grey40",
                    colour = "grey40",
                    con.colour = "grey40",
                    expand = unit(5, "mm"), 
                    con.cap = 0, 
                    data = dat_family %>% 
                      filter(family %in% c("Alopiidae", 
                                           "Lamnidae", 
                                           "Otodontidae"))) +
  annotate("label",
           y = 83, 
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
  scale_x_continuous("Temperature dependancy [log-odds]", 
                     expand = expansion(mult = c(0.1, 0)), 
                     limits = c(-4.5, 1),
                     breaks = c(-4, -2, 0)) +
  scale_y_discrete(NULL, 
                   expand = expansion(mult = c(0.1, 0.2))) +
  guides(size = "none", 
         fill = guide_legend(override.aes = list(size = 2))) +
  theme(legend.position = "none", 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())


# save plot
ggsave(plot_fam, filename = here("figures",
                                   "figure_3.png"),
       width = image_width, height = image_height,
       units = image_units,
       bg = "white", device = ragg::agg_png)
