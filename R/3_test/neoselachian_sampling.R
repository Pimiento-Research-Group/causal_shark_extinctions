# load packages
library(tidyverse)
library(here)
library(deeptime)
library(patchwork)

# plotting configurations
source(here("R", "config_file.R"))


# read data ---------------------------------------------------------------


# load occurrence database
dat_occurrences <- read_rds(here("data",
                                 "fossil_occurrences",
                                 "database_occurrences_15_Apr_2023.rds"))

# get PyRate output
dat_pyrate <- read_delim(here("data",
                              "fossil_occurrences",
                              "combined_10_se_est_species_names.txt")) %>% 
  # estimate fad and lad for each species
  mutate(across(c(ts, te), abs)) %>% 
  group_by(species) %>% 
  summarise(ori_age = mean(ts),
            ext_age = mean(te)) %>% 
  # get superorder
  mutate(accepted_name = str_replace(species, "_", " ")) %>% 
  left_join(dat_occurrences %>% 
              distinct(accepted_name, superorder)) %>% 
  drop_na(superorder)



# visualise ---------------------------------------------------------------

# ranges over time
plot_1 <- dat_pyrate %>%
  filter(superorder != "incertae sedis") %>% 
  mutate(species = fct_reorder(species, ext_age)) %>% 
  ggplot(aes(xmin = ori_age,
             xmax = ext_age, 
             y = species, 
             colour = superorder)) +
  geom_linerange() +
  scale_color_manual(name = NULL, 
                     values = c("#FFBE62", "#EA8778"), 
                     labels = c("Rays", "Sharks")) +
  coord_geo(xlim = c(0, 150), 
            dat = list("epochs", "periods"),
            pos = list("b", "b"),
            alpha = 0.2, 
            height = unit(0.8, "line"), 
            size = list(7/.pt, 10/.pt),
            lab_color = "grey20", 
            color = "grey20", 
            abbrv = list(TRUE, FALSE), 
            fill = "white",
            expand = TRUE, 
            lwd = list(0.4, 0.5)) +
  scale_x_reverse() +
  scale_y_discrete(expand = expansion(mult = c(0.29, 0))) +
  labs(y = "Taxon ranges", 
       x = "Age [myr]") +
  theme(axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        legend.position = c(0.8, 0.8))

# number of occurrences
plot_2 <- dat_occurrences %>%
  mutate(mid_age = ((Max_Ma - Min_Ma)/2) + Min_Ma) %>% 
  filter(superorder %in% c("Batoidea", 
                           "Galeomorphii")) %>% 
  ggplot() +
  geom_density(aes(mid_age, 
                   colour = superorder), 
               linewidth = 0.8) +
  labs(y = "Density", 
       x = NULL) +
  scale_color_manual(name = NULL, 
                     values = c("#FFBE62", "#EA8778")) +
  scale_x_reverse() +
  scale_y_continuous(position = "right", 
                     expand = expansion(mult = c(0, 0.01))) +
  coord_cartesian(xlim = c(150, 0)) +
  theme(legend.position = "none", 
        axis.text = element_blank(), 
        axis.ticks = element_blank()) 



# patch together ----------------------------------------------------------

# add as inset
plot_final <- plot_1 +
  inset_element(plot_2, -0.02, 0, 1.035, 0.26) +
  plot_annotation(tag_level = "a")

# save plot
ggsave(plot_final, filename = here("figures",
                                   "supplement",
                                   "neoselachian_occurrences.png"), 
       width = image_width, height = image_height*1.5,
       units = image_units, 
       bg = "white", device = ragg::agg_png) 
