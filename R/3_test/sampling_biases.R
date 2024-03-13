library(here)
library(tidyverse)
library(maps)

# load 95 bin Phanerozoic time scale based on the stratigraphic stages 
# of Gradstein et al. 2020
load(here("data", 
          "stages.Rdata"))



# outcrop area ------------------------------------------------------------


# fossil data and environmental proxy data on species level
dat_stages <- read_rds(here("data",
                            "processed_merged_data.rds"))


plot_sampling <- dat_stages %>%
  distinct(bin,
           n_units_std,
           shark_collections_std, 
           pbdb_collections_std, 
           outcrop_area_std) %>% 
  left_join(stages %>% 
              select(bin = stg, 
                     age = mid)) %>% 
  select(-c(bin, contains("collections"))) %>% 
  pivot_longer(-age) %>%
  ggplot(aes(age, value, 
             colour = name)) +
  geom_line(linewidth = 0.7) +
  annotate("label", 
           x = c(35, 20), 
           y = c(2.2, 3.1),
           label.size = 0,
           label = c("Rock Units (Macrostrat)", 
                     "Outcrop Area"), 
           size = 12/.pt,
           colour = c("coral4", 
                      "purple4")) +
  scale_color_manual(values = c("coral4", 
                                "purple4"), 
                     name = NULL) +
  scale_x_reverse(name = "Age [myr]") +
  labs(y = "Z-score", 
       title = "Sampling bias", 
       subtitle = "The area that can be sampled for fossils\nvaries widely through time") +
  theme_minimal() +
  theme(legend.position = "none") 


# save plot
ggsave(plot_sampling, filename = here("figures",
                                      "supplement",
                                      "sampling_area_bias.png"),
       width = 183, height = 100,
       units = "mm",
       bg = "white", device = ragg::agg_png)


# per continent and period ------------------------------------------------

dat_cont <- dat_occurrences %>% 
  count(early_period , Continent) %>%
  filter(Continent != "Ocean") %>% 
  filter(early_period != "NA") 

plot_cont <- dat_cont %>%
  mutate(Continent = ordered(Continent, 
                             levels = c("Europe", "North America", 
                                        "Asia", "Africa", 
                                        "South America", "Oceania", 
                                        "Antarctica")), 
         early_period = ordered(early_period, 
                                levels = c("Jurassic", 
                                           "Cretaceous", 
                                           "Paleogene", 
                                           "Neogene", 
                                           "Quaternary"))) %>% 
  ggplot() +
  geom_bar(aes(x = "", y = n, 
               fill = early_period), 
           width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  labs(fill = "Collections", 
       title = "Spatial bias") +
  facet_wrap(~ Continent) +
  theme_void() +
  theme(legend.position = c(0.6, 0.23), 
        legend.key.size = unit(3, "mm"))

# save plot
ggsave(plot_cont, filename = here("figures",
                                      "supplement",
                                      "spatial_sampling_bias.png"),
       width = 183, height = 100,
       units = "mm",
       bg = "white", device = ragg::agg_png)


# occurrences per taxon ---------------------------------------------------

# load occurrence database
dat_occurrences <- read_rds(here("data",
                                 "fossil_occurrences",
                                 "database_occurrences_15_Apr_2023.rds"))


plot_occ <- ggplot() +
  geom_col(aes(rowid, n),
           data = dat_occurrences %>% 
             filter(rank == "species") %>% 
             count(modified_identified_name) %>%
             arrange(desc(n)) %>%
             rowid_to_column(),
           colour = "white",
           fill = "purple4",
           linewidth = 1e-20,
           alpha = 0.8) +
  geom_col(aes(rowid, n),
           data = count(dat_occurrences, genus) %>%
             arrange(desc(n)) %>%
             rowid_to_column(),
           fill = "coral4",
           linewidth = 1e-20,
           alpha = 0.8) +
  annotate("label", 
           x = c(200, 400), 
           y = c(900, 30),
           label.size = 0,
           label = c("Genera", 
                     "Species"), 
           size = 12/.pt,
           colour = c("coral4", 
                      "purple4")) +
  scale_y_continuous(trans = scales::pseudo_log_trans(),
                      breaks = c(1, 10, 100, 1000)) +
  labs(x = NULL, 
       y = NULL, 
       title = "Taxonomic bias", 
       subtitle = "The distribution of occurrences per taxon is heavily skewed") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        panel.grid.major.x = element_blank())

# save plot
ggsave(plot_occ, filename = here("figures",
                                  "supplement",
                                  "taxonomic_bias.png"),
       width = 183, height = 100,
       units = "mm",
       bg = "white", device = ragg::agg_png)

