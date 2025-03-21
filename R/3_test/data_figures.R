library(here)
library(tidyverse)
library(patchwork)
library(ggdist)

# plotting configurations
source(here("R", "config_file.R"))


# stage level -------------------------------------------------------------

# load 95 bin Phanerozoic time scale based on the stratigraphic stages 
# of Gradstein et al. 2020
load(here("data", 
          "stages.Rdata"))

# fossil data and environmental proxy data on species level
dat_stages <- read_rds(here("data",
                            "processed_merged_data.rds"))

plot_stages <- dat_stages %>%
  distinct(
    bin,
    cont_area,
    temp_gat_binned,
    temp_deep_binned,
    d13C_std,
    sr_mean_std,
    sea_level
  ) %>% 
  rename("F) Shelf area [%]" = cont_area, 
         "D) d13C [std]" = d13C_std, 
         "C) Sea level [m]" = sea_level, 
         "E) 87Sr/ 86Sr [std]" = sr_mean_std, 
         "A) Temperature deep [°C]" = temp_deep_binned, 
         "B) Temperature gat [°C]" = temp_gat_binned) %>% 
  pivot_longer(-bin) %>% 
  ggplot(aes(bin, value)) +
  geom_line(colour = "grey20") +
  scale_x_continuous(breaks = c(68, 76, 85, 95), 
                     labels = c(150, 100, 50, 0), 
                     name = "Age [myr]") +
  facet_wrap(~ name, 
             scales = "free") +
  labs(y = NULL) +
  theme(strip.text = element_text(hjust = 0))


# save plot
ggsave(plot_stages, filename = here("figures",
                                    "supplement",
                                    "environmental_data_stages.png"), 
       width = image_width, height = image_height*1.5,
       units = image_units, 
       bg = "white", device = ragg::agg_png)  


# 1myr level -------------------------------------------------------------


# fossil data and environmental proxy data on species level
dat_1myr <- read_rds(here("data",
                          "processed_fossil_data_cenozoic.rds"))

plot_1myr <- dat_1myr %>%
  distinct(
    bin,
    cont_area,
    temp_binned,
    d13C_std,
    diatom_rich_std,
    sea_level
  ) %>% 
  rename("E) Shelf area [%]" = cont_area, 
         "C) d13C [std]" = d13C_std, 
         "B) Sea level [m]" = sea_level,
         "D) Diatom richness [std]" = diatom_rich_std, 
         "B) Global temperature [°C]" = temp_binned) %>% 
  pivot_longer(-bin) %>% 
  ggplot(aes(bin, value)) +
  geom_line(colour = "grey20") +
  scale_x_reverse(name = "Age [myr]") +
  labs(y = NULL) +
  facet_wrap(~ name, 
             scales = "free") +
  theme(strip.text = element_text(hjust = 0))



# save plot
ggsave(plot_1myr, filename = here("figures",
                                  "supplement",
                                  "environmental_data_1myr.png"), 
       width = image_width, height = image_height*1.5,
       units = image_units, 
       bg = "white", device = ragg::agg_png)  


# sampling bias -----------------------------------------------------------

plot_sampling <- dat_stages %>%
  distinct(bin,
           n_units_std,
           shark_collections_std, 
           pbdb_collections_std, 
           outcrop_area_std) %>% 
  left_join(stages %>% 
              select(bin = stg, 
                     age = mid)) %>% 
  select(-bin) %>% 
  rename("C) Marine rock units [std]" = n_units_std, 
         "B) Neoselachian collections [std]" = shark_collections_std, 
         "A) PBDB collections [std]" = pbdb_collections_std,
         "D) Marine outcrop area [std]" = outcrop_area_std) %>% 
  pivot_longer(-age) %>%
  add_column(temp_scale = "Stages") %>% 
  bind_rows(dat_1myr %>% 
              distinct(age = bin, 
                       n_units_std, 
                       shark_collections_std, 
                       pbdb_collections_std) %>% 
              rename("C) Marine rock units [std]" = n_units_std, 
                     "B) Neoselachian collections [std]" = shark_collections_std, 
                     "A) PBDB collections [std]" = pbdb_collections_std) %>%  
              pivot_longer(-age) %>% 
              add_column(temp_scale = "1myr") %>% 
              filter(between(age, 1, 55))) %>% 
  ggplot(aes(age, value, 
             colour = temp_scale)) +
  geom_line() +
  scale_color_manual(values = c("grey10", 
                                "grey40"), 
                     name = NULL) +
  scale_x_reverse(name = "Age [myr]") +
  labs(y = NULL) +
  facet_wrap(~ name, 
             scales = "free") +
  theme(strip.text = element_text(hjust = 0), 
        legend.position = "bottom")



# save plot
ggsave(plot_sampling, filename = here("figures",
                                      "supplement",
                                      "sampling_data.png"), 
       width = image_width, height = image_height*1.5,
       units = image_units, 
       bg = "white", device = ragg::agg_png)



# number of occurrences ---------------------------------------------------

# read in database
dat_occ <- read_rds(here("data", 
                         "fossil_occurrences",
                         "database_occurrences_15_Apr_2023.rds"))

# summarise
plot_occ <- dat_occ %>% 
  count(rank, status) %>% 
  filter(rank %in% c("species", "genus")) %>% 
  filter(status != "NA") %>% 
  mutate(across(where(is.character), str_to_sentence)) %>% 
  ggplot(aes(rank, n, 
             fill = status)) +
  geom_col() +
  scale_fill_manual(values = c(colour_yellow, colour_purple)) +
  labs(x = NULL, 
       y = "# Occurrences", 
       fill = NULL) +
  theme(legend.position = "top")

# same for latitudinal coverage
plot_lat <- dat_occ %>% 
  filter(status != "NA", 
         rank %in% c("species", "genus")) %>% 
  ggplot() +
  geom_linerange(aes(y = paleolat, 
                     xmin = Min_Ma, 
                     xmax = Max_Ma,
                     colour = status), 
                 alpha = 0.6) +
  scale_colour_manual(values = c(colour_yellow, colour_purple)) +
  scale_x_reverse() +
  labs(x = "Age [myr]", 
       y = "Paleolatitude") +
  theme(legend.position = "none") 

# patch together
plot_occ_final <- plot_lat / 
  plot_occ +
  plot_annotation(tag_levels = "A")
  

# save plot
ggsave(plot_occ_final, filename = here("figures",
                                       "supplement",
                                       "nr_occurrences.png"),
       width = image_width, height = image_height*1.5,
       units = image_units,
       bg = "white", device = ragg::agg_png)
