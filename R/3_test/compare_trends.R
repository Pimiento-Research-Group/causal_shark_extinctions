# load packages
library(tidyverse)
library(here)
library(tidybayes)


# plotting configurations
source(here("R", "config_file.R"))

# get original trend estimates
dat_main <- list.files(here("data",
                            "predictions"), 
           full.names = TRUE) %>%
  str_subset("10myr", negate = TRUE) %>% 
  str_subset("genus", negate = TRUE) %>% 
  str_subset("trend") %>% 
  map_df(read_rds) %>% 
  # for some reason shelf area was entered as NA
  replace_na(list(coef_name = "b_shelf_area")) %>% 
  # join parameters together
  mutate(coef_name = as.factor(coef_name), 
         coef_name = fct_collapse(coef_name,
                                  Paleotemperature = grep(":", levels(coef_name), value = TRUE),
                                  Temperature = grep("binned", levels(coef_name), value = TRUE), 
                                  "Sea level" = grep("sea", levels(coef_name), value = TRUE), 
                                  "Shelf area" = grep("shelf", levels(coef_name), value = TRUE), 
                                  Productivity = grep("std", levels(coef_name), value = TRUE), 
         )) %>% 
  add_column(data_source = "main")
  

# get data from 10 myr resolution
dat_10myr <- list.files(here("data",
                             "predictions"), 
           full.names = TRUE) %>%
  str_subset("10myr") %>% 
  str_subset("trend") %>% 
  map_df(read_rds) %>% 
  # for some reason shelf area was entered as NA
  replace_na(list(coef_name = "b_shelf_area"))  %>% 
  # join parameters together
  mutate(coef_name = as.factor(coef_name), 
         coef_name = fct_collapse(coef_name, 
                                  Paleotemperature = grep(":", levels(coef_name), value = TRUE), 
                                  Productivity = grep("std", levels(coef_name), value = TRUE), 
                                  Temperature = grep("binned", levels(coef_name), value = TRUE), 
                                  "Sea level" = grep("sea", levels(coef_name), value = TRUE), 
                                  "Shelf area" = grep("shelf", levels(coef_name), value = TRUE)
         )) %>% 
  add_column(data_source = "10_myr")


# get data from genus resolution
dat_genus <- list.files(here("data",
                             "predictions"), 
                        full.names = TRUE) %>%
  str_subset("genus") %>% 
  str_subset("trend") %>% 
  map_df(read_rds) %>% 
  # for some reason shelf area was entered as NA
  replace_na(list(coef_name = "b_shelf_area"))  %>% 
  # join parameters together
  mutate(coef_name = as.factor(coef_name), 
         coef_name = fct_collapse(coef_name, 
                                  Abundance = c("b_n_genus", "b_abund"), 
                                  "Geographic range" = c("b_geo_dist_std"), 
                                  Paleotemperature = c("b_temp_deep_st:temp_deep_lt2", 
                                                       "b_temp_deep_st:temp_deep_lt4",
                                                       "b_temp_gat_st:temp_gat_lt1"), 
                                  Productivity = c("b_d13C_std", "b_sr_value_std"), 
                                  Temperature = c("b_temp_gat_binned"), 
                                  "Sea level" = c("b_sea_level"), 
                                  "Shelf area" = c("b_shelf_area")
         )) %>% 
  add_column(data_source = "Genus")

# merge
full_join(dat_main, dat_10myr) %>% 
  # full_join(dat_genus) %>% 
  ggplot(aes(y = data_source, coef_val, 
             colour = data_source)) +
  geom_vline(xintercept = 0) +
  stat_pointinterval() +
  facet_wrap(~ coef_name, 
             scales = "free") +
  scale_y_discrete(breaks = NULL) +
  scale_x_continuous(breaks = c(0)) +
  # scale_color_manual(name = NULL,
  #                    values = c("coral", "steelblue", "darkgreen"), 
  #                    labels = c("10 myr", 
  #                               "Genus", 
  #                               "Stages")) +
  labs(y = NULL, 
       x = NULL) +
  theme(panel.border = element_rect(linewidth = 1, 
                                    colour = "grey80", 
                                    fill = NA), 
        legend.position = c(0.8, 0.1))



