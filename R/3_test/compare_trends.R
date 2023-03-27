# load packages
library(tidyverse)
library(here)
library(tidybayes)


# plotting configurations
source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------


# get original trend estimates
dat_main <- list.files(here("data",
                            "predictions"), 
           full.names = TRUE) %>%
  str_subset("10myr", negate = TRUE) %>% 
  str_subset("genus", negate = TRUE) %>% 
  str_subset("cenozoic", negate = TRUE) %>% 
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
  add_column(data_source = "Stages")
  

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
  add_column(data_source = "10myr")


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
                                  Paleotemperature = grep(":", levels(coef_name), value = TRUE), 
                                  Productivity = grep("std", levels(coef_name), value = TRUE), 
                                  Temperature = grep("binned", levels(coef_name), value = TRUE), 
                                  "Sea level" = grep("sea", levels(coef_name), value = TRUE), 
                                  "Shelf area" = grep("shelf", levels(coef_name), value = TRUE)
         )) %>% 
  add_column(data_source = "Genus")


# get data from cenozoic 1myr resolution
dat_ceno <- list.files(here("data",
                            "predictions"), 
                        full.names = TRUE) %>%
  str_subset("cenozoic") %>% 
  str_subset("trend") %>% 
  map_df(read_rds) %>% 
  # join parameters together
  mutate(coef_name = as.factor(coef_name), 
         coef_name = fct_collapse(coef_name, 
                                  Paleotemperature = grep(":", levels(coef_name), value = TRUE), 
                                  Productivity = grep("std", levels(coef_name), value = TRUE), 
                                  Temperature = grep("binned", levels(coef_name), value = TRUE), 
                                  "Sea level" = grep("sea", levels(coef_name), value = TRUE), 
                                  "Shelf area" = grep("cont", levels(coef_name), value = TRUE)
         )) %>% 
  add_column(data_source = "Cenozoic/ 1myr")



# visualise ---------------------------------------------------------------


# create plot
plot_beta <- dat_main %>% 
  # merge
  full_join(dat_10myr) %>% 
  full_join(dat_genus) %>%
  full_join(dat_ceno) %>% 
  filter(coef_name != "Paleotemperature") %>% 
  ggplot(aes(y = data_source, coef_val, 
             colour = data_source)) +
  geom_vline(xintercept = 0) +
  stat_pointinterval() +
  facet_wrap(~ coef_name, 
             scales = "free") +
  scale_y_discrete(breaks = NULL) +
  scale_x_continuous(breaks = c(0)) +
  scale_color_manual(name = NULL,
                     values = c("#571415", "#C01213", 
                                "#DEA268", "#5D7A64")) +
  labs(y = NULL, 
       x = NULL) +
  theme(panel.border = element_rect(linewidth = 1, 
                                    colour = "grey80", 
                                    fill = NA), 
        legend.position = "bottom")


# save plot
ggsave(plot_beta, filename = here("figures",
                                   "beta_plot.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)


