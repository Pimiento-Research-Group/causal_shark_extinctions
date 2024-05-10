# load packages
library(tidyverse)
library(here)
library(tidybayes)
library(patchwork)

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
  add_column(data_source = "Species - Stages")
  

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
  add_column(data_source = "Genera - Stages")


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
                                  Productivity = grep("productivity", levels(coef_name), value = TRUE), 
                                  Temperature = grep("binned", levels(coef_name), value = TRUE), 
                                  "Sea level" = grep("sea", levels(coef_name), value = TRUE), 
                                  "Shelf area" = grep("cont", levels(coef_name), value = TRUE)
         )) %>% 
  add_column(data_source = "Species - Cenozoic subset")



# visualise ---------------------------------------------------------------

# prepare data
dat_plot <- dat_main %>%
  # merge
  full_join(dat_genus) %>%
  full_join(dat_ceno) %>% 
  filter(coef_name != "Paleotemperature")

# create plot
plot_beta <- dat_plot %>%
  ggplot(aes(y = data_source, coef_val, 
             colour = data_source)) +
  geom_vline(xintercept = 0) +
  stat_pointinterval() +
  facet_wrap(~ coef_name, 
             scales = "free") +
  scale_y_discrete(breaks = NULL) +
  scale_x_continuous() +
  scale_color_manual(name = NULL,
                     values = c(colour_coral, 
                                colour_purple, 
                                "#4C634C")) +
  labs(y = NULL, 
       x = NULL) +
  theme(panel.border = element_rect(linewidth = 1, 
                                    colour = "grey80", 
                                    fill = NA), 
        legend.position = "top")



# second plot -------------------------------------------------------------

# averaged effect
plot_smr <- dat_plot %>% 
  ggplot(aes(y = coef_name, coef_val)) +
  geom_vline(xintercept = 0) +
  stat_pointinterval(shape = 21, 
                     point_fill = "white", 
                     point_size = 3) +
  geom_label(aes(label = coef_name, 
                x = median_point, 
                y = coef_name),
             size = 10/.pt,
            position = position_nudge(y = 0.4), 
            data = dat_plot %>% 
              group_by(coef_name) %>% 
              summarise(median_point = median(coef_val, 
                                           na.rm = TRUE))) +
  scale_y_discrete(breaks = NULL) +
  scale_x_continuous() +
  expand_limits(y = c(0, 5)) +
  labs(y = NULL, 
       x = NULL) 


# patch together
plot_final <- plot_beta/ plot_smr + plot_annotation(tag_levels = "A")

# save plot
ggsave(plot_final, filename = here("figures",
                                   "supplement",
                                   "beta_plot.png"), 
       width = image_width, height = image_height*1.5,
       units = image_units, 
       bg = "white", device = ragg::agg_png)


