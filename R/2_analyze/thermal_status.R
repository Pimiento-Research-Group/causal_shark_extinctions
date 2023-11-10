# load packages
library(tidyverse)
library(here)
library(brms)
library(tidybayes)

# plotting configurations
source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

data(stages, package = "divDyn")

# fossil data on species level
dat_fossil <- read_rds(here("data",
                            "processed_merged_data.rds")) %>% 
  # assign thermal status
  mutate(therm_stat = case_when(bin %in% c(81, 87) ~ "Hypothermal",
                                bin %in% c(76, 83, 91) ~ "Hyperthermal",
                                .default = "Background")) %>% 
  filter(therm_stat != "Background")

# fossil data on genus level
dat_genus <- read_rds(here("data",
                           "processed_fossil_data_genus.rds")) %>% 
  # assign thermal status
  mutate(therm_stat= case_when(bin %in% c(81, 87) ~ "Hypothermal",
                               bin %in% c(76, 83, 91) ~ "Hyperthermal",
                               .default = "Background")) %>% 
  filter(therm_stat != "Background")


# fossil data covering cenozoic at higher resolution
dat_ceno <- read_rds(here("data",
                          "processed_fossil_data_cenozoic.rds")) %>% 
  # assign thermal status
  mutate(therm_stat = case_when(bin %in% c(66:72, 34:38) ~ "Hypothermal",
                                 bin %in% c(94:101, 56:62, 12:16) ~ "Hyperthermal",
                                 .default = "Background")) %>% 
  filter(therm_stat != "Background")


# model risk --------------------------------------------------------------

# first stages and species
dat_merged <- dat_fossil

mod_fossil <- brm_logistic("ext_signal ~ therm_stat")

risk_fossil <- tibble(therm_stat = c("Hypothermal",
                                     "Hyperthermal")) %>%
  add_epred_draws(mod_fossil, 
                  ndraws = 1000) %>% 
  group_by(therm_stat) %>% 
  median_qi(.epred) %>% 
  add_column(scale = "Stages - Species")

# stages and genera
dat_merged <- dat_genus

mod_genus <- brm_logistic("ext_signal ~ therm_stat")

risk_genera <- tibble(therm_stat = c("Hypothermal",
                                     "Hyperthermal")) %>% 
  add_epred_draws(mod_genus, 
                  ndraws = 1000) %>% 
  group_by(therm_stat) %>% 
  median_qi(.epred) %>% 
  add_column(scale = "Stages - Genus")


# and 1 myr bins and species over Cenozoic
dat_merged <- dat_ceno

mod_ceno <- brm_logistic("ext_signal ~ therm_stat")

risk_ceno <- tibble(therm_stat = c("Hypothermal",
                                   "Hyperthermal")) %>% 
  add_epred_draws(mod_ceno, 
                  ndraws = 1000) %>% 
  group_by(therm_stat) %>% 
  median_qi(.epred) %>% 
  add_column(scale = "1myr - Species")




# visualise ---------------------------------------------------------------

dat_thermal <- risk_fossil %>% 
  bind_rows(risk_genera) %>% 
  bind_rows(risk_ceno) 
  # %>% select(therm_stat, .epred, scale) %>%
  # pivot_wider(names_from = therm_stat,
  #             values_from = .epred) %>%
  # mutate(perc_inc = (Hypothermal - Hyperthermal)/ Hypothermal) %>% 
  # median_qi(perc_inc)
  

# save data
dat_thermal %>% 
  write_rds(here("data", 
                 "predictions", 
                 "thermal_plot_data.rds"))
  
# plot
plot_thermal <- dat_thermal %>% 
  ggplot(aes(therm_stat)) +
  geom_linerange(aes(ymin = .lower, 
                     ymax = .upper, 
                     colour = scale), 
                 position = position_dodge(width = 0.4)) +
  geom_point(aes(y = .epred, 
                 fill = scale), 
             shape = 21, 
             colour = "white", 
             size = 2, 
             stroke = 0.5, 
             position = position_dodge(width = 0.4)) +
  annotate("text",
           y = 0.26, 
           x = 1.5, 
           size = 8/.pt, 
           label = "+66%\n[+56%, +79%]", 
           colour = "grey40", 
           fontface = "bold") +
  annotate("rect", 
           xmin = 0.8, 
           xmax = 1.2, 
           ymin = -Inf, 
           ymax = Inf, 
           fill = "#C75E6B", 
           alpha = 0.2) + 
  annotate("rect", 
           xmin = 1.8, 
           xmax = 2.2,
           ymin = -Inf, 
           ymax = Inf, 
           fill = "#169199", 
           alpha = 0.2) + 
  labs(y = "Extinction risk [%]", 
       x = NULL) +
  scale_fill_manual(values = c(colour_purple, 
                               colour_coral, 
                               "#4C634C")) +
  scale_colour_manual(values = c(colour_purple, 
                                 colour_coral, 
                                 "#4C634C")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5), 
                     labels = c(0, 25, 50), 
                     limits = c(0, 0.55), 
                     name = "Extinction risk [%]") +
  theme(legend.position = "none", 
        axis.ticks.x = element_blank(), 
        axis.text.x = element_blank())

# save as rds
plot_thermal %>% 
  write_rds(here("data", 
                 "predictions", 
                 "thermal_plot.rds"))
