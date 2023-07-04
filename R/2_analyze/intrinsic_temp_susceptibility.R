library(tidyverse)
library(here)
library(brms)
library(tidybayes)


# plotting configurations
source(here("R", "config_file.R"))



# read data ---------------------------------------------------------------

# fossil data on species level
dat_fossil <- read_rds(here("data",
                            "processed_merged_data.rds")) 

# fossil data on genus level
dat_fossil_genus <-  read_rds(here("data",
                                   "processed_fossil_data_genus.rds"))

# fossil data on species level over Cenozoic
dat_fossil_ceno <- read_rds(here("data",
                                 "processed_fossil_data_genus.rds")) 

# modern red list data
dat_modern <- read_delim(here("data",
                              "iucn",
                              "species_data.txt"), 
                         col_names = FALSE) %>% 
  rename(species = X1, 
         status = X2)



# fit models --------------------------------------------------------------

# start with species level
dat_merged <- dat_fossil 

# fit model
mod_fossil <- brm_logistic("ext_signal ~ temp_gat_binned +  (temp_gat_binned | genus)")

# extract ranks
risk_fossil <- tibble(temp_gat_binned = 10,
       genus = distinct(dat_merged, genus) %>% 
         pull()) %>% 
  add_linpred_draws(mod_fossil) %>% 
  group_by(genus) %>% 
  summarise(median_risk = median(.linpred))

# same for genus level
dat_merged <- dat_fossil_genus

mod_genus <- brm_logistic("ext_signal ~ temp_gat_binned + (temp_gat_binned | genus)")

risk_genus <- tibble(temp_gat_binned = 10,
                     genus = distinct(dat_merged, genus) %>% 
                       pull()) %>% 
  add_linpred_draws(mod_genus) %>% 
  group_by(genus) %>% 
  summarise(median_risk = median(.linpred))

# and for cenozoic species
dat_merged <- dat_fossil_ceno 

mod_ceno <- brm_logistic("ext_signal ~ temp_binned + (temp_binned | genus)")

risk_ceno <- tibble(temp_binned = 10,
                     genus = distinct(dat_merged, genus) %>% 
                       pull()) %>% 
  add_linpred_draws(mod_ceno) %>% 
  group_by(genus) %>% 
  summarise(median_risk = median(.linpred))



# merge data --------------------------------------------------------------

dat_combined <- dat_modern %>% 
  mutate(genus = word(species, 1), 
         status = factor(status,
                         levels = c("DD",
                                    "NE",
                                    "LC", 
                                    "NT", 
                                    "VU", 
                                    "EN", 
                                    "CR"))) %>% 
  filter(!status %in% c("DD", "NE")) %>% 
  left_join(risk_fossil %>%
              filter(genus %in% word(dat_modern$species, 1)) %>% 
              rename(risk_fos = median_risk) %>% 
              mutate(risk_fos = rank(risk_fos))) %>%
  left_join(risk_ceno %>% 
              filter(genus %in% word(dat_modern$species, 1)) %>% 
              rename(risk_ceno = median_risk) %>% 
              mutate(risk_ceno = rank(risk_ceno))) %>%
  left_join(risk_genus %>% 
              filter(genus %in% word(dat_modern$species, 1)) %>% 
              rename(risk_genus = median_risk) %>% 
              mutate(risk_genus = rank(risk_genus))) %>%
  pivot_longer(cols = c(risk_fos, risk_ceno, risk_genus), 
               names_to = "scale", 
               values_to = "median_risk") %>% 
  drop_na(median_risk) 

# save data
dat_combined %>% 
  write_rds(here("data", 
                 "iucn", 
                 "temp_susceptibility.rds"))



# visualise ---------------------------------------------------------------

plot_dep <- dat_combined %>%
  ggplot(aes(status, median_risk)) +
  geom_point(aes(fill = scale),
             position = position_jitter(width = 0.1,
                                        seed = 123),
              alpha = 0.6, 
              shape = 21, 
              colour = "white",
              size = 2.3) +
  stat_halfeye(alpha = 0.5, 
               position = position_nudge(x = 0.1),
               limits = c(1, NA),
               fill = "white",
               point_interval = "mean_qi",
               shape = 21,
               point_alpha = 1,
               point_size = 4,
               point_fill = "white",
               interval_alpha = 0.3, 
               .width = 0.55) +
  geom_smooth(aes(as.numeric(status)),
              position = position_nudge(x = -1.9),
              method = "lm",
              se = FALSE,
              colour = "white",
              linewidth = 3) +
  geom_smooth(aes(as.numeric(status)),
              position = position_nudge(x = -1.9),
              method = "lm",
              se = FALSE,
              colour = colour_yellow,
              linewidth = 1) +
  annotate("text",
           y = 7, 
           x = 5.8, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "High susceptibility") +
  annotate("text",
           y = 58, 
           x = 5.8, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "Low susceptibility") +
  annotate("curve",
           y = 17, 
           yend = 47, 
           x = 5.8, 
           xend = 5.8,
           colour = "grey50", 
           curvature = 0,
           arrow = arrow(length = unit(.2,"cm"), 
                         ends = "both")) +
  annotate("text",
           y = 72, 
           x = 1.2, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "Low risk") +
  annotate("text",
           y = 72, 
           x = 5, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "High risk") +
  annotate("curve",
           y = 72, 
           yend = 72, 
           x = 1.5, 
           xend = 4.7,
           colour = "grey50", 
           curvature = 0,
           arrow = arrow(length = unit(.2,"cm"), 
                         ends = "both")) +
  scale_fill_manual(values = c("#4C634C", 
                               colour_coral, 
                               colour_purple), 
                    labels = c("Stages - Species", 
                               "Stages - Genus", 
                               "1myr - Species"), 
                    name = NULL) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.9))) +
  scale_y_continuous(breaks = c(1, 20, 40, 60)) +
  labs(x = "IUCN red list status", 
       y = "Fossil temperature dependancy [ranked]") +
  theme(legend.position = "bottom") +
  coord_flip()


# save plot
ggsave(plot_dep, filename = here("figures",
                                 "fossil_temp_dependancy.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)     



