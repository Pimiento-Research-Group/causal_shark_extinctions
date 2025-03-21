library(tidyverse)
library(here)
library(brms)
library(tidybayes)
library(patchwork)


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

# temperature data
dat_ipcc <- tibble(filename = list.files(here("data",
                                              "ipcc")) %>%
                     str_remove("tas_global_") %>%
                     str_remove(".csv")) %>% 
  mutate(file_contents = map(list.files(here("data",
                                             "ipcc"),
                                        full.names = TRUE),  
                             ~ read_csv(.x))) %>% 
  filter(filename %in% c("Historical", "SSP1_1_9")) %>% 
  unnest(file_contents) %>% 
  select(year = Year, temp = Mean) %>% 
  filter(year %in% 1990:2020)


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

mod_ceno <- brm_logistic("ext_signal ~ temp_gat_binned + (temp_gat_binned | genus)")

risk_ceno <- tibble(temp_gat_binned = 10,
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
                 "predictions", 
                 "ranked_iucn_risk.rds"))


# fit dependency model ----------------------------------------------------

mod_dep <- brm(
  bf(ext_signal ~ median_risk),
  family = "bernoulli",
  data = dat_combined %>% 
    mutate(ext_signal = if_else(status == "LC",
                                0, 1)),
  seed = 1708,
  iter = 10000,
  warmup = 2000,
  chains = 4,
  cores = 4,
  threads = threading(4),
  backend = "cmdstanr"
  )

# visualise ---------------------------------------------------------------

plot_dep <- tibble(median_risk = 1:70) %>%
  add_epred_draws(mod_dep, 
                  ndraws = 50, 
                  seed = 1) %>% 
  ggplot(aes(median_risk)) +
  stat_dots(data = dat_combined %>% 
              mutate(ext_signal = if_else(status == "LC",
                                          0, 1), 
                     ext_signal_fct = factor(ext_signal)),
            aes(y = ext_signal, 
                side = ifelse(ext_signal == 0, "top", "bottom"), 
                colour = ext_signal_fct),
            scale = 0.89,
            shape = 3, 
            alpha = 0.2) +
  geom_line(aes(y = .epred, 
                group = .draw), 
            alpha = 0.25,
            colour = colour_yellow) +
  scale_color_manual(values = c(colour_grey, "#D55E00")) +
  annotate("text",
           y = -0.07, 
           x = 10, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "Low susceptibility") +
  annotate("text",
           y = -0.07, 
           x = 55, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "High susceptibility") +
  annotate("curve",
           y = -0.07, 
           yend = -0.07, 
           x = 20, 
           xend = 45,
           colour = "grey50", 
           curvature = 0,
           arrow = arrow(length = unit(.2,"cm"), 
                         ends = "both")) +
  scale_y_continuous(labels = function(.x) .x*100, 
                     limits = c(-0.1, 1)) +
  scale_x_continuous(breaks = c(1, 20, 40, 60)) +
  labs(y = "Percentage threatened", 
       x = "Fossil temperature dependancy [ranked]") +
  theme(legend.position = "none")





# increase in risk --------------------------------------------------------

# prepare data
dat_merged <- dat_modern %>% 
  filter(!status %in% c("DD", "NE")) %>% 
  mutate(ext_signal = if_else(status %in% c("LC", "NT"),
                              0, 1))

# fit model to estimate average extinction risk in modern ocean
mod1 <- brm_logistic("ext_signal ~ 1")

# read in change in modern risk due to temperature
dat_pred_modern <- read_rds(here("data",
                                 "logits",
                                 "logit_modern.rds"))
            
# calculate the buffer factor of temperature adaptation
buffer_factor <- dat_ipcc$temp[dat_ipcc$year == 2019] - 
  dat_ipcc$temp[dat_ipcc$year == 1990] * 
  (1-exp(dat_pred_modern$value))+1

# estimate change in risk
dat_risk <- posterior_epred(mod1, ndraws = 1000)[,1] %>% 
  as_tibble() %>% 
  pivot_longer(everything(),
               names_to = "draw", 
               values_to = "ext_risk") %>% 
  mutate(without_temp = ext_risk * buffer_factor) %>% 
  # # calculate contrast
  # mutate(risk_contr = without_temp - ext_risk) %>%
  # median_qi(risk_contr)
  pivot_longer(cols = -draw,
               names_to = "ext_est", 
               values_to = "ext_risk") 


# save data
dat_risk %>% 
  write_rds(here("data", 
                 "predictions", 
                 "change_modern_risk.rds"))

# visualise
plot_risk <- dat_risk %>%
  ggplot(aes(ext_est, ext_risk)) +
  geom_point(position = position_jitter(width = 0.1, 
                                        seed = 123),
             shape = 21, 
             fill = "grey20", 
             colour = "white", 
             alpha = 0.1) +
  stat_pointinterval(position = position_nudge(x = 0.15), 
               limits = c(0, 1),
               shape = 21,
               point_alpha = 1,
               point_fill = "white") +
  annotate("text",
           y = 0.2, 
           x = 1, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "Empirical") +
  annotate("text",
           y = 0.3, 
           x = 2, 
           colour = "grey40", 
           size = 10/.pt, 
           label = "Without adaptation\n(simulated)") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     labels = c("0", "25", "50", "75", "100"), 
                     limits = c(0, 1), 
                     name = "Percentage threatened") +
  scale_x_discrete(labels = c("Empirical", 
                              "Without adaptation"), 
                   name = NULL) +
  annotate("text",
           y = 0.65, 
           x = 1.55, 
           size = 10/.pt, 
           label = "+28pp [22, 34]", 
           colour = colour_yellow, 
           fontface = "bold", 
           angle = 20) +
  theme(axis.ticks.x = element_blank(), 
        axis.text.x = element_blank())
  
  



# save plots --------------------------------------------------------------

plot_full <- plot_risk /
  plot_dep +
  plot_layout(heights = c(1.5, 1)) +
  plot_annotation(tag_levels = "a")

# save plot
ggsave(plot_full, filename = here("figures",
                                 "modern_ext_risk.png"), 
       width = image_width, height = image_height*1.5,
       units = image_units, 
       bg = "white", device = ragg::agg_png)     



