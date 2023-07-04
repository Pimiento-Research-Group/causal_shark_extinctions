dat_fossil -> dat_merged

mod_fossil <- brm_logistic("ext_signal ~ temp_gat_binned +  (temp_gat_binned | genus)")

risk_fossil <- tibble(temp_gat_binned = 10,
       genus = distinct(dat_merged, genus) %>% 
         pull()) %>% 
  add_linpred_draws(mod_fossil) %>% 
  group_by(genus) %>% 
  summarise(median_risk = median(.linpred))

dat_fossil_ceno -> dat_merged
mod_ceno <- brm_logistic("ext_signal ~ temp_binned + (temp_binned | genus)")

risk_ceno <- tibble(temp_binned = 10,
                     genus = distinct(dat_merged, genus) %>% 
                       pull()) %>% 
  add_linpred_draws(mod_ceno) %>% 
  group_by(genus) %>% 
  summarise(median_risk = median(.linpred))

dat_merged <- read_rds(here("data",
                            "processed_fossil_data_genus.rds")) %>% 
  mutate(bin = as.factor(bin))

mod_genus <- brm_logistic("ext_signal ~ temp_gat_binned + (temp_gat_binned | genus)")

risk_genus <- tibble(temp_gat_binned = 10,
                    genus = distinct(dat_merged, genus) %>% 
                      pull()) %>% 
  add_linpred_draws(mod_genus) %>% 
  group_by(genus) %>% 
  summarise(median_risk = median(.linpred))

dat_modern <- read_delim(here("data",
                              "iucn",
                              "species_data.txt"), 
                         col_names = FALSE) %>% 
  rename(species = X1, 
         status = X2)

dat_modern %>% 
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
  drop_na(median_risk) %>% 
  ggplot(aes(status, median_risk)) +
  geom_point(aes(fill = scale),
             position = position_jitter(width = 0.05,
                                        seed = 123),
              alpha = 0.7, 
              shape = 21, 
              colour = "white",
              size = 3) +
  stat_slab(alpha = 0.5, 
               position = position_nudge(x = 0.1), 
            fill = "grey93", 
            slab_colour = "grey60", 
            slab_linewidth = 0.6) +
  geom_smooth(aes(as.numeric(status)),
              position = position_nudge(x = -1.9),
              method = "lm", 
              se = FALSE, 
              colour = "white", 
              linewidth = 4.5) +
  geom_smooth(aes(as.numeric(status)),
              position = position_nudge(x = -1.9),
              method = "lm", 
              se = FALSE, 
              colour = colour_yellow, 
              linewidth = 1.5) +
  scale_fill_manual(values = c("#4C634C",
                               colour_coral,
                               colour_purple)) +
  theme(legend.position = "none")


  