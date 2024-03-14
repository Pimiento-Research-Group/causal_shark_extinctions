library(tidyverse)
library(here)
library(patchwork)

# plotting configurations
source(here("R", "config_file.R"))

# read data ---------------------------------------------------------------

# type of threat
dat_threat <- read_rds(here("data",
                            "iucn",
                            "threat_per_spp.rds")) %>% 
  # extract category of threat
  mutate(threat = word(code, 1, sep = "\\.")) %>% 
  select(-c(url_id, name, invasive))
  



# barplot -----------------------------------------------------------------

# percentage affected 

# set up vector
vec_prop <- vector("double", 11)

for (i in 1:11) {
  vec_prop[[i]] <- dat_threat %>% 
    group_by(scientific_name) %>% 
    distinct(threat) %>% 
    filter(threat == i) %>% 
    {{ nrow(.) / length(unique(dat_threat$scientific_name))}} 
}

# make to dataframe
dat_prop <- vec_prop %>% 
  enframe(name = "threat", 
          value = "perc_af") %>%
  # https://www.iucnredlist.org/resources/threat-classification-scheme
  mutate(threat_str = case_when(
    threat == 1 ~ "Residential & commercial development", 
    threat == 2 ~ "Agriculture & aquaculture",
    threat == 3 ~ "Energy production & mining",
    threat == 4 ~ "Transportation & service corridors",
    threat == 5 ~ "Biological resource use",
    threat == 6 ~ "Human intrusions & disturbance",
    threat == 7 ~ "Natural system modifications",
    threat == 8 ~ "Invasive & other problematic species, genes & diseases",
    threat == 9 ~ "Pollution",
    threat == 10 ~ "Geological events",
    threat == 11 ~ "Climate change & severe weather",
  ))

# plot
plot_prop <- dat_prop %>% 
  mutate(threat_str = fct_reorder(threat_str, perc_af)) %>% 
  ggplot(aes(threat_str, perc_af, 
             fill = threat == 11)) +
  geom_col() +
  labs(y = "Percentage threatened by", 
       x = NULL) +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", 
                                "50", 
                                "100")) +
  scale_fill_manual(values = c("grey", "#138086ff")) +
  theme(axis.text.x = element_text(colour = "grey20", 
                                   size = 10, 
                                   angle = 40, 
                                   hjust = 1),
         legend.position = "none")



# categories --------------------------------------------------------------

# timing
plot_timing <- dat_threat %>% 
  filter(threat == 11) %>% 
  count(timing) %>% 
  mutate(fraction = n / sum(n), 
         ymax = cumsum(fraction), 
         ymin = lag(ymax, default = 0)) %>% 
  ggplot(aes(ymax = ymax, ymin = ymin,
             xmax = 4, xmin = 3,
             fill = timing)) +
  geom_rect() +
  annotate("text", 
           label = "Timing", 
           x = 0.75, 
           y = 0.8, 
           size = 13/.pt) +
  scale_fill_manual(values = c("#074d65", "#138086ff")) +
  coord_polar(theta = "y") + 
  xlim(c(0.5, 4)) +
  theme_void() +
  theme(legend.position = "none")

# severity
plot_severity <- dat_threat %>%
  filter(threat == 11) %>% 
  count(severity) %>% 
  mutate(severity = ordered(severity, 
                            levels = c("Slow, Significant Declines", 
                                       "Causing/Could cause fluctuations", 
                                       "Negligible declines", 
                                       "Unknown"))) %>% 
  drop_na() %>% 
  mutate(fraction = n / sum(n), 
         ymax = cumsum(fraction), 
         ymin = lag(ymax, default = 0)) %>% 
  ggplot(aes(ymax = ymax, ymin = ymin,
             xmax = 4, xmin = 3,
             fill = severity)) +
  geom_rect() +
  annotate("text", 
           label = "Severity", 
           x = 0.75, 
           y = 0.8, 
           size = 13/.pt) +
  scale_fill_manual(values = scales::seq_gradient_pal("#138086ff", "#23e0be", "Lab")(seq(0, 1, length.out = 4))) +
  coord_polar(theta = "y") + 
  xlim(c(0.5, 4)) +
  theme_void() +
  theme(legend.position = "none")

# scope
plot_scope <- dat_threat %>%
  filter(threat == 11) %>% 
  count(scope) %>% 
  drop_na() %>% 
  mutate(scope = ordered(scope, 
                         levels = c(
                           "Whole (>90%)",
                           "Majority (50-90%)",
                           "Minority (<50%)",
                           "Unknown"))) %>% 
  mutate(fraction = n / sum(n), 
         ymax = cumsum(fraction), 
         ymin = lag(ymax, default = 0)) %>% 
  ggplot(aes(ymax = ymax, ymin = ymin,
             xmax = 4, xmin = 3,
             fill = scope)) +
  geom_rect() +
  annotate("text", 
           label = "Scope", 
           x = 0.75, 
           y = 0.8, 
           size = 13/.pt) +
  scale_fill_manual(values = scales::seq_gradient_pal("#138086ff", "#508e74", "Lab")(seq(0, 1, length.out = 4))) +
  coord_polar(theta = "y") + 
  xlim(c(0.5, 4)) +
  theme_void() +
  theme(legend.position = "none")
  

# patch together -----------------------------------------------------------

plot_final <- plot_prop + 
  inset_element(plot_timing, 
                0.85, 0.6, 0.6, 0.85) +
  inset_element(plot_severity, 
                0.6, 0.6, 0.35, 0.85) +
  inset_element(plot_scope, 
                0.35, 0.6, 0.1, 0.85)

# save plot
ggsave(plot_final, filename = here("figures",
                                  "figure_4.svg"),
       width = image_width, height = image_height*1.5,
       units = image_units,
       bg = "white")
