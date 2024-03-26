# load packages
library(tidyverse)
library(here)
library(patchwork)
library(ggdist)
library(deeptime)
library(readxl)


# plotting configurations
source(here("R", "config_file.R"))



# load data  ------------------------------------------------------------

# get predictions from different scales/ data sets
set.seed(123)
dat_pred <- paste0("pred_temperature_", 
                   c("ceno", 
                     "genus", 
                     "stage"), 
                   ".rds") %>% 
  map_df(~ read_rds(here("data",
                         "predictions",
                         .x)) %>% 
           group_by(temperature) %>% 
           slice_sample(n = 100) %>% 
           ungroup() %>% 
           add_column(dataset = .x) %>% 
           mutate(dataset = str_remove_all(dataset, 
                                           "pred_temperature_"), 
                  dataset = str_remove_all(dataset, 
                                           ".rds")))


# predicted logit-values
dat_pred_full <- read_rds(here("data",
                               "logits",
                               "logit_full.rds"))

# load temperature data
dat_temp <- # Phanerozoic deep-ocean and global average temperature from Scotese et al 2021 
  read_xlsx(here("data",
                 "raw",
                 "temperature",
                 "scotese_et_al_2021.xlsx")) %>% 
  # clean up colnames
  select(age = Age,
         temp = "Deep Ocean") %>% 
  #Cenozoic deep-ocean temperature  used with equation 7afrom Cramer et al 2011
  bind_rows(read_xlsx(here("data",
                           "raw",
                           "temperature",
                           "cramer_et_al_2011_7a.xlsx")) %>%
              # clean up colnames
              select(age = Age,
                     temp = Temperature))

# family level logit values (temperature dependancy)
dat_family <- read_rds(here(here("data",
                                 "logits", 
                                 "logit_family_full.rds"))) %>% 
  group_by(family) %>% 
  median_qi(value) %>% 
  # add superorder
  left_join(read_rds(here("data",
                          "fossil_occurrences",
                          "database_occurrences_01_Nov_2023.rds")) %>% 
              drop_na(family, superorder) %>% 
              filter(order != "incertae sedis") %>% 
              count(superorder, family))

# load 95 bin Phanerozoic time scale based on the stratigraphic stages 
# of Gradstein et al. 2020
load(here("data", 
          "stages.Rdata"))

# time of extinction per taxa
dat_ext <- read_delim(here("data", 
                           "fossil_occurrences", 
                           "combined_10_se_est_species_names.txt")) %>% 
  select(species, te) %>% 
  mutate(te = abs(te))


# type of modern IUCN threat
dat_threat <- read_rds(here("data",
                            "iucn",
                            "threat_per_spp.rds")) %>% 
  # extract category of threat
  mutate(threat = word(code, 1, sep = "\\.")) %>% 
  select(-c(url_id, name, invasive))


# temperature-extinction relationship -------------------------------------

# visualise trends
plot_temp <- dat_pred %>%
  mutate(dataset = case_when(
    dataset == "ceno" ~ "1myr", 
    dataset == "genus" ~ "Genus", 
    dataset == "stage" ~ "Stages"
  ), 
  dataset = factor(dataset, 
                   levels = c("Stages", 
                              "Genus", 
                              "1myr"))) %>% 
  ggplot(aes(temperature,
             value,
             group = interaction(nr_draw, dataset),
             colour = dataset)) +
  geom_line(alpha = 0.2) +
  scale_colour_manual(
    values = c("#4C634C",
               colour_coral, 
               colour_purple),
    labels = c("Species - Stages", 
               "Genera - Stages", 
               "Species - Cenozoic subset"),
    name = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.9, 
                                                   linewidth = 0.6))) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  coord_cartesian(ylim = c(0, 1), 
                  xlim = c(0, 25)) +
  labs(y = "Extinction Risk [%]", 
       x = "Global Temperature [°C]") +
  theme(legend.position = c(0.75, 0.75))


plot_temp

# logit through time ------------------------------------------------------

# create plot
plot_logit <- dat_pred_full %>%
  filter(data_set %in% c("Stages", 
                         "Genus", 
                         "1myr")) %>% 
  mutate(data_set = factor(data_set, 
                           levels = c("Stages", 
                                      "Genus", 
                                      "1myr"))) %>% 
  ggplot(aes(myr, value, group = data_set)) +
  annotate("rect", 
           xmin = stages$bottom[c(81, 87)], 
           xmax = stages$top[c(81, 87)], 
           ymin = -Inf, 
           ymax = Inf, 
           fill = "#169199", 
           alpha = 0.2) + 
  annotate("rect", 
           xmin = stages$bottom[c(76, 83, 91)], 
           xmax = stages$top[c(76, 83, 91)], 
           ymin = -Inf, 
           ymax = Inf, 
           fill = "#C75E6B", 
           alpha = 0.2) + 
  annotate("text",
           y = 1.8, 
           x = 112, 
           colour = "#C75E6B", 
           size = 9/.pt, 
           label = "Hyperthermal", 
           alpha = 0.9) +
  annotate("text",
           y = 1.8, 
           x = 83, 
           colour = "#169199", 
           size = 9/.pt, 
           label = "Hypothermal", 
           alpha = 0.9) +
  geom_hline(yintercept = 0, 
             linetype = "dotted", 
             colour = colour_grey) +
  geom_linerange(aes(xmin = xmin,
                     xmax = xmax),
                 colour = "grey75") +
  geom_linerange(aes(ymin = .lower,
                     ymax = .upper),
                 colour = "grey75") +
  geom_point(aes(fill = data_set), 
             size = 2, 
             shape = 21, 
             colour = "white") +
  scale_x_reverse() +
  scale_fill_manual(values = c("#4C634C", 
                               colour_coral, 
                               colour_purple), 
                    labels = c("Species - Stages", 
                               "Genera - Stages", 
                               "Species - Cenozoic subset"), 
                    name = NULL) +
  scale_y_continuous(limits = c(-8.9, 2), 
                     breaks = seq(-8, 2, by = 2)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5))) +
  coord_geo(xlim = c(150, 0), 
            dat = list("epochs", "periods"),
            pos = list("b", "b"),
            alpha = 0.2, 
            height = unit(0.8, "line"), 
            size = list(5/.pt, 9/.pt),
            lab_color = "grey20", 
            color = "grey50", 
            abbrv = list(TRUE, FALSE), 
            fill = "white",
            expand = TRUE, 
            lwd = list(0.1, 0.2)) +
  labs(x = "Million years", 
       y = "Temperature dependancy\n[log-odds]") +
  theme(legend.position = "none") 


plot_logit




# major extinctions -------------------------------------------------------

# extract time of major extinctions from literature
dat_ext_time <- tibble(age = c((95.95 - 83.35) / 2 + 83.35,
                               (73.55 - 71.75) / 2 + 71.75,
                               (66.95 - 65.75) / 2 + 65.75,
                               (38.25 - 33.55) / 2 + 33.55,
                               19,
                               (3.75 - 0) / 2 + 0), 
                       ext_event = paste0("event", 1:6)) %>% 
  # add equally-spaced buffer to calculate number of extinctions
  mutate(age_min = age - 3, 
         age_max = age + 3)

# number of extinctions per event
dat_ext_event <- dat_ext_time %>%
  # create sequence for fuzzy matching
  mutate(age_seq = map2(age_min, 
                        age_max, 
                        ~seq(.x, .y, by = 0.0001))) %>% 
  select(ext_event, te = age_seq) %>% 
  unnest(te) %>% 
  # match with extinction
  right_join(mutate(dat_ext, 
                    te = round(te, 1))) %>% 
  drop_na(ext_event) %>% 
  count(ext_event) %>% 
  # get age of event
  left_join(dat_ext_time %>% 
              select(age_ext = age, ext_event)) %>% 
  # get closest temperature value of event
  left_join(dat_temp, join_by(closest(age_ext >= age)))


# calculate change in temperature to previous bin
dat_change <- dat_ext_event %>%
  mutate(age_seq = map(age, 
                       ~seq(.x+10, .x, by = -0.0001))) %>% 
  unnest(age_seq) %>% 
  select(ext_event, age = age_seq) %>% 
  right_join(dat_temp) %>% 
  drop_na(ext_event) %>% 
  group_by(ext_event) %>% 
  nest() %>% 
  mutate(temp_change_model = map(data,
                                 ~ lm(temp ~ age, data = .x)),
         temp_change = map_dbl(temp_change_model,
                               ~ -coef(.x)[2])) %>% 
  select(ext_event, temp_change) %>% 
  full_join(select(dat_ext_event, 
                   ext_event, n, temp, age)) %>% 
  select(ext_event, age, n, temp, temp_change) %>% 
  ungroup()



# create plot
plot_cons <-dat_temp %>%
  filter(between(age, 0, 150)) %>% 
  ggplot(aes(x = age, y = temp)) +
  geom_line(aes(colour = temp), 
            linewidth = 1) +
  geom_point(aes(size = n, 
                 shape = temp_change >= 0), 
             data = dat_change, 
             colour = "grey20", 
             fill = alpha("white", 0.4)) +
  geom_text(aes(label = n),
            position = position_nudge(x = c(3, 1.5, 4, 
                                            -7, 3, 3), 
                                      y = c(2, 2.2, 2, -0.7, 2.1, 3)),
            size = 10/.pt, 
            colour = "grey50",
            data = dat_change) +
  annotate("text", 
           x = 115, 
           y = 8, 
           label = "# Extinctions", 
           colour = "grey70", 
           size = 10/.pt) +
  scale_shape_manual(values = c(25, 24), 
                     name = NULL, 
                     labels = c("Cooling", 
                                "Warming")) +
  scale_color_gradient2(high = "coral", 
                        mid = "#97B9C1",
                        low = "#698BAB", 
                        midpoint = 13, 
                        guide = "none") +
  scale_size(range = c(2, 7), 
             guide = "none") +
  scale_x_reverse() +  
  scale_y_continuous(limits = c(-1, 23), 
                     breaks = seq(0, 20, by = 5)) +
  labs(x = "Million years", 
       y = "Global Temperature [°C]") +
  coord_geo(xlim = c(150, 0), 
            dat = list("epochs", "periods"),
            pos = list("b", "b"),
            alpha = 0.2, 
            height = unit(0.8, "line"), 
            size = list(5/.pt, 9/.pt),
            lab_color = "grey20", 
            color = "grey50", 
            abbrv = list(TRUE, FALSE), 
            fill = "white",
            expand = TRUE, 
            lwd = list(0.1, 0.2)) +
  theme(legend.position = c(0.2, 0.2), 
        legend.text = element_text(colour = "grey40"))

plot_cons


# modern threats ----------------------------------------------------------


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
    threat == 1 ~ "Habitat loss", 
    threat == 5 ~ "Exploitation",
    threat == 9 ~ "Pollution",
    threat == 11 ~ "Climate change",
    .default = "other"
  ))

# plot
plot_prop <- dat_prop %>%
  filter(threat_str != "other") %>% 
  mutate(threat_str = fct_reorder(threat_str, perc_af)) %>% 
  ggplot(aes(threat_str, perc_af, 
             fill = threat == 11)) +
  geom_col() +
  geom_text(aes(label = threat_str, 
                colour = threat_str), 
            angle = 90, 
            size = 8/.pt, 
            position = position_nudge(y = c(0.3, -0.5, 0.21, 0.35))) +
  labs(y = "Modern Threats", 
       x = NULL) +
  scale_y_continuous(breaks = c(0, 0.5, 1), 
                     labels = c("0", 
                                "50", 
                                "100")) +
  scale_fill_manual(values = c("grey", "#138086ff")) +
  scale_colour_manual(values = c("grey", "#138086ff", "grey", "white")) +
  theme(axis.text.x = element_blank(),
        legend.position = "none", 
        panel.background = element_rect(colour = "grey40", 
                                        fill = NULL))

plot_prop

# patch together and save -------------------------------------------------

plot_full <- plot_temp /
  plot_logit /
  (plot_cons +
  inset_element(plot_prop, 0.7, 0.65, 1, 1)) +
  plot_annotation(tag_levels = "a")


# save plot
ggsave(plot_full, filename = here("figures",
                                  "figure_2.svg"),
       width = image_width, height = image_height*2,
       units = image_units,
       bg = "white")


