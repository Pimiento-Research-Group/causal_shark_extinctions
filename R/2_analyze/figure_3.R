library(tidyverse)
library(here)
library(readxl)
library(deeptime)


# plotting configurations
source(here("R", "config_file.R"))

# read data ---------------------------------------------------------------

# load data on time of extinction per taxa
dat_ext <- read_delim(here("data", 
                           "fossil_occurrences", 
                           "combined_10_se_est_species_names.txt")) %>% 
  select(species, te) %>% 
  mutate(te = abs(te))


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


# preprocess data ---------------------------------------------------------


# extract time of major extinctions from literature
dat_ext_time <- tibble(age = c((95.95 - 83.35) / 2 + 83.35,
                               (73.55 - 71.75) / 2 + 71.75,
                               (66.95 - 65.75) / 2 + 65.75,
                               (57.05 - 55.65) / 2 + 55.65,
                               (38.25 - 33.55) / 2 + 33.55,
                               19,
                               (3.75 - 0) / 2 + 0), 
                       ext_event = paste0("event", 1:7)) %>% 
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


# visualise consistency ---------------------------------------------------

# create plot
plot_cons <- dat_temp %>%
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
                                            1, -7, 3, 3), 
                                      y = c(2, 2.2, 2.2, 2, -0.7, 2.1, 3)),
            size = 10/.pt, 
            colour = "grey50",
            data = dat_change) +
  annotate("text", 
           x = 15, 
           y = 17.5, 
           label = "# Extinctions", 
           colour = "grey70", 
           size = 10/.pt) +
  annotate("curve", 
           x = 25, xend = 4,
           y = 16, yend = 16, 
           curvature = 0, 
           colour = "grey80") +
  annotate("curve", 
           x = 25, xend = -2,
           y = 16, yend = 7, 
           curvature = -0.2, 
           arrow = arrow(length = unit(.2,"cm")),
           colour = "grey80") +
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



# save plot
ggsave(plot_cons, filename = here("figures",
                                   "figure_3.png"),
       width = image_width, height = image_height,
       units = image_units,
       bg = "white", device = ragg::agg_png)
