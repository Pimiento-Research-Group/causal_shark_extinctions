# load packages
library(tidyverse)
library(here)
library(readxl)
library(deeptime)

# plotting configurations
source(here("R", "config_file.R"))


# sea level ---------------------------------------------------------------

### Jurassic to recent sea level from miller et al 2005 ###
# load data
dat_sealevel_full <- read_xlsx(here("data",
                                    "raw",
                                    "sea_level",
                                    "miller_et_al_2005.xlsx")) %>% 
  # clean up column names
  rename(sea_level = "Sea level best estimate", 
         age = "Age for best estimate (MA)")  %>% 
  # calculate 1 myr moving average
  mutate(bin = cut(age, breaks = seq(180, 0, by = -1))) %>% 
  group_by(bin) %>% 
  mutate(sea_level_mean = mean(sea_level))

# save
dat_sealevel_full %>% 
  write_rds(here("data", 
                 "sealevel_full.rds"))

# visualize
plot_sealevel_full <- dat_sealevel_full %>%
  ggplot(aes(age, sea_level)) +
  geom_hline(yintercept = 0, 
             colour = "grey70", 
             linetype = "dashed") +
  geom_point(shape = 21, 
             alpha = 0.3, 
             colour = "grey10", 
             fill = "grey50") +
  geom_line(aes(y = sea_level_mean), 
            colour = "#2c7c94", 
            linewidth = 1.3) +
  scale_x_reverse() +
  coord_geo(xlim = c(180, 0), 
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Global mean sea level [m]", 
       title = "Miller et al. 2005", 
       subtitle = "1 myr moving average")


  
### cenozoic sea level from miller et al 2020 ###
dat_sealevel_ceno <- read_xlsx(here("data",
                                    "raw",
                                    "sea_level",
                                    "miller_et_al_2020.xlsx")) %>% 
  # clean up column names
  rename(sea_level = "Sea Level (m) Smoothed", 
         age = "Age (ka)")  %>% 
  # transform to million years
  mutate(age = age/1000) 

# save
dat_sealevel_ceno %>% 
  write_rds(here("data", 
                 "sealevel_ceno.rds"))

# visualize
plot_sealevel_ceno <- dat_sealevel_ceno %>%
  ggplot(aes(age, sea_level)) +
  geom_hline(yintercept = 0, 
             colour = "grey70", 
             linetype = "dashed") +
  geom_line(colour = "#2c7c94", 
            linewidth = 1.3) +
  scale_x_reverse() +
  coord_geo(xlim = c(65, 0), 
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Global mean sea level [m]", 
       title = "Miller et al. 2020", 
       subtitle = "Smoothed to 20ka") +
  theme(plot.title = element_text())




# productivity ------------------------------------------------------------


### Cenozoic delta13C values from Westerhold et al 2020 ###
# load data
dat_13C_ceno <- read_xlsx(here("data",
                                "raw",
                                "productivity",
                                "westerhold_et_al_2020.xlsx"),
                           sheet = "Table S34",
                           skip = 1) %>% 
  # remove redundant columns
  select(age = age_tuned, fine_loess = ISOBENd13cLOESSsmooth, 
         coarse_loess = ISOBENd13cLOESSsmoothLongTerm)

# save
dat_13C_ceno %>% 
  write_rds(here("data", 
                 "d13C_ceno.rds"))

# visualize
plot_13C_ceno <- dat_13C_ceno %>%
  ggplot(aes(age, fine_loess)) +
  geom_hline(yintercept = 0, 
             colour = "grey70", 
             linetype = "dashed") +
  geom_line(colour = "#c4aa23", 
            linewidth = 1.3) +
  geom_line(aes(y = coarse_loess), 
            colour = "#fbe45b", 
            linewidth = 1.3) +
  scale_x_reverse() +
  coord_geo(xlim = c(70, 0), 
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Benthic d13C [‰]", 
       title = "Westerhold et al. 2020", 
       subtitle = "Smoothed to 20ka (dark yellow) and 1myr (light yellow)") +
  theme(plot.title = element_text())


### Phanerozoic delta13C values from Veizer and Prokoph ###
# load data
dat_13C_full <- read_csv(here("data",
                               "raw",
                               "productivity",
                               "veizer_prokoph_2015.csv")) %>% 
  # clean up column names
  select(age = "gts2012", 
         d13C)  %>% 
  filter(age <= 180) %>% 
  drop_na(age, d13C) %>% 
  # calculate 5 myr moving average
  mutate(bin = cut(age, breaks = seq(180, -5, by = -5))) %>% 
  group_by(bin) %>% 
  mutate(d13C_mean = mean(d13C))

# save
dat_13C_full %>% 
  write_rds(here("data", 
                 "13C_full.rds"))

# visualize
plot_13C_full <- dat_13C_full %>%
  ggplot(aes(age, d13C)) +
  geom_hline(yintercept = 0, 
             colour = "grey70", 
             linetype = "dashed") +
  geom_point(shape = 21, 
             alpha = 0.3, 
             colour = "grey10", 
             fill = "grey50") +
  geom_line(aes(y = d13C_mean), 
            colour = "#fbe45b", 
            linewidth = 1.3) +
  scale_x_reverse() +
  coord_geo(xlim = c(180, 0), 
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Benthic d13C [‰]", 
       title = "Veizer and Prokoph 2015", 
       subtitle = "5 myr moving average") +
  theme(plot.title = element_text())



### Phanerozoic 87Sr/86Sr values from McArthur, Howarth, Shields 2012 ###
# load data
dat_SR_full <- read_csv(here("data",
                              "raw",
                              "productivity",
                              "mcarthur_howarth_shields_2012.csv"), 
                        col_names = FALSE) %>% 
  # clean up column names
  select(age = X1, 
         sr_value = X2) 
# save
dat_SR_full %>% 
  write_rds(here("data", 
                 "SR_full.rds"))

# visualize
plot_SR_full <- dat_SR_full %>%
  ggplot(aes(age, sr_value)) +
  geom_hline(yintercept = 0, 
             colour = "grey70", 
             linetype = "dashed") +
  geom_line(colour = "#fbe45b", 
            linewidth = 1.3) +
  scale_x_reverse() +
  coord_geo(xlim = c(140, 0), 
            ylim = c(0.7060, 0.7095),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "87Sr / 86Sr", 
       title = "Ccarthur, Howarth, Shields 2012", 
       subtitle = "Lowess fit curve") +
  theme(plot.title = element_text())



# outcrop area ------------------------------------------------------------


### Phanerozoic rock units from the macrostrat database ###
# load data
dat_outcrop_full <- read_csv(here("data",
                                  "raw",
                                  "outcrop_area",
                                  "macrostrat_24_11_2022.csv")) %>% 
  # select only marine environments
  filter(!str_detect(.$environ, "non-marine")) %>% 
  mutate(mean_age = (t_age+b_age)/2) %>% 
  select(unit_id, t_age, b_age, mean_age, environ) 

# save
dat_outcrop_full %>% 
  write_rds(here("data", 
                 "outcrop_full.rds"))

# visualize
plot_outcrop_full <- dat_outcrop_full %>%
  # bin to 5 myr
  mutate(bin = cut(mean_age, breaks = seq(70, 0, by = -5))) %>% 
  count(bin) %>% 
  drop_na(bin) %>% 
  mutate(age = seq(2.5, 67.5, by = 5)) %>% 
  ggplot(aes(age, n)) +
  geom_line(linewidth = 1.5, colour = "#a6d0c8") +
  scale_x_reverse() +
  coord_geo(xlim = c(70, 0), 
            ylim = c(0, 3200),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Marine rock units", 
       title = "Macrostrat", 
       subtitle = "Binned into 5 myr") +
  theme(plot.title = element_text())


### Phanerozoic rock units from the macrostrat database ###
# load data
dat_marine_units <- read_csv(here("data",
                                  "raw",
                                  "outcrop_area",
                                  "macrostrat_24_11_2022.csv")) %>% 
  # select only marine environments
  filter(!str_detect(.$environ, "non-marine")) %>% 
  mutate(mean_age = (t_age+b_age)/2) %>% 
  select(unit_id, t_age, b_age, mean_age, environ) 

# save
dat_marine_units %>% 
  write_rds(here("data", 
                 "marine_units_full.rds"))

# visualize
plot_marine_units <- dat_marine_units %>%
  # bin to 5 myr
  mutate(bin = cut(mean_age, breaks = seq(70, 0, by = -5))) %>% 
  count(bin) %>% 
  drop_na(bin) %>% 
  mutate(age = seq(2.5, 67.5, by = 5)) %>% 
  ggplot(aes(age, n)) +
  geom_line(linewidth = 1.5, colour = "#a6d0c8") +
  scale_x_reverse() +
  coord_geo(xlim = c(70, 0), 
            ylim = c(0, 3200),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Marine rock units", 
       title = "Macrostrat", 
       subtitle = "Binned into 5 myr") +
  theme(plot.title = element_text())


### Phanerozoic outcrop area from Wall, Ivany, Wilkinson 2009 ###
# load data
dat_outcrop_full <- read_xlsx(here("data",
                                  "raw",
                                  "outcrop_area",
                                  "wall_ivany_wilkinson_2009.xlsx")) %>% 
  # clean up colnames
  rename(age = "age (ma)", 
         area = "cumul_area (10^6 km^2)")

# save
dat_outcrop_full %>% 
  write_rds(here("data", 
                 "outcrop_full.rds"))

# visualize
plot_outcrop_full <- dat_outcrop_full %>%
  ggplot(aes(age, area)) +
  geom_line(linewidth = 1.5, colour = "#a6d0c8") +
  scale_x_reverse() +
  coord_geo(xlim = c(160, 0), 
            ylim = c(0, 6),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Marine outcrop area [10^6 km^2", 
       title = "Wall, Ivany, Wilkinson 2009") +
  theme(plot.title = element_text())



# shelf area --------------------------------------------------------------



# temperature -------------------------------------------------------------


