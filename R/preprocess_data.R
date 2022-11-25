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


### Ceonozoic delta13C values from Westerhold et al 2020 ###
# load data
dat_prod_cen <- read_xlsx(here("data",
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


chronosphere::datasets() %>% View()
dat %>% 
  filter(!str_detect(.$environ, "non-marine")) %>% 
  mutate(mean_age = (t_age+b_age)/2) %>% 
  ggplot(aes(mean_age)) +
  stat_density(bounds = c(0, 150), 
               geom = "line")

str_detect(dat$environ, "non-marine")
