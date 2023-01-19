# load packages
library(tidyverse)
library(here)
library(readxl)


# bin proxy data ----------------------------------------------------------

### sea level
dat_sealevel_full <- read_xlsx(here("data",
                                    "raw",
                                    "sea_level",
                                    "miller_et_al_2005.xlsx")) %>% 
  # clean up column names
  rename(sea_level = "Sea level best estimate", 
         age = "Age for best estimate (MA)")  %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(sea_level_mean = mean(sea_level)) 





### productivity
#d13C
dat_13C_full <- read_csv(here("data",
                              "raw",
                              "productivity",
                              "veizer_prokoph_2015.csv")) %>% 
  # clean up column names
  select(age = "gts2012", 
         d13C)  %>% 
  filter(age <= 180) %>% 
  drop_na(age, d13C) %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(d13C_mean = mean(d13C)) 

#Sr87/86
dat_SR_full <- read_csv(here("data",
                             "raw",
                             "productivity",
                             "mcarthur_howarth_shields_2012.csv"), 
                        col_names = FALSE) %>% 
  # clean up column names
  select(age = X1, 
         sr_value = X2) %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(sr_value_mean = mean(sr_value)) %>% 
  add_row(tibble(bin = 15, 
                 sr_value_mean = mean(read_csv(here("data",
                                                    "raw",
                                                    "productivity",
                                                    "mcarthur_howarth_shields_2012.csv"),
                                               col_names = FALSE)$X2)))



### outcrop area

# rock units
dat_marine_units_full <- read_csv(here("data",
                                       "raw",
                                       "outcrop_area",
                                       "macrostrat_24_11_2022.csv")) %>% 
  # select only marine environments
  filter(!str_detect(.$environ, "non-marine")) %>% 
  # bin to 10myr
  mutate(age = (t_age+b_age)/2, 
         bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  select(bin, unit_id, age, environ) %>% 
  count(bin) %>% 
  # add bins with zero counts
  rename(n_units = n)


# outcrop area
dat_outcrop <- read_xlsx(here("data",
                              "raw",
                              "outcrop_area",
                              "wall_ivany_wilkinson_2009.xlsx")) %>% 
  # clean up colnames
  rename(age = "age (ma)", 
         outcrop_area = "cumul_area (10^6 km^2)") %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  # use natural spline for interpolation
  { spline(.$bin, .$outcrop_area, 
           xout = 1:15) } %>% 
  as_tibble() %>% 
  rename(bin = x, outcrop_area = y) 



### shelf area
dat_cont_area_full <- read_csv(here("data",
                                    "raw",
                                    "shelf_area",
                                    "kocsis_scotese_2021.csv")) %>% 
  # clean up colnames
  rename(cont_area = "shelf-rgeos") %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(150, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(cont_area = mean(cont_area))


### temperature
dat_temp_full <- read_xlsx(here("data",
                                "raw",
                                "temperature",
                                "scotese_et_al_2021.xlsx")) %>% 
  # clean up colnames
  select(age = Age,
         temp_gat = GAT,
         temp_deep = "Deep Ocean") %>% 
  # bin to 10myr
  mutate(bin = cut(age, breaks = seq(250, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(temp_gat_binned = mean(temp_gat), 
            temp_deep_binned = mean(temp_deep))

### paleotemperature
# calculate lagged temperatures trends
dat_paleotemp <- dat_temp_full %>% 
  mutate(temp_gat_st = temp_gat_binned - lead(temp_gat_binned), 
         temp_gat_lt1 = lead(temp_gat_st), 
         temp_gat_lt2 = lead(temp_gat_st, n = 2), 
         temp_gat_lt3 = lead(temp_gat_st, n = 3), 
         temp_gat_lt4 = lead(temp_gat_st, n = 4), 
         # same for deep ocean temperature
         temp_deep_st = temp_deep_binned - lead(temp_deep_binned), 
         temp_deep_lt1 = lead(temp_deep_st), 
         temp_deep_lt2 = lead(temp_deep_st, n = 2), 
         temp_deep_lt3 = lead(temp_deep_st, n = 3), 
         temp_deep_lt4 = lead(temp_deep_st, n = 4)) %>% 
  # add missing bin
  fill(contains("temp"), .direction = "downup") %>% 
  filter(between(bin, 1, 15)) %>% 
  select(-contains("binned"))

### merge and combine full datasets
dat_proxy <- dat_13C_full %>% 
  full_join(dat_cont_area_full) %>% 
  full_join(dat_marine_units_full) %>% 
  full_join(dat_outcrop) %>% 
  full_join(dat_paleotemp) %>% 
  full_join(dat_sealevel_full) %>% 
  full_join(dat_SR_full) %>% 
  full_join(dat_temp_full %>% 
              filter(between(bin, 1, 15))) %>% 
  rename(sea_level = sea_level_mean, 
         d13C = d13C_mean, 
         sr_value = sr_value_mean)
  


# save data file
dat_proxy %>% 
  write_rds(here("data", 
                 "processed_proxy_data_10myr.rds"))




