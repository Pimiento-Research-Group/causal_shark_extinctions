# load packages
library(tidyverse)
library(here)
library(readxl)


# plotting configurations
source(here("R", "config_file.R"))



# stage data --------------------------------------------------------------

# load 95 bin Phanerozoic time scale based on the stratigraphic stages 
# of Gradstein et al. 2020
load(here("data", 
          "stages.Rdata"))

# clean up for binning
dat_stages <- stages %>% 
  as_tibble() %>% 
  select(stg, bottom, top)


# sea level ---------------------------------------------------------------

### Jurassic to recent sea level from miller et al 2005 ###
# load data
dat_sealevel <- read_xlsx(here("data",
                               "raw",
                               "sea_level",
                               "miller_et_al_2005.xlsx")) %>% 
  # clean up column names
  rename(sea_level = "Sea level best estimate", 
         age = "Age for best estimate (MA)")  %>% 
  # bin to stages
  mutate(bin = 95 - cut(age, breaks = dat_stages$bottom,
                        include.lowest = TRUE,
                        labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(sea_level = mean(sea_level)) %>% 
  ungroup() %>% 
  select(bin, sea_level) %>% 
  filter(between(bin, 69, 94)) 



# productivity ------------------------------------------------------------


### Phanerozoic delta13C values from Veizer and Prokoph ###
# load data
dat_13C <- read_csv(here("data",
                         "raw",
                         "productivity",
                         "veizer_prokoph_2015.csv")) %>% 
  # clean up column names
  select(age = "gts2012", 
         d13C)  %>% 
  filter(age <= 180) %>% 
  drop_na(age, d13C) %>% 
  # bin to stages
  mutate(bin = 95 - cut(age, breaks = dat_stages$bottom,
                        include.lowest = TRUE,
                        labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(d13C = mean(d13C)) %>% 
  ungroup() %>% 
  select(bin, d13C) %>% 
  filter(between(bin, 69, 94)) 



### Phanerozoic 87Sr/86Sr values from McArthur, Howarth, Shields 2012 ###
# load data
dat_SR <- read_csv(here("data",
                        "raw",
                        "productivity",
                        "mcarthur_howarth_shields_2012.csv"), 
                        col_names = FALSE) %>% 
  # clean up column names
  select(age = X1, 
         sr_value = X2) %>% 
  # bin to stages
  mutate(bin = 95 - cut(age, breaks = dat_stages$bottom,
                        include.lowest = TRUE,
                        labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(sr_mean = mean(sr_value)) %>% 
  ungroup() %>% 
  select(bin, sr_mean) %>% 
  filter(bin >= 69) %>%
  complete(bin = 69:94) %>%
  fill(sr_mean, .direction = "downup")



# outcrop area ------------------------------------------------------------


### Phanerozoic rock units from the macrostrat database ###
# load data
dat_marine_units <- read_csv(here("data",
                                  "raw",
                                  "outcrop_area",
                                  "macrostrat_24_11_2022.csv")) %>% 
  # select only marine environments
  filter(!str_detect(.$environ, "non-marine")) %>% 
  # bin to stages
  mutate(age = (t_age+b_age)/2, 
         bin = 95 - cut(age, breaks = dat_stages$bottom,
                        include.lowest = TRUE,
                        labels = FALSE)) %>% 
  drop_na(bin) %>% 
  select(bin, unit_id, age, environ) %>% 
  count(bin, 
        name = "n_units") %>% 
  arrange(bin) %>% 
  filter(bin >= 69)

### Phanerozoic outcrop area from Wall, Ivany, Wilkinson 2009 ###
# load data
dat_outcrop <- read_xlsx(here("data",
                              "raw",
                              "outcrop_area",
                              "wall_ivany_wilkinson_2009.xlsx")) %>% 
  # clean up colnames
  rename(age = "age (ma)", 
         outcrop_area = "cumul_area (10^6 km^2)") %>% 
  # bin to stages
  mutate(bin = 95 - cut(age, breaks = dat_stages$bottom,
                        include.lowest = TRUE,
                        labels = FALSE)) %>% 
  drop_na(bin) %>% 
  select(bin, outcrop_area, age) %>% 
  # use natural spline for interpolation
  { spline(.$bin, .$outcrop_area, 
           xout = 69:95) } %>% 
  as_tibble() %>% 
  rename(bin = x, outcrop_area = y) %>% 
  filter(between(bin, 69, 94))


# shelf area --------------------------------------------------------------

### Phanerozoic flooded continental area as proportion of earths surface area from Kocsis and Scotese 2021  ###
# load data
dat_cont_area <- read_csv(here("data",
                                    "raw",
                                    "shelf_area",
                                    "kocsis_scotese_2021.csv")) %>% 
  # clean up colnames
  rename(cont_area = "shelf-rgeos") %>% 
  # bin to stages
  mutate(bin = 95 - cut(age, breaks = dat_stages$bottom,
                        include.lowest = TRUE,
                        labels = FALSE)) %>% 
  drop_na(bin) %>% 
  select(bin, cont_area) %>% 
  filter(bin >= 69) %>%
  complete(bin = 69:94) %>%
  fill(cont_area, .direction = "downup")



# temperature -------------------------------------------------------------


### Phanerozoic deep-ocean and global average temperature from Scotese et al 2021 ###
# load data
dat_temp <- read_xlsx(here("data",
                                "raw",
                                "temperature",
                                "scotese_et_al_2021.xlsx")) %>% 
  # clean up colnames
  select(age = Age,
         temp_gat = GAT,
         temp_deep = "Deep Ocean") %>% 
  # bin to stages
  mutate(bin = 95 - cut(age, breaks = dat_stages$bottom,
                        include.lowest = TRUE,
                        labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(temp_gat_binned = mean(temp_gat), 
         temp_deep_binned = mean(temp_deep)) %>% 
  ungroup() %>% 
  select(bin, temp_gat_binned, temp_deep_binned) %>% 
  filter(between(bin, 69, 94)) %>% 
  complete(bin = 69:94) %>%
  fill(c(temp_gat_binned, temp_deep_binned), .direction = "downup")




# paleotemperature --------------------------------------------------------

# calculate lagged temperatures trends
dat_paleotemp <- dat_temp %>% 
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
  filter(bin >= 69) %>%
  complete(bin = 69:95) %>% 
  fill(contains("temp"), .direction = "downup") %>% 
  select(-contains("binned")) %>% 
  filter(between(bin, 69, 94)) 


# extinction signal -------------------------------------------------------

# get those species that survive until now
spec_modern <- read_delim(here("data",
                                   "fossil_occurrences",
                                   "combined_10_se_est_species_names.txt")) %>% 
  filter(te == 0) %>% 
  pull(species)

# load PyRate estimates
dat_ext <- read_delim(here("data",
                           "fossil_occurrences",
                           "combined_10_se_est_species_names.txt")) %>% 
  # estimate fad and lad for each species
  mutate(across(c(ts, te), abs)) %>% 
  group_by(species) %>% 
  summarise(ori_age = mean(ts),
            ext_age = mean(te)) %>% 
  # bin fad and lad to stages
  mutate(bin_ori = 95 - cut(ori_age, breaks = dat_stages$bottom,
                            include.lowest = TRUE,
                            labels = FALSE), 
         bin_ext = 96 - cut(ext_age, breaks = dat_stages$top,
                            include.lowest = TRUE,
                            labels = FALSE)) %>% 
  drop_na(bin_ori, bin_ext) %>% 
  # fill in duration bins
  mutate(bin_occ = map2(.x = bin_ori, 
                        .y = bin_ext,
                        .f = ~ seq(.x, .y, by = 1))) %>% 
  select(species, bin_occ, bin_ext) %>% 
  unnest(bin_occ) %>% 
  # create extinction signal
  group_by(species) %>% 
  mutate(ext_signal = if_else(bin_occ == bin_ext, 1, 0), 
         # assign 0 to those species that survive until the modern
         ext_signal = if_else(species %in% spec_modern, 
                              0, ext_signal)) %>% 
  # clean up
  ungroup() %>% 
  select(-bin_ext, 
         bin = bin_occ) %>% 
  mutate(accepted_name = str_replace(species, "_", " ")) %>%
  select(-species) 


# sampling effort ---------------------------------------------------------

# load occurrence database
dat_occurrences <- read_rds(here("data",
                                 "fossil_occurrences",
                                 "database_occurrences_13_Sept_2024.rds"))

# bin the occurrences to stages
dat_occ_binned <- dat_occurrences %>% 
  mutate(bin_min = 95 - cut(Min_Ma, breaks = dat_stages$bottom,
                            include.lowest = TRUE,
                            labels = FALSE), 
         bin_max = 96 - cut(Max_Ma, breaks = dat_stages$top,
                            include.lowest = TRUE,
                            labels = FALSE)) %>% 
  # select entries, where the early and late interval fields indicate 
  # the same stg or the late_intervar field is empty
  filter(bin_min == bin_max | is.na(bin_max)) %>% 
  # in these entries, use the stg indicated by the early_interval
  select(bin = bin_min, everything(), -bin_max) %>% 
  # clean up
  drop_na(bin) %>%
  # get the age estimate for the corresponding bin
  left_join(dat_stages %>% 
              select(bin = stg))

# estimate sampling effort by the number of collections per bin
dat_sampling <- dat_occ_binned %>% 
  distinct(bin, collection_no) %>% 
  count(bin, 
        name = "shark_collections") %>% 
  # add bins with zero counts
  complete(bin = 69:94, fill = list(shark_collections = 0))

# alternatively the number of collections per bin from the pbdb

# download pbdb collections from bin 72 (Hauterivian) to bin 95 (Holocene)
# via the api call
dat_pbdb <- read_csv("https://paleobiodb.org/data1.2/colls/list.csv?interval=Hauterivian,Holocene")

# bin the data
dat_pbdb_sampling <- dat_pbdb %>% 
  mutate(bin_min = 95 - cut(min_ma, breaks = dat_stages$bottom,
                            include.lowest = TRUE,
                            labels = FALSE), 
         bin_max = 96 - cut(max_ma, breaks = dat_stages$top,
                            include.lowest = TRUE,
                            labels = FALSE)) %>% 
  # select entries, where the early and late interval fields indicate 
  # the same stg or the late_intervar field is empty
  filter(bin_min == bin_max | is.na(bin_max)) %>% 
  # in these entries, use the stg indicated by the early_interval
  select(bin = bin_min, everything(), -bin_max) %>% 
  # clean up
  drop_na(bin) %>% 
  # count collections per bin 
  distinct(bin, collection_no) %>% 
  count(bin, 
        name = "pbdb_collections") %>% 
  # add bins with zero counts
  complete(bin = 69:94, fill = list(pbdb_collections = 0))



# preservation potential --------------------------------------------------

# get preservation rate from the PyRate output
dat_preservation_raw <- read_delim(here("data",
                                        "raw",
                                        "fossil_occurrences",
                                        "combined_10_mcmc.log"))
# bring in right format
dat_preservation <- dat_preservation_raw %>% 
  select(mean_q, contains("TS")) %>% 
  pivot_longer(cols = contains("TS"), 
               names_to = "species", 
               values_to = "age") %>% 
  # bin to stages
  mutate(bin = 95 - cut(age, breaks = dat_stages$bottom,
                        include.lowest = TRUE,
                        labels = FALSE)+1, 
         .before = 1) %>% 
  # summarise
  drop_na(bin) %>% 
  group_by(bin, species) %>% 
  summarise(mean_q = mean(mean_q)) %>% 
  ungroup() %>% 
  # clean up species names for joining
  mutate(species = str_remove(species, "_TS"), 
         accepted_name = str_replace(species, "_", " ")) %>% 
  select(bin, accepted_name, mean_q) %>% 
  arrange(accepted_name) 


# merge and combine -------------------------------------------------------

# combine all datasets to one
dat_full <- list(dat_ext,
                 dat_preservation,
                 dat_sealevel,
                 dat_13C,
                 dat_SR,
                 dat_marine_units,
                 dat_outcrop,
                 dat_cont_area,
                 dat_temp,
                 dat_paleotemp,
                 dat_sampling,
                 dat_pbdb_sampling) %>% 
  reduce(full_join)  %>% 
  # remove missing values
  left_join(dat_preservation %>%
              group_by(accepted_name) %>%
              summarise(mean_q_av = mean(mean_q))) %>% 
  mutate(mean_q = if_else(is.na(mean_q), mean_q_av, mean_q)) %>% 
  drop_na(everything()) %>% 
  # get taxonomy
  left_join(dat_occurrences %>% 
              distinct(accepted_name, genus, family, order)) %>% 
  select(order, family, genus, species = accepted_name,
         everything())  


# scale predictors --------------------------------------------------------


# bring predictors on meaningful scales to enable comparisons
dat_merged <- dat_full %>%
  mutate(across(c(d13C, sr_mean, # scale all productivity parameters
                  n_units, outcrop_area, # same for outcrop parameters
                  mean_q, # preservation rate
                  pbdb_collections, shark_collections), #  outcrop parameters
                ~ scale(.x)[,1], 
                .names = "{col}_std")) %>% 
  select(bin, order, family, genus, species,
         ext_signal, 
         sea_level, 
         cont_area, 
         contains("temp"), 
         contains("std")) 

# save dataset
dat_merged %>% 
  write_rds(here("data",
                 "processed_merged_data.rds"))


