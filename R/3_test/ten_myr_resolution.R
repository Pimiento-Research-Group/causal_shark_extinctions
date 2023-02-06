# load packages
library(chronosphere)
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
  



# fossil data -------------------------------------------------------------


# load PyRate estimates
dat_species <- read_delim(here("data",
                               "fossil_occurrences",
                               "all_species_10_Grj_se_est_species_names.txt")) %>% 
  # estimate fad and lad for each species
  pivot_longer(cols = contains("ts"), 
               names_to = "origination", 
               values_to = "origination_age") %>% 
  pivot_longer(cols = contains("te"), 
               names_to = "extinction", 
               values_to = "extinction_age") %>% 
  group_by(species) %>% 
  summarise(ori_age = mean(origination_age),
            ext_age = mean(extinction_age)) %>% 
  # bin to 10myr
  mutate(bin_ori = cut(ori_age, breaks = seq(250, 0, by = -10),
                       include.lowest = TRUE,
                       labels = FALSE), 
         bin_ext = cut(ext_age, breaks = seq(250, 0, by = -10),
                       include.lowest = TRUE,
                       labels = FALSE)) %>% 
  drop_na(bin_ori, bin_ext) %>% 
  # fill in duration bins
  mutate(bin_occ = map2(.x = bin_ori, 
                        .y = bin_ext,
                        .f = ~ seq(.x, .y, by = -1))) %>% 
  select(species, bin_occ, bin_ext) %>% 
  unnest(bin_occ) %>% 
  # create extinction signal
  group_by(species) %>% 
  mutate(ext_signal = if_else(bin_occ == bin_ext, 1, 0)) %>% 
  # clean up
  ungroup() %>% 
  select(-bin_ext, 
         bin = bin_occ)



# load occurrence database
dat_occurrences <- read_rds(here("data",
                                 "fossil_occurrences",
                                 "database_occurrences_10_Jan_2023.rds"))

# bin the occurrences to 10 myr
dat_occ_binned <- dat_occurrences %>% 
  mutate(bin_min = cut(Min_Ma, breaks = seq(250, 0, by = -10),
                       include.lowest = TRUE,
                       labels = FALSE), 
         bin_max = cut(Max_Ma, breaks = seq(250, 0, by = -10),
                       include.lowest = TRUE,
                       labels = FALSE)) %>% 
  # select entries, where the early and late interval fields indicate 
  # the same stg or the late_interval field is empty
  filter(bin_min == bin_max | is.na(bin_max)) %>% 
  # in these entries, use the stg indicated by the early_interval
  select(bin = bin_min, everything(), -bin_max) %>% 
  # clean up
  drop_na(bin) %>% 
  # get the age estimate for the corresponding bin
  left_join(tibble(bin = 25:1, age = seq(245, 5, by = -10)))



# reconstruct paleo-coordinates
# reformat data
dat_coords <- dat_occ_binned %>% 
  distinct(Longitude, Latitude, age) %>% 
  # prepare for reconstruct function
  mutate(age = round(age)) %>% 
  select(Longitude, Latitude, age) %>% 
  as.data.frame() 


# preallocate vector
coord_vector <- data.frame("paleolong" = rep(NA, nrow(dat_coords)), 
                           "paleolat" = rep(NA, nrow(dat_coords)))

# iterate over rows
for (i in 1:nrow(dat_coords)) {
  
  long_lat <- dat_coords[i, 1:2]
  age <- dat_coords[i, 3]
  suppressMessages(
    tryCatch({
      coord_vector[i,] <- reconstruct(long_lat, age = age, chunk = 1)
    }, error = function(e){})
  )
  cat(paste0(i, "  "))
}

# merge back to original data
dat_occ_binned_coord <- dat_occ_binned %>% 
  distinct(Longitude, Latitude, age) %>% 
  select(Longitude, Latitude, age) %>% 
  bind_cols(coord_vector) %>% 
  full_join(dat_occ_binned) %>% 
  drop_na(paleolong, paleolat)  

# save data
dat_occ_binned_coord %>% 
  write_rds(here("data",
                 "fossil_occurrences",
                 "binned_rotated_occurrences_10myr.rds"))



# calculate latitude range in absolute degrees
dat_range_lat <- dat_occ_binned_coord %>% 
  # select only species
  filter(rank == "species") %>% 
  group_by(modified_identified_name) %>% 
  summarise(min_lat = min(paleolat), 
            max_lat = max(paleolat)) %>% 
  mutate(range_lat = abs(min_lat - max_lat)) %>% 
  select(modified_identified_name, range_lat)

# calculate great circle distance via geodesic distance between two 
# points specified by radian latitude/longitude using the
# Spherical Law of Cosines (slc)
dat_dist <- dat_occ_binned_coord %>% 
  # select only species
  filter(rank == "species") %>%
  group_by(modified_identified_name, bin) %>% 
  # get maximum distance
  summarise(min_lat = min(paleolat), 
            max_lat = max(paleolat), 
            min_long = min(paleolong), 
            max_long = max(paleolong)) %>% 
  ungroup() %>% 
  # convert to radians
  mutate(across(min_lat:max_long, ~ .x*pi/180)) %>% 
  # use earth radius and slc
  mutate(geo_dist = acos(sin(min_lat)*sin(max_lat) + cos(min_lat)*cos(max_lat) * cos(max_long-min_long)) * 6371) %>% 
  drop_na(geo_dist) %>% 
  select(modified_identified_name, bin, 
         geo_dist)


# calculate latitudinal preference of species
dat_latitude <- dat_occ_binned_coord %>% 
  # select only species
  filter(rank == "species") %>% 
  group_by(modified_identified_name, bin) %>% 
  summarise(latitude_pref = mean(paleolat)) %>% 
  ungroup()


# number of occurrences
dat_abund <- dat_occ_binned %>% 
  # select only species
  filter(rank == "species") %>% 
  count(bin, modified_identified_name) 

# number of species within genera
dat_abund_genus <- dat_occ_binned %>% 
  # select only species
  filter(rank == "species") %>% 
  count(bin, genus) 


# estimate sampling effort by the number of collections per bin
dat_sampling <- dat_occ_binned %>% 
  distinct(bin, collection_no) %>% 
  count(bin) %>% 
  # add bins with zero counts
  mutate(bin = factor(bin, levels = 1:14)) %>% 
  complete(bin, fill = list(n = 0))

# alternatively the number of collections per bin from the pbdb

# download pbdb collections from bin 72 (Hauterivian) to bin 95 (Holocene)
# via the api call
dat_pbdb <- read_csv("https://paleobiodb.org/data1.2/colls/list.csv?interval=Hauterivian,Holocene")

# bin the data
dat_pbdb_sampling <- dat_pbdb %>% 
  mutate(bin_min = cut(min_ma, breaks = seq(250, 0, by = -10),
                       include.lowest = TRUE,
                       labels = FALSE), 
         bin_max = cut(max_ma, breaks = seq(250, 0, by = -10),
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
  count(bin) %>% 
  # add bins with zero counts
  mutate(bin = factor(bin, levels = 1:14)) %>% 
  complete(bin, fill = list(n = 0))



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
  # bin to 10 myr
  mutate(bin = cut(age, breaks = seq(250, 0, by = -10),
                   include.lowest = TRUE,
                   labels = FALSE)-1, 
         .before = 1) %>% 
  # summarise
  drop_na(bin) %>% 
  group_by(bin, species) %>% 
  summarise(mean_q = mean(mean_q)) %>% 
  ungroup() %>% 
  # clean up species names for joining
  mutate(species = str_remove(species, "_TS"), 
         modified_identified_name = str_replace(species, "_", " ")) %>% 
  select(bin, modified_identified_name, mean_q) %>% 
  arrange(modified_identified_name)


  
# merge and combine -------------------------------------------------------

# combine all datasets to one
dat_full <- dat_range_lat %>% 
  # geographic range
  full_join(dat_dist) %>% 
  full_join(dat_abund) %>% 
  # latitudinal preference
  full_join(dat_latitude) %>% 
  # abundance
  rename(abund = n) %>% 
  # preservation rate
  left_join(dat_preservation) %>% 
  # sampling
  full_join(dat_pbdb_sampling %>% 
              mutate(bin = as.numeric(as.character(bin)))) %>% 
  rename(pbdb_collections = n) %>% 
  full_join(dat_sampling %>% 
              mutate(bin = as.numeric(as.character(bin))) %>% 
              rename(shark_collections = n)) 

# add species extinction signal
dat_full %>% 
  left_join(dat_species %>%
              mutate(modified_identified_name = str_replace(species, "_", " ")) %>%
              select(-species) %>% 
              filter(modified_identified_name %in% dat_full$modified_identified_name)) %>%
  # make missing values explicit
  drop_na(ext_signal) %>% 
  replace_na(list(geo_dist = 0, range_lat = 0,  
                  latitude_pref = mean(dat_latitude$latitude_pref))) %>%
  left_join(dat_preservation %>%
              group_by(bin) %>%
              summarise(mean_q_av = mean(mean_q))) %>% 
  mutate(mean_q = if_else(is.na(mean_q), mean_q_av, mean_q)) %>% 
  # get taxonomy
  left_join(dat_occ_binned %>% 
              distinct(modified_identified_name, genus, family, order)) %>% 
  # add genus counts
  left_join(dat_abund_genus %>% rename(n_genus = n)) %>% 
  select(order, family, genus, species = modified_identified_name,
         everything()) %>% 
  # save dataset
  write_rds(here("data",
                 "processed_fossil_data_10myr.rds"))







