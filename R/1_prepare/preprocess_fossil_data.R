# load packages
library(chronosphere)
library(tidyverse)
library(here)


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
  select(stg, bottom, mid, top)

# load species data ---------------------------------------------------------------


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
  # bin fad and lad to stages
  mutate(bin_ori = 95 - cut(ori_age, breaks = dat_stages$bottom,
                            include.lowest = TRUE,
                            labels = FALSE), 
         bin_ext = 95 - cut(ext_age, breaks = dat_stages$bottom,
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
  mutate(ext_signal = if_else(bin_occ == bin_ext, 1, 0)) %>% 
  # clean up
  ungroup() %>% 
  select(-bin_ext, 
         bin = bin_occ)

# save data
dat_species %>% 
  write_rds(here("data", 
                 "species_extinction_signal.rds"))



# occurrence data ---------------------------------------------------------

# load occurrence database
dat_occurrences <- read_rds(here("data",
                                 "fossil_occurrences",
                                 "database_occurrences_10_Jan_2023.rds"))

# bin the occurrences to stages
dat_occ_binned <- dat_occurrences %>% 
  mutate(bin_min = 95 - cut(Min_Ma, breaks = dat_stages$bottom,
                        include.lowest = TRUE,
                        labels = FALSE), 
         bin_max = 95 - cut(Max_Ma, breaks = dat_stages$bottom,
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
              select(bin = stg, age = mid))



# reconstruct paleo-coordinates -------------------------------------------

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
                 "binned_rotated_occurrences.rds"))
  


# geographic range --------------------------------------------------------


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
  group_by(modified_identified_name) %>% 
  # get maximum distance
  summarise(min_lat = min(paleolat), 
            max_lat = max(paleolat), 
            min_long = min(paleolong), 
            max_long = max(paleolong)) %>% 
  # convert to radians
  mutate(across(min_lat:max_long, ~ .x*pi/180)) %>% 
  # use earth radius and slc
  mutate(geo_dist = acos(sin(min_lat)*sin(max_lat) + cos(min_lat)*cos(max_lat) * cos(max_long-min_long)) * 6371) %>% 
  drop_na(geo_dist) %>% 
  select(modified_identified_name, 
         geo_dist)



# latitude ----------------------------------------------------------------

# calculate latitudinal preference of species
dat_latitude <- dat_occ_binned_coord %>% 
  # select only species
  filter(rank == "species") %>% 
  group_by(modified_identified_name) %>% 
  summarise(latitude_pref = mean(paleolat))
  

# calculate abundance -----------------------------------------------------

dat_abund <- dat_occ_binned %>% 
  # select only species
  filter(rank == "species") %>% 
  count(bin, modified_identified_name) 
  


# sampling effort ---------------------------------------------------------

# estimate sampling effort by the number of collections per bin
dat_sampling <- dat_occ_binned %>% 
  distinct(bin, collection_no) %>% 
  count(bin) %>% 
  # add bins with zero counts
  mutate(bin = factor(bin, levels = 95:69)) %>% 
  complete(bin, fill = list(n = 0))

# alternatively the number of collections per bin from the pbdb

# download pbdb collections from bin 72 (Hauterivian) to bin 95 (Holocene)
# via the api call
dat_pbdb <- read_csv("https://paleobiodb.org/data1.2/colls/list.csv?interval=Hauterivian,Holocene")

# bin the data
dat_pbdb_sampling <- dat_pbdb %>% 
  mutate(bin_min = 95 - cut(min_ma, breaks = dat_stages$bottom,
                          include.lowest = TRUE,
                          labels = FALSE), 
       bin_max = 95 - cut(max_ma, breaks = dat_stages$bottom,
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
  mutate(bin = factor(bin, levels = 95:69)) %>% 
  complete(bin, fill = list(n = 0))



# preservation potential --------------------------------------------------

# get preservation rate from the PyRate output
dat_preservation <- read_rds(here("data", 
                              "fossil_occurrences", 
                              "preservation_rate_per_species.rds")) %>% 
  # use only occurrences from species that could get binned to stages
  filter(modified_identified_name %in% dat_abund$modified_identified_name)



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
  full_join(dat_preservation) %>% 
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
  select(order, family, genus, species = modified_identified_name,
         everything()) %>% 
  # save dataset
  write_rds(here("data",
                 "processed_fossil_data.rds"))

  
