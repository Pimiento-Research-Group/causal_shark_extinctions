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


# reconstruct paleo-coordinates
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

