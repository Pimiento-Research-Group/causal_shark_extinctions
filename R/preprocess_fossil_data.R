# load packages
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
  select(stg, bottom, top)

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

list.files(here("data"), pattern = "*full*")

divDyn:::stag