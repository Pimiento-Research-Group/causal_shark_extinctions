# load packages
library(tidyverse)
library(here)


# plotting configurations
source(here("R", "config_file.R"))


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
            ext_age = mean(extinction_age))
