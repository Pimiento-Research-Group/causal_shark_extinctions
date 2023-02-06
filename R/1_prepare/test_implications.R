# load packages
library(tidyverse)
library(here)
library(dagitty)
library(ggm)

# plotting configurations
source(here("R", "config_file.R"))



# load data  ------------------------------------------------------------

# fossil data on species level
dat_fossil <- read_rds(here("data",
                            "processed_fossil_data.rds"))

# environmental proxy data
dat_proxy <- read_rds(here("data",
                           "processed_proxy_data.rds"))



# scale predictors --------------------------------------------------------

# bring predictors on meaningful scales to enable comparisons
dat_proxy <- dat_proxy %>% 
  mutate(across(c(d13C, sr_value, # scale all productivity parameters
                  n_units, outcrop_area), # same for outcrop parameters
                ~ scale(.x)[,1], 
                .names = "{col}_std"))

# same for fossil data
dat_fossil <- dat_fossil %>% 
  # add zero genus counts
  replace_na(list(n_genus = 0)) %>% 
  mutate(across(c(range_lat, geo_dist, # scale geographic range parameters
                  mean_q, # and preservation rate
                  pbdb_collections, shark_collections), # and sampling effort
                ~ scale(.x)[,1], 
                .names = "{col}_std"), 
         latitude_pref_abs = abs(latitude_pref)) # use absolute latitude 



# combine data sources ----------------------------------------------------

# merge
dat_merged <- left_join(dat_proxy,
                        dat_fossil,
                        by = "bin") %>% 
  drop_na(species) %>% 
  replace_na(list(abund = 0)) %>% 
  select(bin, order, family, genus, species,
         latitude_pref_abs,
         n_genus,
         abund, 
         ext_signal, 
         sea_level, 
         cont_area, 
         contains("temp"), 
         contains("std")) 

# save dataset
dat_merged %>% 
  write_rds(here("data", 
                 "processed_merged_data.rds"))

# change column names to fit with the directed acyclic graph
dat_dag <- dat_merged %>% 
  transmute(latitude = latitude_pref_abs,
            productivity = d13C_std,
            "outcrop area" = outcrop_area_std,
            "sea level" = sea_level,
            "shelf area" = cont_area,
            "geographic range" = geo_dist_std,
            "temperature" = temp_deep_binned,
            paleotemperature = temp_deep_lt1, 
            "extinction risk" = ext_signal, 
            "taxonomic identity" = as.integer(as.factor(genus)), 
            abundance = abund, 
            "sampling effort" = shark_collections_std, 
            "preservation potential" = mean_q_std)



# directed acyclic graph ---------------------------------------------

# load the graph 
dag <- downloadGraph("dagitty.net/m_UM7hV")

# get testable implications of that model
implied_conditions <- impliedConditionalIndependencies(dag)

# set up dataframe to save results of tests
cor_val <- vector("double", length = length(implied_conditions))

for (i in 1:length(implied_conditions)) {
  
  if (length(implied_conditions[[i]]$Z) == 0) {
    
    cor.output <- pcor(c(implied_conditions[[i]]$X,
                         implied_conditions[[i]]$Y), 
                       var(dat_dag))
    
  } else {
    
    cor.output <- pcor(c(implied_conditions[[i]]$X,
                         implied_conditions[[i]]$Y,
                         implied_conditions[[i]]$Z), 
                       var(dat_dag)) 
    
  }
  
  cor_val[i] <- cor.output
  
}

implied_conditions[cor_val >= 0.5]
