# load packages
library(tidyverse)
library(here)


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
  mutate(across(c(range_lat, geo_dist, # scale geographic range parameters
                  mean_q, # and preservation rate
                  pbdb_collections, shark_collections), # and sampling effort
                ~ scale(.x)[,1], 
                .names = "{col}_std"), 
         latitude_pref_abs = abs(latitude_pref)) # use absolute latitude 



# combine data sources ----------------------------------------------------

# merge
left_join(dat_proxy, 
          dat_fossil, 
          by = "bin") %>%
  drop_na(species) %>% 
  replace_na(list(abund = 0)) %>% 
  select(bin, order, family, genus, species,
         ext_signal, 
         sea_level, 
         cont_area, 
         contains("temp"), 
         contains("std")) %>%
  ggplot(aes(x = temp_deep_binned, y = ext_signal)) +
  geom_point() +
  stat_smooth(method = "glm",
              method.args = list(family = binomial))
  ggplot(aes(temp_gat_binned)) +
  geom_density()

