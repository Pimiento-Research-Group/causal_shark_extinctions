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
         abund, 
         ext_signal, 
         sea_level, 
         cont_area, 
         contains("temp"), 
         contains("std")) 


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
dat_test <- tibble(cor_val = vector("double", # Student's t-test statistic
                                    length = length(implied_conditions)), 
                   pval = vector("double", # the p-value, assuming a two-sided alternative
                                 length = length(implied_conditions)))

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
  
  
  test.output <- pcor.test(cor.output, 
                           length(implied_conditions[[i]]$Z), 
                           n = nrow(dat_dag))
  
  dat_test[i, "cor_val"] <- cor.output
  
  dat_test[i, "pval"] <- pluck(test.output, "pvalue")
  
}

dat_test %>% 
  ggplot(aes(cor_val)) +
  geom_density()

dat_test <- dat_test %>% 
  # apply bonferroni correction 
  mutate(pval_cor = pval * length(implied_conditions), 
         p_sign = if_else(pval_cor <= 0.05, "sign", "non-sign")) 

implied_conditions[dat_test[, "p_sign"] == "sign"]

implied_conditions[abs(dat_test[, "cor_val"]) >= 0.5]
