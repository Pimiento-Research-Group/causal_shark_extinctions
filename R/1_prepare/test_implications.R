# load packages
library(tidyverse)
library(here)
library(dagitty)
library(ggm)

# plotting configurations
source(here("R", "config_file.R"))



# load data  ------------------------------------------------------------


# save dataset
dat_merged <- read_rds(here("data", 
                 "processed_merged_data.rds"))

# change column names to fit with the directed acyclic graph
dat_dag <- dat_merged %>% 
  transmute(productivity = d13C_std,
            "outcrop area" = outcrop_area_std,
            "sea level" = sea_level,
            "shelf area" = cont_area,
            "temperature" = temp_deep_binned,
            paleotemperature = temp_deep_lt1, 
            "extinction risk" = ext_signal, 
            "taxonomic identity" = as.integer(as.factor(genus)), 
            "sampling effort" = shark_collections_std, 
            "preservation potential" = mean_q_std) %>% 
  drop_na(everything())

# create subsets for bootstrapping
dat_dag_list <- replicate(5, slice_sample(dat_dag,
                                          n = nrow(dat_dag),
                                          replace = TRUE),
                          simplify = FALSE)



# directed acyclic graph ---------------------------------------------

# load the graph 
dag <- downloadGraph("dagitty.net/mjiV5Qf")

# get testable implications of that model
implied_conditions <- impliedConditionalIndependencies(dag)

# set up function to get partial correlation estimates
part_cor <- function(data_set) {
  
  # vector to save results into 
  cor_val <- vector("double", length = length(implied_conditions))
  
  for (i in 1:length(implied_conditions)) {
    
    # simple correlation case
    if (length(implied_conditions[[i]]$Z) == 0) {
      
      cor.output <- pcor(c(implied_conditions[[i]]$X,
                           implied_conditions[[i]]$Y), 
                         var(data_set))
      
    # partical correlation case  
    } else {
      
      cor.output <- pcor(c(implied_conditions[[i]]$X,
                           implied_conditions[[i]]$Y,
                           implied_conditions[[i]]$Z), 
                         var(data_set)) 
      
    }
    
    cor_val[i] <- cor.output
    
  }
  
  enframe(cor_val)
  
}

dat_dag_list %>% 
  map(part_cor) %>% 
  reduce(full_join) %>% 
  group_by(name) %>% 
  summarise(mean_cl_normal(value)) 

implied_conditions[abs(cor_val) >= 0.3]

