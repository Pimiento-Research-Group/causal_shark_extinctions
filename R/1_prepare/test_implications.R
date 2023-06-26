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
dat_dag_list <- replicate(2, slice_sample(dat_dag,
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



# visualise ---------------------------------------------------------------


plot_cor <- dat_dag_list %>% 
  map(part_cor) %>% 
  reduce(full_join) %>% 
  group_by(name) %>% 
  summarise(mean_cl_normal(value)) %>% 
  add_column(imp_cond = implied_conditions %>%
               map_chr(as.character)) %>% 
  # use the same notation as in Figure 1
  mutate(imp_cond = str_replace_all(imp_cond, "txni", "I"), 
         imp_cond = str_replace_all(imp_cond, "prsp", "P"),
         imp_cond = str_replace_all(imp_cond, "smpe", "Se"),
         imp_cond = str_replace_all(imp_cond, "prdc", "N"),
         imp_cond = str_replace_all(imp_cond, "shla", "As"),
         imp_cond = str_replace_all(imp_cond, "otca", "Ao"),
         imp_cond = str_replace_all(imp_cond, "slvl", "L"),
         imp_cond = str_replace_all(imp_cond, "tmpr", "T"),
         imp_cond = str_replace_all(imp_cond, "pltm", "Tp"),
         imp_cond = str_replace_all(imp_cond, "extr", "E"),
         imp_cond = fct_reorder(imp_cond, desc(y))) %>% 
  ggplot() +
  geom_vline(xintercept = 0, 
             linetype = "dotted", 
             colour = "grey20") +
  geom_vline(xintercept = c(-0.3, 0.3), 
             linetype = "dashed", 
             colour = "coral3") +
  geom_linerange(aes(y, imp_cond, 
                     xmin = ymin, xmax = ymax)) +
  geom_point(aes(y, imp_cond), 
             size = 4, 
             shape = 21, 
             fill = "white", 
             colour = "grey10") +
  labs(y = NULL, 
       x = "Partial correlation")



# save plot
ggsave(plot_cor,
       filename = here("figures",
                       "supplement",
                       "conditional_dependencies.png"),
       width = image_width, height = image_height*1.5,
       units = image_units,
       bg = "white", device = ragg::agg_png)


