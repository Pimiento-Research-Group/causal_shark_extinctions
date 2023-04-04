library(iucnsim)
library(reticulate)
library(rredlist)
library(here)
library(tidyverse)


# plotting configurations
source(here("R", "config_file.R"))

# python function for simulations
reticulate::source_python(here("data", 
                               "iucn", 
                               "iucn_sim.py"))


# load data  ------------------------------------------------------------

# fossil data and environmental proxy data on species level
dat_merged <- read_rds(here("data",
                            "processed_merged_data.rds"))

# assign unique key
iucn_key = "1fecbeea639ba430f60510af483c5d4b282e3f097aa98e4613003c9903970df6"

# get IUCN history file
iucn_history_file <- get_iucn_history(reference_group = "Chondrichthyes",
                                      reference_rank = "class",
                                      iucn_key = iucn_key,
                                      outdir = here("data", 
                                                    "iucn"))

# get binomial species list
species_list <- read_rds(here("data",
                             "fossil_occurrences",
                             "database_occurrences_10_Jan_2023.rds")) %>% 
  filter(status == "extant", 
         rank == "species") %>% 
  distinct(modified_identified_name) %>% 
  mutate(nr_words = str_count(modified_identified_name, "\\W+")) %>% 
  filter(nr_words == 1) %>% 
  pull(modified_identified_name) 


# get most recent status for each taxon in target species list
extant_taxa_current_status <- get_most_recent_status_target_species(species_list = species_list,
                                                                   iucn_history_file =
                                                                     iucn_history_file,
                                                                   iucn_key =
                                                                     iucn_key,
                                                                   outdir =
                                                                     here("data",
                                                                          "iucn"))


dat_possible_extinct <- get_possibly_extinct_iucn_info(iucn_history_file,
                                                       outdir = here("data",
                                                                     "iucn"))

# estimate status transition rates
transition_rates_out <- estimate_transition_rates(extant_taxa_current_status,
                                                 iucn_history_file,
                                                 outdir = here("data",
                                                               "iucn"),
                                                 extinction_probs_mode = 0,
                                                 possibly_extinct_list =
                                                   dat_possible_extinct,
                                                 rate_samples = 100)


# simulate extinction times
future_sim_output <- run_future_sim(transition_rates_out,
                                    here("data",
                                         "iucn"),
                                   n_years = 200,
                                   n_sim = 1e4)

# extract extinction times
extinction_times <- future_sim_output[[1]]

# define mode
get_mode <- function(x) {
  ux <- na.omit(unique(x))
  ux[which.max(tabulate(match(x, ux)))]
}

# iterate through species and summarise simulations
dat_ext_times <- vector("double", nrow(extinction_times))

for (i in 1:nrow(extinction_times)) {
  
  dat_ext_times[[i]] <- extinction_times[[i]] %>%
    unlist() %>%
    get_mode()
  
}

# assign to dataframe
tibble(species = rownames(extinction_times), 
       ext_time = dat_ext_times) %>% 
  # round extinction times
  mutate(ext_time = round(ext_time, 0) %>% 
           as.integer()) %>% 
  # save results
  write_rds(here("data", 
                 "iucn_extinction_times.rds"))
