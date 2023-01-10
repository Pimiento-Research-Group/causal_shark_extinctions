# load packages
library(tidyverse)
library(here)


# plotting configurations
source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

# species level extinction signal
dat_species <- read_rds(here("data",
                             "species_extinction_signal.rds"))

# environmental data on stage level spanning the full time range
# set up function
format_proxy_data <- function(data_file, 
                              column_name) {
  read_rds(here("data",
                "proxy_data", 
                data_file)) %>% 
    group_by(bin) %>% 
    summarise({{ column_name }} := mean({{ column_name }})) %>%
    arrange(bin) %>%
    filter(bin >= min(dat_species$bin)) %>%
    complete(bin = 69:95) %>%
    fill({{ column_name }}, .direction = "downup")
  
}

# apply function
dat_proxy <- map2(
  .x = c(
    "13C_full.rds",
    "cont_area_full.rds",
    "diatom_full.rds",
    "marine_units_full.rds",
    "sealevel_full.rds",
    "SR_full.rds",
    "temp_full.rds",
    "temp_full.rds"
  ),
  .y = c(
    as.symbol("d13C"),
    as.symbol("area"),
    as.symbol("div_mean"),
    as.symbol("n_units"),
    as.symbol("sea_level"),
    as.symbol("sr_value"),
    as.symbol("temp_gat"),
    as.symbol("temp_deep")
  ),
  .f = ~ format_proxy_data(data_file = .x ,
                           column_name = {{ .y }})) %>%
  reduce(full_join) 


lm(area ~ temp_gat + sea_level, data = dat_proxy) %>% 
  summary()

  
read_rds(here("data",
              "proxy_data",
              "temp_full.rds")) %>%
  group_by(bin) %>% 
  summarise(n = mean(n)) %>% 
  arrange(bin) %>% 
  filter(bin >= min(dat_species$bin)) %>% 
  complete(bin = 69:95) %>% 
  fill(n, .direction = "downup") 
  
  


list.files(here("data", 
                "proxy_data"), 
           pattern = "*full*") 
  