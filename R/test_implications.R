# load packages
library(tidyverse)
library(here)


# plotting configurations
source(here("R", "config_file.R"))



# species level extinction signal -----------------------------------------

# load data 
dat_species <- read_rds(here("data",
                             "species_extinction_signal.rds"))



# proxy data --------------------------------------------------------------


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
    "marine_units_full.rds",
    "outcrop_full.rds",
    "sealevel_full.rds",
    "SR_full.rds",
    "temp_full.rds",
    "temp_full.rds"
  ),
  .y = c(
    as.symbol("d13C"),
    as.symbol("cont_area"),
    as.symbol("n_units"),
    as.symbol("outcrop_area"),
    as.symbol("sea_level"),
    as.symbol("sr_value"),
    as.symbol("temp_gat"),
    as.symbol("temp_deep")
  ),
  .f = ~ format_proxy_data(data_file = .x ,
                           column_name = {{ .y }})) %>%
  reduce(full_join) 




# scale predictors --------------------------------------------------------

dat_proxy %>% 
  mutate(n_units_log = log(n_units),
         outcrop_area_std = scale(outcrop_area)[,1]) %>% 
  ggplot(aes(outcrop_area_std)) +
  geom_density()


lm(area ~ temp_gat + sea_level, data = dat_proxy) %>% 
  confint()

  
read_rds(here("data",
              "proxy_data",
              "outcrop_full.rds")) %>%
  group_by(bin) %>% 
  summarise(area = mean(area)) %>% 
  arrange(bin) %>% 
  filter(bin >= min(dat_species$bin)) %>% 
  complete(bin = 69:95) %>% 
  fill(area, .direction = "downup") 
  
  


list.files(here("data", 
                "proxy_data"), 
           pattern = "*full*") 
  