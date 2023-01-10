# load packages
library(tidyverse)
library(here)


# plotting configurations
source(here("R", "config_file.R"))


# load data ---------------------------------------------------------------

# species level extinction signal
dat_species <- read_rds(here("data",
                             "species_extinction_signal.rds"))
