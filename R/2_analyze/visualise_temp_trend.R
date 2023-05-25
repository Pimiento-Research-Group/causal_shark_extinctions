# load packages
library(tidyverse)
library(here)


# plotting configurations
source(here("R", "config_file.R"))



# load data  ------------------------------------------------------------

# get predictions from different scales/ data sets
dat_merged <- read_rds(here("data",
                            "processed_merged_data.rds"))