# load packages
library(tidyverse)
library(here)
library(readxl)


# plotting configurations
source(here("R", "config_file.R"))


# stage data --------------------------------------------------------------

# load 95 bin Phanerozoic time scale based on the stratigraphic stages 
# of Gradstein et al. 2020
load(here("data", 
          "stages.Rdata"))

# clean up for binning
dat_stages <- stages %>% 
  as_tibble() %>% 
  select(stg, bottom, top)



### Phanerozoic deep-ocean and global average temperature from Scotese et al 2021 ###
# load data
dat_temp <- read_xlsx(here("data",
                           "raw",
                           "temperature",
                           "scotese_et_al_2021.xlsx")) %>% 
  # clean up colnames
  select(age = Age,
         temp_gat = GAT,
         temp_deep = "Deep Ocean") %>% 
  # bin to stages
  mutate(bin = 95 - cut(age, breaks = dat_stages$bottom,
                        include.lowest = TRUE,
                        labels = FALSE)) %>% 
  drop_na(bin) %>% 
  group_by(bin) %>% 
  summarise(temp_gat_binned = mean(temp_gat), 
            temp_deep_binned = mean(temp_deep)) %>% 
  ungroup() %>% 
  select(bin, temp_gat_binned, temp_deep_binned) %>% 
  filter(between(bin, 69, 94)) %>% 
  complete(bin = 69:94) %>%
  fill(c(temp_gat_binned, temp_deep_binned), .direction = "downup")


# extinction signal -------------------------------------------------------

# get those species that survive until now
spec_modern <- read_delim(here("data",
                               "fossil_occurrences",
                               "combined_10_se_est_species_names.txt")) %>% 
  filter(te == 0) %>% 
  pull(species)

# load PyRate estimates
dat_pyr <- read_delim(here("data",
                           "fossil_occurrences",
                           "combined_10_se_est_species_names.txt")) %>% 
  # estimate fad and lad for each species
  mutate(across(c(ts, te), abs)) %>% 
  group_by(species) %>% 
  summarise(ori_age = mean(ts),
            ext_age = mean(te)) %>% 
  # bin fad and lad to stages
  mutate(bin_ori = 95 - cut(ori_age, breaks = dat_stages$bottom,
                            include.lowest = TRUE,
                            labels = FALSE), 
         bin_ext = 96 - cut(ext_age, breaks = dat_stages$top,
                            include.lowest = TRUE,
                            labels = FALSE)) %>% 
  drop_na(bin_ori, bin_ext) %>% 
  select(-contains("age")) %>% 
  filter(!species %in% spec_modern)


# combine -----------------------------------------------------------------

# merge
dat_pyr %>% 
  left_join(dat_temp %>% 
              select(bin_ori = bin, 
                     ori_temp = temp_gat_binned)) %>% 
  left_join(dat_temp %>% 
              select(bin_ext = bin, 
                     ext_temp = temp_gat_binned)) %>% 
  mutate(temp_bin_ori = 30 - cut(ori_temp, breaks = 15:30, 
                            labels = FALSE), 
         temp_ext_ori = 30 - cut(ext_temp, breaks = 15:30, 
                                 labels = FALSE)) %>% 
  drop_na() %>% 
  count(temp_bin_ori, temp_ext_ori) %>% 
  ggplot(aes(temp_bin_ori, temp_ext_ori)) +
  geom_point(aes(size = n), 
             shape = 21, 
             position = position_jitter(width = 0.5, 
                                        height = 0.5)) +
  geom_smooth(method = "lm") +
  coord_cartesian(xlim = c(15, 30), 
                  ylim = c(15, 30)) +
  labs(x = "Temperature at Origination", 
       y = "Temperature at Extinction")
