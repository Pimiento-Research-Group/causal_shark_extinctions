# load packages
library(tidyverse)
library(here)
library(readxl)
library(deeptime)

# plotting configurations
source(here("R", "config_file.R"))


# sea level ---------------------------------------------------------------

### Jurassic to recent sea level from miller et al 2005 ###
# load data
dat_sealevel_full <- read_xlsx(here("data",
                                    "raw",
                                    "sea_level",
                                    "miller_et_al_2005.xlsx")) %>% 
  # clean up column names
  rename(sea_level = "Sea level best estimate", 
         age = "Age for best estimate (MA)")  %>% 
  # calculate 1 myr moving average
  mutate(bin = cut(age, breaks = seq(180, 0, by = -1), 
                   include.lowest = TRUE)) %>% 
  group_by(bin) %>% 
  mutate(sea_level_mean = mean(sea_level)) %>% 
  ungroup()

# save
dat_sealevel_full %>% 
  write_rds(here("data", 
                 "proxy_data",
                 "sealevel_full.rds"))

# visualize
plot_sealevel_full <- dat_sealevel_full %>%
  ggplot(aes(age, sea_level)) +
  geom_hline(yintercept = 0, 
             colour = "grey70", 
             linetype = "dashed") +
  geom_point(shape = 21, 
             alpha = 0.3, 
             colour = "grey10", 
             fill = "grey50") +
  geom_line(aes(y = sea_level_mean), 
            colour = "#2c7c94", 
            linewidth = 1.3) +
  scale_x_reverse() +
  coord_geo(xlim = c(180, 0), 
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Global mean sea level [m]", 
       title = "Miller et al. 2005", 
       subtitle = "1 myr moving average")


  
### cenozoic sea level from miller et al 2020 ###
dat_sealevel_ceno <- read_xlsx(here("data",
                                    "raw",
                                    "sea_level",
                                    "miller_et_al_2020.xlsx")) %>% 
  # clean up column names
  rename(sea_level = "Sea Level (m) Smoothed", 
         age = "Age (ka)")  %>% 
  # transform to million years
  mutate(age = age/1000) 

# save
dat_sealevel_ceno %>% 
  write_rds(here("data", 
                 "proxy_data",
                 "sealevel_ceno.rds"))

# visualize
plot_sealevel_ceno <- dat_sealevel_ceno %>%
  ggplot(aes(age, sea_level)) +
  geom_hline(yintercept = 0, 
             colour = "grey70", 
             linetype = "dashed") +
  geom_line(colour = "#2c7c94", 
            linewidth = 1.3) +
  scale_x_reverse() +
  coord_geo(xlim = c(65, 0), 
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Global mean sea level [m]", 
       title = "Miller et al. 2020", 
       subtitle = "Smoothed to 20ka") 



# productivity ------------------------------------------------------------


### Cenozoic delta13C values from Westerhold et al 2020 ###
# load data
dat_13C_ceno <- read_xlsx(here("data",
                                "raw",
                                "productivity",
                                "westerhold_et_al_2020.xlsx"),
                           sheet = "Table S34",
                           skip = 1) %>% 
  # remove redundant columns
  select(age = age_tuned, fine_loess = ISOBENd13cLOESSsmooth, 
         coarse_loess = ISOBENd13cLOESSsmoothLongTerm)

# save
dat_13C_ceno %>% 
  write_rds(here("data", 
                 "proxy_data",
                 "d13C_ceno.rds"))

# visualize
plot_13C_ceno <- dat_13C_ceno %>%
  ggplot(aes(age, fine_loess)) +
  geom_hline(yintercept = 0, 
             colour = "grey70", 
             linetype = "dashed") +
  geom_line(colour = "#c4aa23", 
            linewidth = 0.8) +
  geom_line(aes(y = coarse_loess), 
            colour = "#fbe45b", 
            linewidth = 1.3) +
  scale_x_reverse() +
  coord_geo(xlim = c(70, 0), 
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Benthic d13C [‰]", 
       title = "Westerhold et al. 2020", 
       subtitle = "Smoothed to 20ka (dark yellow) and 1myr (light yellow)") 


### Phanerozoic delta13C values from Veizer and Prokoph ###
# load data
dat_13C_full <- read_csv(here("data",
                               "raw",
                               "productivity",
                               "veizer_prokoph_2015.csv")) %>% 
  # clean up column names
  select(age = "gts2012", 
         d13C)  %>% 
  filter(age <= 180) %>% 
  drop_na(age, d13C) %>% 
  # calculate 5 myr moving average
  mutate(bin = cut(age, breaks = seq(180, -5, by = -5))) %>% 
  group_by(bin) %>% 
  mutate(d13C_mean = mean(d13C)) %>% 
  ungroup()

# save
dat_13C_full %>% 
  write_rds(here("data", 
                 "proxy_data",
                 "13C_full.rds"))

# visualize
plot_13C_full <- dat_13C_full %>%
  ggplot(aes(age, d13C)) +
  geom_hline(yintercept = 0, 
             colour = "grey70", 
             linetype = "dashed") +
  geom_point(shape = 21, 
             alpha = 0.3, 
             colour = "grey10", 
             fill = "grey50") +
  geom_line(aes(y = d13C_mean), 
            colour = "#fbe45b", 
            linewidth = 1.3) +
  scale_x_reverse() +
  coord_geo(xlim = c(180, 0), 
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Benthic d13C [‰]", 
       title = "Veizer and Prokoph 2015", 
       subtitle = "5 myr moving average") 



### Phanerozoic 87Sr/86Sr values from McArthur, Howarth, Shields 2012 ###
# load data
dat_SR_full <- read_csv(here("data",
                              "raw",
                              "productivity",
                              "mcarthur_howarth_shields_2012.csv"), 
                        col_names = FALSE) %>% 
  # clean up column names
  select(age = X1, 
         sr_value = X2) 

# save
dat_SR_full %>% 
  write_rds(here("data", 
                 "proxy_data",
                 "SR_full.rds"))

# visualize
plot_SR_full <- dat_SR_full %>%
  ggplot(aes(age, sr_value)) +
  geom_hline(yintercept = 0, 
             colour = "grey70", 
             linetype = "dashed") +
  geom_line(colour = "#fbe45b", 
            linewidth = 1.3) +
  scale_x_reverse() +
  coord_geo(xlim = c(140, 0), 
            ylim = c(0.7060, 0.7095),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "87Sr / 86Sr", 
       title = "Mcarthur, Howarth, Shields 2012", 
       subtitle = "Lowess fit curve") 


### Phanerozoic diatom diversity from the Neptune database ###
# load data
dat_diatom_full <- read_tsv(here("data",
                                 "raw",
                                 "productivity",
                                 "diatom_2022-11-27_09-58-09.csv")) %>% 
  # clean up column names
  select(age = "Age (Ma) Gradstein et al. 2012", 
         taxon_id = "Resolved Taxon ID") 

# rarefaction based on this paper: https://www.nature.com/articles/nature07435
# number of rarefaction iteration
nr_iter <- 1000
# number of samples
nr_samples <- 96
# pre-allocate list
dat_diatom_list <- vector(mode = "list", 
                          length = nr_iter)
# iterate
for (i in 1:nr_iter) {
  
  dat_diatom_list[[i]] <- dat_diatom_full %>% 
    # bin data
    mutate(bin = cut(age, breaks = seq(90, 0, by = -1))) %>% 
    group_by(bin) %>% 
    # rarefaction
    slice_sample(n = 96, replace = TRUE) %>% 
    # calculate diversity
    summarise(count = n_distinct(taxon_id)) %>% 
    # add bins with zero counts
    complete(bin, fill = list(count = 0)) %>% 
    ungroup() %>% 
    # add numeric age for plotting
    add_column(age = seq(0.5, 89.5), 
               rare_id = i)
  
}

# summarize
dat_diatom_full <- dat_diatom_list %>% 
  bind_rows() %>% 
  group_by(bin, age) %>% 
  summarise(div_mean = mean(count), 
            div_sd = sd(count)) %>% 
  mutate(div_low = div_mean - 1.96*div_sd, 
         div_high = div_mean + 1.96*div_sd)

# save
dat_diatom_full %>% 
  write_rds(here("data", 
                 "proxy_data",
                 "diatom_full.rds"))

# visualize
plot_diatom_full <- dat_diatom_full %>%
  ggplot(aes(age, div_mean, 
             ymin = div_low, 
             ymax = div_high)) +
  geom_ribbon(fill = "#c4aa23", 
              alpha = 0.6) +
  geom_line(colour = "#fbe45b", 
            linewidth = 1.3) +
  scale_x_reverse() +
  coord_geo(xlim = c(70, 0), 
            ylim = c(-1, 70),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Diatom diversity", 
       title = "Neptune database", 
       subtitle = "Rarefied (n = 96) and cleaned data,\nbinned to one myr") 



# outcrop area ------------------------------------------------------------


### Phanerozoic rock units from the macrostrat database ###
# load data
dat_outcrop_full <- read_csv(here("data",
                                  "raw",
                                  "outcrop_area",
                                  "macrostrat_24_11_2022.csv")) %>% 
  # select only marine environments
  filter(!str_detect(.$environ, "non-marine")) %>% 
  mutate(mean_age = (t_age+b_age)/2) %>% 
  select(unit_id, t_age, b_age, mean_age, environ) 

# save
dat_outcrop_full %>% 
  write_rds(here("data", 
                 "proxy_data",
                 "outcrop_full.rds"))

# visualize
plot_outcrop_full <- dat_outcrop_full %>%
  # bin to 5 myr
  mutate(bin = cut(mean_age, breaks = seq(170, 0, by = -2))) %>% 
  count(bin) %>% 
  # add bins with zero counts
  complete(bin, fill = list(n = 0)) %>% 
  mutate(age = seq(1, 171, by = 2)) %>% 
  ggplot(aes(age, n)) +
  geom_line(linewidth = 1.5, colour = "#a6d0c8") +
  scale_x_reverse() +
  coord_geo(xlim = c(160, 0), 
            ylim = c(0, 3200),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Marine rock units", 
       title = "Macrostrat", 
       subtitle = "Binned into 2 myr") 




### Phanerozoic outcrop area from Wall, Ivany, Wilkinson 2009 ###
# load data
dat_outcrop_full <- read_xlsx(here("data",
                                  "raw",
                                  "outcrop_area",
                                  "wall_ivany_wilkinson_2009.xlsx")) %>% 
  # clean up colnames
  rename(age = "age (ma)", 
         area = "cumul_area (10^6 km^2)")

# save
dat_outcrop_full %>% 
  write_rds(here("data", 
                 "outcrop_full.rds"))

# visualize
plot_outcrop_full <- dat_outcrop_full %>%
  ggplot(aes(age, area)) +
  geom_line(linewidth = 1.5, colour = "#a6d0c8") +
  scale_x_reverse() +
  coord_geo(xlim = c(160, 0), 
            ylim = c(0, 6),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Marine outcrop area [10^6 km^2]", 
       title = "Wall, Ivany, Wilkinson 2009") 



# shelf area --------------------------------------------------------------

### Phanerozoic flooded continental area as proportion of earths surface area from Kocsis and Scotese 2021  ###
# load data
dat_cont_area_full <- read_csv(here("data",
                                    "raw",
                                    "shelf_area",
                                    "kocsis_scotese_2021.csv")) %>% 
  # clean up colnames
  rename(area = "shelf-rgeos")

# save
dat_cont_area_full %>% 
  write_rds(here("data", 
                 "proxy_data",
                 "cont_area_full.rds"))

# visualize
plot_cont_area_full <- dat_cont_area_full %>%
  ggplot(aes(age, area)) +
  geom_line(linewidth = 1.5, colour = "#343c24") +
  scale_x_reverse() +
  coord_geo(xlim = c(160, 0), 
            ylim = c(0, 0.25),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Proportion of Earth`s surface", 
       title = "Kocsis and Scotese 2021", 
       subtitle = "Flooded continental area")



### Cenozoic flooded continental area (10^6 km^2) from Miller et al 2005  ###
# load data
dat_cont_area_ceno <- read_csv(here("data",
                                    "raw",
                                    "shelf_area",
                                    "miller_et_al_2005.csv"),
                               col_names = FALSE) %>% 
  # clean up colnames
  rename(age = X2, 
         area = X1)

# save
dat_cont_area_ceno%>% 
  write_rds(here("data", 
                 "proxy_data",
                 "cont_area_ceno.rds"))

# visualize
plot_cont_area_ceno <- dat_cont_area_ceno %>%
  ggplot(aes(age, area)) +
  geom_line(linewidth = 1.5, colour = "#343c24") +
  scale_x_reverse() +
  coord_geo(xlim = c(100, 0), 
            ylim = c(0, 30),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Flooded continental area [10^6 km^2]", 
       title = "Miller et al 2005")



# temperature -------------------------------------------------------------


### Cenozoic deep-ocean temperature  used with equation 7afrom Cramer et al 2011 ###
# load data
dat_temp_ceno <- read_xlsx(here("data",
                                  "raw",
                                  "temperature",
                                  "cramer_et_al_2011_7a.xlsx")) %>% 
  # clean up colnames
  select(age = Age,
         temp = Temperature,
         temp_low = "Temperature min",
         temp_high = "Temperature max", 
         # smoothed long term
         temp_long = "Temperature (long)",
         temp_long_low = "Temperature min (long)",
         temp_long_high = "Temperature max (long)")

# save
dat_temp_ceno %>% 
  write_rds(here("data", 
                 "proxy_data",
                 "temp_ceno.rds"))

# visualize
plot_temp_ceno <- dat_temp_ceno %>%
  ggplot(aes(age, temp_long, 
             ymin = temp_long_low, 
             ymax = temp_long_high)) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             colour = "grey30") +
  geom_ribbon(fill = "#a65852", 
              alpha = 0.3) +
  geom_line(aes(y = temp), 
            linewidth = 1.5, colour = "grey30") +
  geom_line(linewidth = 1.5, colour = "#a65852") +
  scale_x_reverse() +
  coord_geo(xlim = c(70, 0), 
            # ylim = c(0, 0.25),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Deep-ocean temperature [°C]", 
       title = "Cramer et al 2011", 
       subtitle = "First equation\nSmoothed long-term in red")



### Phanerozoic deep-ocean and global average temperature from Scotese et al 2021 ###
# load data
dat_temp_full <- read_xlsx(here("data",
                                  "raw",
                                  "temperature",
                                  "scotese_et_al_2021.xlsx")) %>% 
  # clean up colnames
  select(age = Age,
         temp_gat = GAT,
         temp_deep = "Deep Ocean")

# save
dat_temp_full %>% 
  write_rds(here("data", 
                 "proxy_data",
                 "temp_full.rds"))

# visualize
plot_temp_full <- dat_temp_full %>%
  ggplot(aes(age, temp_gat)) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             colour = "grey30") +
  geom_line(aes(y = temp_deep), 
            linewidth = 1.5, colour = "coral") +
  geom_line(linewidth = 1.5, colour = "#a65852") +
  scale_x_reverse() +
  coord_geo(xlim = c(160, 0), 
            ylim = c(-5, 35),
            height = unit(1.2, "line"), 
            size = 4,
            alpha = 0.4) +
  labs(x = "Age [myr]", 
       y = "Temperature [°C]", 
       title = "Scotese et al 2021", 
       subtitle = "Global average in red\nDeep-ocean in orange")



# save plots --------------------------------------------------------------

# get all plots in a list
plt_list <- mget(ls()[str_detect(ls(), "plot")])

# get their names
plt_names <- ls()[str_detect(ls(), "plot")] %>% 
  str_remove("plot_") %>% 
  str_c(".png")

# save via walk
walk2(.x = plt_list, 
      .y = plt_names, 
      .f = ~ ggsave(.x, filename = here("figures",
                                        "data_plots",
                                        .y), 
                    width = image_width, height = image_height, units = image_units, 
                    bg = "white", device = ragg::agg_png))
