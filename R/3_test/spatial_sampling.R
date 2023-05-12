# load packages
library(chronosphere)
library(sf)
library(tidyverse)
library(here)

# plotting configurations
source(here("R", "config_file.R"))



# load data ---------------------------------------------------------------


# fossil data and environmental proxy data on species level
dat_occ <-  read_rds(here("data",
                          "fossil_occurrences",
                          "database_occurrences_15_Apr_2023.rds"))


# gplates -----------------------------------------------------------------

# download gplates continental plate model

# see what resolution and version we have
datasets(dat = "paleomap")

# get the data
dat_maps <- fetch("paleomap",
                  var = "paleocoastlines")


# assign to nearest map ---------------------------------------------------

# see what resolution the map is
map_res <- layers(dat_maps) %>%
  str_extract(".+?(?=_)") %>%
  str_extract("(\\d)+") %>%
  unique() %>%
  as.numeric()


# bin fossil data
dat_binned <- dat_occ %>%
  # bin fad and lad to stages
  mutate(bin_low = cut(Max_Ma,
                       breaks = map_res-2.5,
                       labels = map_res[1:length(map_res)-1]),
         bin_high = cut(Min_Ma,
                        breaks = map_res-2.5,
                        labels = map_res[1:length(map_res)-1]))  %>%
  filter(bin_low == bin_high) %>%
  select(everything(), bin = bin_low,
         -c(Max_Ma, Min_Ma, bin_high)) %>%
  mutate(bin = as.numeric(as.character(bin)))



# visualise ---------------------------------------------------------------

# set up list for plot
plot_list <- list()


# loop through maps and bins
plot_list <- unique(dat_binned$bin) %>%
  sort() %>%
  map(~ dat_maps[[paste0("X",
                          .x,
                          "Ma_CM_v7")]] %>%
         as("sf") %>%
         rownames_to_column("id") %>%
         ggplot() +
         geom_sf(colour = "white",
                 fill = "grey70") +
         theme_minimal() +
         geom_point(aes(geometry = geometry),
                    colour = "grey20",
                    fill = alpha("#44978C", 0.4),
                    position = position_jitter(width = 3,
                                               height = 3,
                                               seed = 123),
                    shape = 21,
                    size = 3,
                    stroke = 0.3,
                    stat = "sf_coordinates",
                    data = dat_binned %>%
                      drop_na(paleolon, paleolat) %>%
                      filter(bin == .x) %>%
                      st_as_sf(coords = c("paleolon", "paleolat"),
                               crs = st_crs("WGS84"))) +
         labs(y = NULL,
              x = NULL,
              fill = NULL, title = paste(.x, "myr")) +
         theme(legend.position = "bottom") +
         coord_sf(ylim = c(-90, 90),
                  datum = st_crs("WGS84"))

  )



# save plots
for(i in 1:length(plot_list)){
  ggsave(plot = plot_list[[i]],
         file = here("figures",
                     "maps",
                     paste("file",i,".png",sep="")),
         bg = "white")
}

# create gif
list.files(here("figures",
                "maps"),
           pattern = "*png",
           full.names = TRUE) %>%
  as_tibble() %>%
  mutate(ord_id = str_sub(value, -6, -1),
         ord_id = str_extract(ord_id, "(\\d)+"),
         ord_id = as.numeric(ord_id)) %>%
  arrange(ord_id) %>%
  pull(value) %>%
  gifski::gifski(.,
                 gif_file = here("figures",
                                 "animation_sampling.gif"), 
                 width = 800, height = 600, delay = 1)

