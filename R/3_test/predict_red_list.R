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


# visualisations ----------------------------------------------------------



# # visualise transition rates
# future_sim_output[[3]] %>%
#   py_to_r() %>%
#   as_tibble() %>%
#   mutate(year = year + 2023, 
#          status = factor(status, 
#                          levels = rev(c("LC", "NT", "VU", 
#                                     "EN", "CR", "EX")))) %>% 
#   ggplot(aes(year, status_count, 
#              fill = status, 
#              colour = status)) +
#   geom_col(width = 1) +
#   scale_color_manual(values = colorspace::lighten(c("grey10", "#DF5436", "#E79609",
#                                                     "#E7BA29", "lightgreen", "#4C634C"), 0.4)) +
#   scale_fill_manual(values = colorspace::lighten(c("grey10", "#DF5436", "#E79609",
#                                                    "#E7BA29", "lightgreen", "#4C634C"), 0.4)) +
#   scale_x_continuous(breaks = c(2023, 2100, 2200), 
#                      name = "Year") +
#   theme(legend.position = "top")
  

# IUCN history
# read in dataset
dat_status <- read_delim(here("data",
                              "iucn",
                              "CHONDRICHTHYES_iucn_history.txt"))

# plot
plot_status <- dat_status %>% 
  pivot_longer(-species, 
               names_to = "year", 
               values_to = "status") %>% 
  mutate(status = case_when(
    status == "LR/nt" ~ "NT", 
    status == "LR/lc" ~ "LC", 
    status == "LR/cd" ~ "NT", 
    .default = status)) %>% 
  filter(!status %in% c("K", "I")) %>% 
  drop_na() %>% 
  count(year, status) %>% 
  mutate(status = factor(status, 
                         levels = c("LC", 
                                    "NT", 
                                    "VU", 
                                    "EN", 
                                    "CR", 
                                    "DD"))) %>% 
  ggplot(aes(year, 
             n, 
             fill = status)) +
  geom_col(position = position_dodge()) +
  scale_fill_manual(values = c("#2eb774", "#E7BA29", "#E79609", 
                               "#DF5436", "darkred", "grey70"), 
                    name = "Red list status") +
  labs(y = "Number of species", 
       x = NULL) +
  theme(legend.position = c(0.8, 0.8), 
        legend.key.size = unit(3, "mm"))


# save
ggsave(plot_status, 
       filename = here("figures",
                       "supplement", 
                       "redlist_status_history.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)


# number of extinctions
testi <- read_table(here("data",
                "iucn", 
                "te_all_species.txt"), 
           col_names = FALSE) %>% 
  mutate(species = paste(X1, X2), 
         .before = 1) %>% 
  select(-c(X1, X2)) %>% 
  pivot_longer(-species,
               values_to = "te_time") %>% 
  select(-name) %>% 
  drop_na(te_time) %>% 
  group_by(species) %>% 
  mutate(te_mean = mean(te_time)) %>% 
  ungroup() %>% 
  mutate(species = fct_reorder(species, te_mean)) 
  
p1 <- testi %>%
  ggplot(aes(te_time, species)) +
  tidybayes::stat_pointinterval(colour = colour_mint, 
                                shape = 4, 
                                point_colour = colour_coral, 
                                point_size = 0.8, 
                                alpha = 0.06, 
                                point_alpha = 0.1, 
                                point_interval = "median_qi",) +
  scale_y_discrete(expand = c(0.07, 0), 
                   labels = NULL, breaks = NULL) +
  scale_x_continuous(labels = NULL, breaks = NULL) +
  labs(y = NULL, 
       x = NULL) +
  theme(plot.margin = unit(c(-10, 0, 5, 10), "mm"))

  
p2 <- future_sim_output[[3]] %>%
  py_to_r() %>%
  as_tibble() %>%
  ggplot(aes(year, EX)) +
  geom_line(colour = "white", 
            linewidth = 1.5, 
            alpha = 0.9) +
  geom_line(colour = colour_coral, 
            linewidth = 1, 
            alpha = 0.9) +
  scale_x_continuous(labels = function(b) b + 2023) +
  labs(x = "Year", 
       y = "Number of Extinctions")
  
# patch together
p3 <- p1 + inset_element(p2, -0.07, -0.09, 1, 1)  
  

# save
ggsave(p3, filename = here("figures",
                           "supplement",
                           "simulated_extinctions.png"), 
         width = image_width, height = image_height,
         units = image_units, 
         bg = "white", device = ragg::agg_png)


# red list status
p4 <- dat_status %>% 
  pivot_longer(-species, 
               names_to = "year", 
               values_to = "status") %>% 
  distinct(species, status) %>% 
  mutate(status = case_when(
    status == "LR/nt" ~ "NT", 
    status == "LR/lc" ~ "LC", 
    status == "LR/cd" ~ "NT", 
    .default = status)) %>% 
  filter(!status %in% c("K", "I")) %>% 
  drop_na() %>% 
  count(status) %>% 
  mutate(status = factor(status,
                         levels = c("LC",
                                    "NT",
                                    "VU",
                                    "EN",
                                    "CR",
                                    "DD"))) %>% 
  arrange(status) %>% 
  mutate(fraction = n / sum(n), 
         ymax = cumsum(fraction), 
         ymin = c(0, lag(ymax)[-1]), 
         labelPosition = (ymax + ymin) / 2, 
         label = paste0(status, "\n", 
                        as.character(n))) %>% 
  ggplot(aes(
    ymax = ymax,
    ymin = ymin,
    xmax = 4,
    xmin = 3,
    fill = status
  )) +
  geom_rect() +
  geom_label(x = 3.5, 
             aes(y = labelPosition, label = label), 
             size = 10/.pt) +
  scale_fill_manual(values = c("#2eb774", "#E7BA29", "#E79609", 
                               "#DF5436", "darkred", "grey70")) +
  coord_polar(theta = "y") +
  xlim(c(2, 4)) +
  theme_void() +
  theme(legend.position = "none")

# save
ggsave(p4, filename = here("figures",
                           "supplement",
                           "nr_red_list_status.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)

