# load packages
library(tidyverse)
library(here)
library(dagitty)
library(brms)
library(tidybayes)
library(patchwork)

# plotting configurations
source(here("R", "config_file.R"))



# load data  ------------------------------------------------------------

# fossil data and environmental proxy data on species level
dat_merged <- read_rds(here("data",
                            "processed_merged_data.rds")) %>% 
  drop_na(order) %>% 
  mutate(superorder = fct_collapse(order,
                              Sharks = c("Carcharhiniformes",
                                         "Orectolobiformes",
                                         "Lamniformes",
                                         "Squaliformes",
                                         "Hexanchiformes",
                                         "Echinorhiniformes",
                                         "Heterodontiformes",
                                         "Pristiophoriformes",
                                         "Squatiniformes",
                                         "Synechodontiformes"),
                              Rays = c("Myliobatiformes",
                                       "Rajiformes",
                                       "Rhinopristiformes",
                                       "Torpediniformes"))) 


# read iucn data
dat_iucn <- read_delim(here("data",
                            "iucn",
                            "CHONDRICHTHYES_iucn_history.txt"))

# reformat
dat_iucn <- dat_iucn %>% 
  pivot_longer(cols = -c(species), 
               names_to = "year", 
               values_to = "iucn_status") %>% 
  drop_na(iucn_status) %>% 
  mutate(ext_signal = if_else(iucn_status %in% c("LR/lc", "LC",
                                                 "LR / cd", "NT" ,
                                                 "LR / nt"),
                              0, 1))

# read temperature data
dat_ipcc <- tibble(filename = list.files(here("data",
                                              "ipcc")) %>%
                     str_remove("tas_global_") %>%
                     str_remove(".csv")) %>% 
  mutate(file_contents = map(list.files(here("data",
                                             "ipcc"),
                                        full.names = TRUE),  
                             ~ read_csv(.x))) %>% 
  filter(filename %in% c("Historical", "SSP1_1_9")) %>% 
  unnest(file_contents) %>% 
  select(year = Year, temp = Mean) %>% 
  filter(year %in% unique(dat_iucn$year))

# reformat
dat_ipcc <- dat_ipcc %>% 
  mutate(temp_st = temp - lead(temp), 
         temp_lt1 = lead(temp_st), 
         temp_lt2 = lead(temp_st, n = 2), 
         temp_lt3 = lead(temp_st, n = 3), 
         temp_lt4 = lead(temp_st, n = 4)) %>% 
  fill(contains("temp"), .direction = "downup") %>% 
  mutate(year = as.character(year))


# merge together
dat_modern <- left_join(dat_iucn, dat_ipcc)

# add order to modern data
dat_modern <- dat_modern %>% 
  mutate(genus = word(species, 1)) %>% 
  left_join(dat_merged %>% 
              distinct(superorder, order, genus)) %>% 
  drop_na(superorder, order, genus)


# total effect adjustment sets --------------------------------------------


# load the graph 
dag <- downloadGraph("dagitty.net/mjiV5Qf")


# get adjustments sets for the total effect
adjustmentSets(dag, 
               exposure = "taxonomic identity",
               outcome = "extinction risk", 
               effect = "total")

# no need to account for confounders


# fossil data --------------------------------------------------------------


# fit final model
mod1 <- brm_logistic("ext_signal ~ order")

# extract predictions 

# set up grid to average over
dat_pred_fossil <- dat_merged %>% 
  distinct(order, superorder) %>% 
  # select only those that have modern equivalents
  filter(order %in% dat_modern$order) %>% 
  # add posterior draws
  add_epred_draws(mod1, 
                  ndraws = 1000) %>% 
  ungroup()
  
# summarise
dat_pred_fossil_smr <- dat_pred_fossil %>% 
  group_by(superorder, order) %>% 
  median_qi(.epred) %>% 
  # add number of samples
  full_join(count(dat_merged, order)) %>% 
  # clean up
  drop_na(superorder) %>% 
  select(-c(.width, .point, .interval))

# visualise
plot_fossil <- dat_pred_fossil_smr %>%
  ggplot(aes(.epred, fct_reorder(order, .epred))) +
  geom_linerange(aes(xmin = .lower, 
                     xmax = .upper), 
                 colour = "grey30") +
  geom_point(aes(size = n, 
                 fill = superorder), 
             shape = 21, 
             colour = "grey20", 
             alpha = 0.8) +
  scale_y_discrete(position = "right") +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  scale_fill_manual(name = NULL, 
                    values = c("#EA8778", "#FFBE62"), 
                    limits = c("Sharks", "Rays")) +
  scale_size(name = NULL, 
             breaks = c(100, 250, 500)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(y = NULL, 
       x = "Extinction Risk [%]") +
  theme(legend.position = c(0.8, 0.3), 
        legend.text = element_text(colour = "grey50", size = 10)) 




# modern data ------------------------------------------------------

# repeat for modern data
dat_merged <- dat_modern

# fit final model
mod2 <- brm_logistic("ext_signal ~ order")

# extract predictions 

# set up grid to average over
dat_pred_modern <- dat_merged %>% 
  distinct(order, superorder) %>% 
  # add posterior draws
  add_epred_draws(mod2, 
                  ndraws = 1000) %>% 
  ungroup()

# summarise
dat_pred_modern_smr <- dat_pred_modern %>% 
  group_by(superorder, order) %>% 
  median_qi(.epred) %>% 
  # add number of samples
  full_join(count(dat_merged, order)) %>% 
  # clean up
  drop_na(superorder) %>% 
  select(-c(.width, .point, .interval))


# visualise
plot_modern <- dat_pred_modern_smr %>%
  ggplot(aes(.epred, fct_reorder(order, .epred))) +
  geom_linerange(aes(xmin = .lower, 
                     xmax = .upper), 
                 colour = "grey30") +
  geom_point(aes(size = n, 
                 fill = superorder), 
             shape = 21, 
             colour = "grey20", 
             alpha = 0.8) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  scale_fill_manual(name = NULL, 
                    values = c("#EA8778", "#FFBE62"), 
                    limits = c("Sharks", "Rays")) +
  scale_size(name = NULL, 
             breaks = c(100, 250, 500)) +
  coord_cartesian(xlim = c(0, 1)) +
  labs(y = NULL, 
       x = "Extinction Risk [%]") +
  theme(legend.position = "none", 
        legend.text = element_text(colour = "grey50", size = 10)) 



# scatter plot ------------------------------------------------------------

plot_scatter <- dat_pred_fossil_smr %>%
  left_join(dat_pred_modern_smr %>% 
              rename(.epred_mod = .epred, 
                     .lower_mod = .lower, 
                     .upper_mod = .upper, 
                     n_mod = n)) %>% 
  ggplot(aes(.epred, .epred_mod)) +
  stat_smooth(method = "lm", 
              se = FALSE, 
              colour = "grey70", 
              fullrange = TRUE, 
              linetype = "dashed") +
  geom_linerange(aes(xmin = .lower,
                     xmax = .upper), 
                 colour = "grey50") +
  geom_linerange(aes(ymin = .lower_mod,
                     ymax = .upper_mod), 
                 colour = "grey50") +
  geom_point(aes(fill = superorder), 
             shape = 21, 
             colour = "grey20", 
             size = 2.5) +
  scale_x_continuous(breaks = c(0.2, 0.4, 0.6), 
                     labels = c("20", "40", "60"), 
                     name = "Fossil extinction risk [%]") +
  scale_y_continuous(breaks = c(0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("20", "40", "60", "80", "100"), 
                     name = "Modern extinction risk [%]", 
                     limits = c(0.05, 1.05), 
                     position = "right") +
  scale_fill_manual(name = NULL, 
                    values = c("#EA8778", "#FFBE62"), 
                    limits = c("Sharks", "Rays")) +
  theme(legend.position = "none") 
  


# rank plot ---------------------------------------------------------------

# loop through samples and calculate ranks

# first for fossil ranks

# first set up matrix to save results
rank_vec <- matrix(nrow = 13, ncol = 1000)

# loop through samples
for (i in 1:1000) {
  rank_vec[, i] <- dat_pred_fossil %>% 
    filter(.draw == i) %>% 
    mutate(r_rank = dense_rank(.epred)) %>% 
    pull(r_rank) 
}

# extract results for each order
dat_ranks_fossil <- 1:13 %>% 
  map_df(~ rank_vec[.x, ] %>% 
           tabulate() %>% 
           desc() %>% 
           dense_rank() %>% 
           { . %in% c(1, 2, 3) } %>% 
           which() %>% 
           enframe() %>% 
           pivot_wider(names_from = name, 
                       values_from = value)) %>% 
  rename(xmin = 1, 
         x = 2, 
         xmax = 3) %>% 
  bind_cols(dat_pred_fossil %>% 
              filter(.draw == 1) %>% 
              select(superorder, order))


# and for modern ranks

# first set up matrix to save results
rank_vec <- matrix(nrow = 13, ncol = 1000)

# loop through samples
for (i in 1:1000) {
  rank_vec[, i] <- dat_pred_modern %>% 
    filter(.draw == i) %>% 
    mutate(r_rank = dense_rank(.epred)) %>% 
    pull(r_rank) 
}

# extract results for each order
dat_ranks_modern <- 1:13 %>%
  map_df(~ rank_vec[.x, ] %>% 
           tabulate() %>% 
           desc() %>% 
           dense_rank() %>% 
           { . %in% c(1, 2, 3) } %>% 
           which() %>% 
           enframe() %>% 
           pivot_wider(names_from = name, 
                       values_from = value)) %>% 
  rename(ymin = 1, 
         y = 2, 
         ymax = 3) %>% 
  bind_cols(dat_pred_modern %>% 
              filter(.draw == 1) %>% 
              select(superorder, order))

# visualise
plot_rank <- dat_ranks_fossil %>% 
  left_join(dat_ranks_modern) %>% 
  ggplot(aes(x, y)) +
  stat_smooth(method = "lm", 
              se = FALSE, 
              colour = "grey70", 
              fullrange = TRUE, 
              linetype = "dashed") +
  geom_linerange(aes(ymin = ymin, 
                     ymax = ymax), 
                 colour = "grey50") +
  geom_linerange(aes(xmin = xmin, 
                     xmax = xmax), 
                 colour = "grey50") +
  geom_point(aes(fill = superorder), 
             shape = 21, 
             colour = "grey20", 
             size = 2.5) +
  scale_x_continuous(expand = c(0.15, 0.15), 
                     breaks = 1:13, 
                     name = "Modern rank") +
  scale_y_continuous(expand = c(0.15, 0.15), 
                     breaks = 1:13, 
                     name = "Fossil rank") +
  scale_fill_manual(name = NULL, 
                    values = c("#EA8778", "#FFBE62"), 
                    limits = c("Sharks", "Rays")) +
  theme(legend.position = "none") 

  

# merge plots -------------------------------------------------------------


# patch plots together
plot_final <- plot_rank + 
  plot_fossil +
  plot_modern + 
  plot_scatter &
  theme(legend.position = 'none')

# save plot
ggsave(plot_final, filename = here("figures",
                                  "taxonomic_agreement.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)  
