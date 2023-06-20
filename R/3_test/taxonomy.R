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
  # abbreviate order for nicer plotting
  mutate(order = str_replace_all(order, "formes", ".")) %>%
  ggplot(aes(.epred, fct_reorder(order, .epred))) +
  geom_linerange(aes(xmin = .lower, 
                     xmax = .upper), 
                 colour = "grey30") +
  geom_point(aes(size = n, 
                 fill = superorder), 
             shape = 21, 
             colour = "grey20", 
             alpha = 0.8) +
  geom_text(aes(x = .epred, 
                y = fct_reorder(order, .epred), 
                label = order), 
            position = position_nudge(x = c(rep(0.23, 12), 
                                            -0.25)),
            size = 8/.pt, 
            colour = "grey40") +
  scale_fill_manual(name = NULL, 
                    values = c("#EA8778", "#FFBE62"), 
                    limits = c("Sharks", "Rays")) +
  scale_size(name = NULL, 
             breaks = c(100, 250, 500)) +
  scale_y_discrete(breaks = dat_pred_fossil_smr %>%
                     mutate(order = str_replace_all(order, 
                                                    "formes",
                                                    ".")) %>% 
                     arrange(.epred) %>% 
                     pull(order) %>% 
                     .[seq(1, 13, by = 2)]) +
  coord_cartesian(xlim = c(0, 0.6)) +
  labs(y = NULL, 
       x = NULL) +
  theme(legend.position = "none", 
        legend.text = element_text(colour = "grey50", size = 10), 
        axis.text = element_blank(), 
        axis.ticks = element_blank()) 




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
  # abbreviate order for nicer plotting
  mutate(order = str_replace_all(order, "formes", ".")) %>%
  ggplot(aes(fct_reorder(order, .epred), 
             .epred)) +
  geom_linerange(aes(ymin = .lower, 
                     ymax = .upper), 
                 colour = "grey30") +
  geom_point(aes(size = n, 
                 fill = superorder), 
             shape = 21, 
             colour = "grey20", 
             alpha = 0.8) +
  geom_text(aes(y = .epred, 
                x = fct_reorder(order, .epred), 
                label = order), 
            position = position_nudge(x = -0.4, 
                                      y = c(rep(-0.28, 6), 
                                            0.28, 
                                            rep(-0.28, 6))),
            size = 8/.pt, 
            angle = 90, 
            colour = "grey40") +
  scale_x_discrete(breaks = dat_pred_modern_smr %>%
                     mutate(order = str_replace_all(order, 
                                                    "formes",
                                                    ".")) %>% 
                     arrange(.epred) %>% 
                     pull(order) %>% 
                     .[seq(1, 13, by = 2)]) +
  scale_fill_manual(name = NULL, 
                    values = c("#EA8778", "#FFBE62"), 
                    limits = c("Sharks", "Rays")) +
  scale_size(name = NULL, 
             breaks = c(100, 250, 500)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(y = NULL, 
       x = NULL) +
  theme(legend.position = "none", 
        legend.text = element_text(colour = "grey50", size = 10), 
        axis.text = element_blank(), 
        axis.ticks = element_blank()) 



# scatter plot ------------------------------------------------------------

# calculate correlation
cor_pear <- cor.test(dat_pred_fossil_smr$.epred, 
                      dat_pred_modern_smr$.epred,
                      method = "pearson") 

# create label
dat_cor_pear <- paste0("r = ",
                       round(cor_pear$estimate, 2),
                       " [",
                       round(cor_pear$conf.int[[1]], 1),
                       ", ",
                       round(cor_pear$conf.int[[2]], 1),
                       "]") 

# visualise
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
  annotate("label", 
           x = 0.47, 
           y = 0.01, 
           label.size = 0,
           colour = "grey20", 
           size = 8/.pt,
           label = dat_cor_pear) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6), 
                     labels = c("0", "20", "40", "60"), 
                     name = "Fossil extinction risk [%]", 
                     limits = c(0, 0.6)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100"), 
                     name = "Modern extinction risk [%]", 
                     limits = c(0, 1), 
                     position = "right") +
  scale_fill_manual(name = NULL, 
                    values = c("#EA8778", "#FFBE62"), 
                    limits = c("Sharks", "Rays")) +
  theme(legend.position = "none") 
  


# rank plot ---------------------------------------------------------------

# iterate through samples and calculate ranks with uncertainty

# first for fossil ranks

# extract results for each order
dat_ranks_fossil <- dat_pred_fossil %>% 
  group_by(order) %>%
  summarise(sd = sd(.epred)) %>% 
  left_join(dat_pred_fossil_smr %>% 
              mutate(y_rank = dense_rank(.epred)) %>% 
              select(superorder, order, y_rank)) %>% 
  mutate(y_lwr = qnorm(0.025, y_rank, sd),
         y_lwr = floor(y_lwr),
         y_upr = qnorm(0.975, y_rank, sd), 
         y_upr = ceiling(y_upr)) %>% 
  select(-sd)


# and for modern ranks

# extract results for each order
dat_ranks_modern <- dat_pred_modern %>% 
  group_by(order) %>%
  summarise(sd = sd(.epred)) %>% 
  left_join(dat_pred_modern_smr %>% 
              mutate(x_rank = dense_rank(.epred)) %>% 
              select(superorder, order, x_rank)) %>% 
  mutate(x_lwr = qnorm(0.025, x_rank, sd),
         x_lwr = floor(x_lwr),
         x_upr = qnorm(0.975, x_rank, sd), 
         x_upr = ceiling(x_upr)) %>% 
  select(-sd)

# calculate correlation
cor_kend <- dat_ranks_fossil %>%
  left_join(dat_ranks_modern) %>% 
  { cor.test(.$x_rank, .$y_rank, method = "kendall") }

# get function for tau confidence interval from Fieller et al. (1957)
tau.ci <- function(tau, N, conf.level = 0.95, correct='fieller') {
  if(correct=='none') tau.se <- 1/(N - 3)^0.5
  if(correct=='fieller') tau.se <- (0.437/(N - 4))^0.5
  moe <- qnorm(1 - (1 - conf.level)/2) * tau.se
  zu <- atanh(tau) + moe
  zl <- atanh(tau) - moe
  tanh(c(zl, zu))
}

# apply function
cor_kend_ci <- tau.ci(tau = cor_kend$estimate, N = nrow(dat_ranks_fossil))

# create label
dat_cor_kend <- paste0("tau = ",
                       round(cor_kend$estimate, 2),
                       " [",
                       round(cor_kend_ci[[1]], 1),
                       ", ",
                       round(cor_kend_ci[[2]], 1),
                       "]") 



# visualise
plot_rank <- dat_ranks_fossil %>%
  left_join(dat_ranks_modern) %>% 
  ggplot(aes(x_rank, y_rank)) +
  stat_smooth(method = "lm", 
              se = FALSE, 
              colour = "grey70", 
              fullrange = TRUE, 
              linetype = "dashed") +
  geom_linerange(aes(ymin = y_lwr, 
                     ymax = y_upr), 
                 colour = "grey50") +
  geom_linerange(aes(xmin = x_lwr, 
                     xmax = x_upr), 
                 colour = "grey50") +
  geom_point(aes(fill = superorder), 
             shape = 21, 
             colour = "grey20", 
             size = 2.5) +
  annotate("label", 
           x = 3, 
           y = 13, 
           label.size = 0,
           colour = "grey20", 
           size = 8/.pt,
           label = dat_cor_kend) +
  scale_x_continuous(breaks = seq(1, 13, by = 2), 
                     name = "Modern rank") +
  scale_y_continuous(breaks = seq(1, 13, by = 2), 
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
  plot_scatter +
  plot_annotation(tag_levels = "a") &
  theme(legend.position = 'none') 

# save plot
ggsave(plot_final, filename = here("figures",
                                  "taxonomic_agreement.png"), 
       width = image_width, height = image_height,
       units = image_units, 
       bg = "white", device = ragg::agg_png)  
