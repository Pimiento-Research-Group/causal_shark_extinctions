# load packages
library(tidyverse)
library(here)
library(patchwork)
library(ggdist)
library(deeptime)
library(readxl)


# plotting configurations
source(here("R", "config_file.R"))



# load data  ------------------------------------------------------------

# get predictions from different scales/ data sets
set.seed(123)
dat_pred <- paste0("pred_temperature_", 
                   c("ceno", 
                     "genus", 
                     "stage"), 
                   ".rds") %>% 
  map_df(~ read_rds(here("data",
                         "predictions",
                         .x)) %>% 
           group_by(temperature) %>% 
           slice_sample(n = 100) %>% 
           ungroup() %>% 
           add_column(dataset = .x) %>% 
           mutate(dataset = str_remove_all(dataset, 
                                           "pred_temperature_"), 
                  dataset = str_remove_all(dataset, 
                                           ".rds")))


# predicted logit-values
dat_pred_full <- read_rds(here("data",
                               "logits",
                               "logit_full.rds"))

# load temperature data
dat_temp <- # Phanerozoic deep-ocean and global average temperature from Scotese et al 2021 
  read_xlsx(here("data",
                 "raw",
                 "temperature",
                 "scotese_et_al_2021.xlsx")) %>% 
  # clean up colnames
  select(age = Age,
         temp = "Deep Ocean") %>% 
  #Cenozoic deep-ocean temperature  used with equation 7afrom Cramer et al 2011
  bind_rows(read_xlsx(here("data",
                           "raw",
                           "temperature",
                           "cramer_et_al_2011_7a.xlsx")) %>%
              # clean up colnames
              select(age = Age,
                     temp = Temperature))

# family level logit values (temperature dependancy)
dat_family <- read_rds(here(here("data",
                                 "logits", 
                                 "logit_family_full.rds"))) %>% 
  group_by(family) %>% 
  median_qi(value) %>% 
  # add superorder
  left_join(read_rds(here("data",
                          "fossil_occurrences",
                          "database_occurrences_01_Nov_2023.rds")) %>% 
              drop_na(family, superorder) %>% 
              filter(order != "incertae sedis") %>% 
              count(superorder, family))

# load 95 bin Phanerozoic time scale based on the stratigraphic stages 
# of Gradstein et al. 2020
load(here("data", 
          "stages.Rdata"))


# temperature-extinction relationship -------------------------------------

# visualise trends
plot_temp <- dat_pred %>%
  mutate(dataset = case_when(
    dataset == "ceno" ~ "1myr", 
    dataset == "genus" ~ "Genus", 
    dataset == "stage" ~ "Stages"
  ), 
  dataset = factor(dataset, 
                   levels = c("Stages", 
                              "Genus", 
                              "1myr"))) %>% 
  ggplot(aes(temperature,
             value,
             group = interaction(nr_draw, dataset),
             colour = dataset)) +
  geom_line(alpha = 0.2) +
  scale_colour_manual(
    values = c("#4C634C",
               colour_coral, 
               colour_purple),
    labels = c("Species - Stages", 
               "Genera - Stages", 
               "Species - Cenozoic subset"),
    name = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(alpha = 0.9, 
                                                   linewidth = 0.6))) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
                     labels = c("0", "20", "40", "60", "80", "100")) +
  coord_cartesian(ylim = c(0, 1), 
                  xlim = c(0, 25)) +
  labs(y = "Extinction Risk [%]", 
       x = "Global Temperature [°C]") +
  theme(legend.position = c(0.75, 0.75))


plot_temp

# logit through time ------------------------------------------------------

# create plot
plot_logit <- dat_pred_full %>%
  filter(data_set != "Future") %>% 
  mutate(data_set = factor(data_set, 
                           levels = c("Stages", 
                                      "Genus", 
                                      "1myr", 
                                      "Modern"))) %>% 
  ggplot(aes(myr, value, group = data_set)) +
  annotate("rect", 
           xmin = stages$bottom[c(81, 87)], 
           xmax = stages$top[c(81, 87)], 
           ymin = -Inf, 
           ymax = Inf, 
           fill = "#169199", 
           alpha = 0.2) + 
  annotate("rect", 
           xmin = stages$bottom[c(76, 83, 91)], 
           xmax = stages$top[c(76, 83, 91)], 
           ymin = -Inf, 
           ymax = Inf, 
           fill = "#C75E6B", 
           alpha = 0.2) + 
  annotate("text",
           y = 1.8, 
           x = 112, 
           colour = "#C75E6B", 
           size = 9/.pt, 
           label = "Hyperthermal", 
           alpha = 0.9) +
  annotate("text",
           y = 1.8, 
           x = 83, 
           colour = "#169199", 
           size = 9/.pt, 
           label = "Hypothermal", 
           alpha = 0.9) +
  geom_hline(yintercept = 0, 
             linetype = "dotted", 
             colour = colour_grey) +
  geom_linerange(aes(xmin = xmin,
                     xmax = xmax),
                 colour = "grey75") +
  geom_linerange(aes(ymin = .lower,
                     ymax = .upper),
                 colour = "grey75") +
  geom_point(aes(fill = data_set), 
             size = 2, 
             shape = 21, 
             colour = "white") +
  scale_x_reverse() +
  scale_fill_manual(values = c("#4C634C", 
                               colour_coral, 
                               colour_purple, 
                               "#FF5F1F"), 
                    labels = c("Species - Stages", 
                               "Genera - Stages", 
                               "Species - Cenozoic subset", 
                               "Species - Modern"), 
                    name = NULL) +
  scale_y_continuous(limits = c(-8.9, 2), 
                     breaks = seq(-8, 2, by = 2)) +
  guides(fill = guide_legend(override.aes = list(size = 2.5))) +
  coord_geo(xlim = c(150, 0), 
            dat = list("epochs", "periods"),
            pos = list("b", "b"),
            alpha = 0.2, 
            height = unit(0.8, "line"), 
            size = list(5/.pt, 9/.pt),
            lab_color = "grey20", 
            color = "grey50", 
            abbrv = list(TRUE, FALSE), 
            fill = "white",
            expand = TRUE, 
            lwd = list(0.1, 0.2)) +
  labs(x = "Million years", 
       y = "Temperature dependancy\n[log-odds]") +
  theme(legend.position = "none") 


plot_logit



# logit per family --------------------------------------------------------


  
plot_fam <- dat_family %>%
  mutate(endo = if_else(family %in% c("Alopiidae",
                                      "Lamnidae",
                                      "Otodontidae"), 
                        "Endothermic", 
                        "Ectothermic")) %>% 
  ggplot(aes(y = family, 
             x = value)) +
  geom_vline(xintercept = 0, 
             linetype = "dotted", 
             colour = colour_grey) +
  geom_linerange(aes(xmin = .lower, 
                     xmax = .upper), 
                 alpha = 0.3, 
                 linewidth = 0.1) +
  geom_point(aes(fill = endo), 
             size = 2.5, 
             shape = 21, 
             colour = "white") +
  scale_fill_manual(values = c("grey60", 
                               "#FFBE62"), 
                    name = NULL, 
                    labels = rev) +
  scale_x_continuous("Temperature dependancy [log-odds]", 
                     expand = expansion(mult = c(0.1, 0)), 
                     breaks = c(-4, -2, 0)) +
  scale_y_discrete("Families", 
                   expand = expansion(mult = c(0.2, 0.3))) +
  guides(size = "none", 
         fill = guide_legend(nrow = 1, 
                             reverse = TRUE)) +
  theme(legend.position = c(0.5, 0.95), 
        axis.ticks.y = element_blank(), 
        axis.text.y = element_blank())

plot_fam

# patch together and save -------------------------------------------------

plot_full <- plot_temp /
  plot_logit /
  plot_fam +
  plot_annotation(tag_levels = "a")


# save plot
ggsave(plot_full, filename = here("figures",
                                  "figure_2.png"),
       width = image_width, height = image_height*2,
       units = image_units,
       bg = "white", device = ragg::agg_png)


