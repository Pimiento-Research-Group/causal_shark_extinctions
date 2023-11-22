library(tidyverse) 
library(here)
library(flextable)
library(officer)
library(ggdist)

# load data ---------------------------------------------------------------

# family level logit values (temperature dependancy)
dat_order <- read_rds(here(here("data", 
                                 "logits", 
                                 "logit_order.rds")))

# family level logit values (temperature dependancy)
dat_family <- read_rds(here(here("data", 
                      "logits", 
                      "logit_family_full.rds")))


# order level table -------------------------------------------------------

# generate nice looking flextable
table_order <- dat_order %>% 
  select(order, logit_val, 
         .lower, .upper) %>% 
  mutate(across(where(is.numeric), ~round(.x, 1)), 
         conf_int = paste0("[", 
                           .lower, 
                           ", ", 
                           .upper, 
                           "]")) %>% 
  select(-c(.lower, .upper)) %>% 
  flextable() %>% 
  theme_vanilla() %>% 
  # set column names
  compose(j = 1, part = "header", 
          value = as_paragraph("Order")) %>%
  compose(j = 2, part = "header", 
          value = as_paragraph("Log-odds")) %>% 
  compose(j = 3, part = "header", 
          value = as_paragraph("95% CI"))


# family level table -------------------------------------------------------

# generate nice looking flextable
table_family <- dat_family %>%
  group_by(family) %>% 
  median_qi(value) %>% 
  select(family, value, 
         .lower, .upper) %>% 
  mutate(across(where(is.numeric), ~round(.x, 1)), 
         conf_int = paste0("[", 
                           .lower, 
                           ", ", 
                           .upper, 
                           "]")) %>% 
  select(-c(.lower, .upper)) %>% 
  flextable() %>% 
  theme_vanilla() %>% 
  # set column names
  compose(j = 1, part = "header", 
          value = as_paragraph("Family")) %>%
  compose(j = 2, part = "header", 
          value = as_paragraph("Log-odds")) %>% 
  compose(j = 3, part = "header", 
          value = as_paragraph("95% CI"))


# percentage change -------------------------------------------------------

# calculate percentage increase of selected families compared to endotherms
dat_perc <- dat_family %>%
  group_by(scale) %>% 
  nest() %>% 
  mutate(outcome = map(data, ~ .x %>% 
                         group_by(metab = if_else(!family %in% c("Lamnidae", "Otodontidae", "Alopiidae"),
                                                  "meso", family)) %>% 
                         nest() %>% 
                         ungroup() %>% 
                         mutate(value_sample = map(data, 
                                                   ~ slice_sample(.x, n =  10000, replace = TRUE) %>% 
                                                     select(value))) %>% 
                         select(-data) %>% 
                         pivot_wider(names_from = metab, values_from = value_sample) %>% 
                         mutate(oto = map2(Otodontidae, meso, 
                                           ~ (.x - .y)/.y),
                                oto = map(oto, median_qi),
                                lamni = map2(Lamnidae, meso, 
                                             ~ (.x - .y)/.y), 
                                lamni = map(lamni, median_qi), 
                                alo = map2(Alopiidae, meso, 
                                           ~ (.x - .y)/.y), 
                                alo = map(alo, median_qi))  %>% 
                         select(oto, lamni, alo) %>% 
                         unnest(cols = c(oto, lamni, alo), 
                                names_sep = "_") %>% 
                         select(where(is.double)) %>% 
                         pivot_longer(cols = everything(),
                                      names_to = c("family", "pointval"),
                                      names_sep = "_") %>% 
                         pivot_wider(names_from = pointval, 
                                     values_from = value))) %>% 
  select(scale, outcome) %>% 
  unnest(outcome) %>% 
  mutate(across(where(is.numeric), ~round(.x, 1)*100), 
         conf_int = paste0("[", 
                           .lower, 
                           ", ", 
                           .upper, 
                           "]")) %>% 
  select(-c(.lower, .upper, .width)) %>% 
  ungroup() %>% 
  left_join(tibble(family = c("oto",
                              "lamni",
                              "alo"), 
                   family_cl = c("Otodontidae", 
                                 "Lamnidae", 
                                 "Alopiidae"))) %>%
  left_join(tibble(scale = c("genus", 
                             "stages", 
                             "ceno"), 
                   scale_cl = c("Genera - Stages",
                                "Species - Stages", 
                                "Species - Cenozoic subset"))) %>% 
  select(family_cl, value, conf_int, scale_cl) 

# make it a nice table
table_perc <- dat_perc %>% 
  flextable() %>% 
  theme_vanilla() %>% 
  # set column names
  compose(j = 1, part = "header", 
          value = as_paragraph("Family")) %>%
  compose(j = 2, part = "header", 
          value = as_paragraph("Percentage Change")) %>% 
  compose(j = 3, part = "header", 
          value = as_paragraph("95% CI")) %>% 
  compose(j = 4, part = "header", 
          value = as_paragraph("Scale")) %>% 
  merge_v(j = 4) %>% 
  fix_border_issues()

# create word document ----------------------------------------------------

# open docx-file and add flextable
my_doc <- read_docx() %>% 
  body_add_flextable(table_order) %>% 
  body_add_break() %>% 
  body_add_flextable(table_family, pos = "after") %>% 
  body_add_break() %>% 
  body_add_flextable(table_perc, pos = "after") 

# convert to word file/ add input to empty docx
print(my_doc, target = here("figures",
                            "supplement", 
                            "supplemental_tables.docx"))

