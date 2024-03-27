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
  arrange(logit_val) %>% 
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
  arrange(value) %>% 
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



# create word document ----------------------------------------------------

# open docx-file and add flextable
my_doc <- read_docx() %>% 
  body_add_flextable(table_order) %>% 
  body_add_break() %>% 
  body_add_flextable(table_family, pos = "after") 

# convert to word file/ add input to empty docx
print(my_doc, target = here("figures",
                            "supplement", 
                            "supplemental_tables.docx"))

