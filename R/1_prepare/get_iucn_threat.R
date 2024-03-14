library(tidyverse)
library(jsonlite)
library(httr)
library(here)

# get species names and ids in the iucn red list --------------------------

# link to the API output as a JSON file
token <- "1fecbeea639ba430f60510af483c5d4b282e3f097aa98e4613003c9903970df6"
url_json <- paste0("https://apiv3.iucnredlist.org/api/v3/comp-group/getspecies/sharks_and_rays?token=", 
                   token)

# get the raw json into R
raw_json <- GET(url_json) %>%
  content()

# create the dataframe and tidy it up
ex_output <- pluck(raw_json, "result") %>%
  enframe() %>%
  unnest_wider(value) 


# get the threat level of each species via id -----------------------------

# create list tibble, this can take quite some time
dat_threat <- ex_output %>% 
  select(taxonid, scientific_name, category) %>% 
  mutate(url_id = paste0("https://apiv3.iucnredlist.org/api/v3/threats/species/id/",
                         taxonid ,
                         "/?token=",
                         token)) %>% 
  mutate(thread_tibble = map(url_id, 
                             ~ GET(.x) %>% 
                               content() %>% 
                               pluck("result") %>% 
                               enframe() %>% 
                               unnest_wider(value), 
                             .progress = TRUE))

# clean up 
dat_threat %>% 
  unnest(thread_tibble) %>% 
  # save
  write_rds(here("data", 
                 "iucn", 
                 "threat_per_spp.rds"))

