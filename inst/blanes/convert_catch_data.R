library(tidyverse)
library(googlesheets)
params <- readRDS("inst/blanes/params.rds")

size_frequency_sheet <- gs_title("size-frequency")
size_frequency <- gs_read(size_frequency_sheet)

catch <- gs_title("size-frequency") %>% 
    gs_read() %>% 
    pivot_longer(-cm, 
                 names_to = "latin.name", 
                 values_to = "catch",
                 values_drop_na = TRUE) %>% 
    separate(cm, c("start", "end"), 
             sep = "-",
             convert = TRUE) %>% 
    mutate(length = start - 0.1,
           dl = end - length) %>% 
    left_join(params@species_params, by = "latin.name") %>% 
    select(species, latin.name, length, dl, catch) %>% 
    filter(!is.na(species))
catch
saveRDS(catch, "inst/blanes/catch.rds")