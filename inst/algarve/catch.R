# Script to bin the individual weight observations in fish_sizes_sp.csv
# into binned catch data in catch.rds

library(tidyverse)
params <- readRDS("inst/algarve/params_36.rds")
fish_sizes_sp <- read_csv("inst/algarve/fish_sizes_sp.csv")

# The species names are in the first row
names(fish_sizes_sp) <- as.character(slice(fish_sizes_sp, 1))
# Hack due to the fact that Dourada was misspelled in species data frame
names(fish_sizes_sp)[1] <- "Dourada "

catch <- fish_sizes_sp %>% 
    slice(-1) %>%  # Get rid of row containing species names
    mutate_all(as.numeric) %>% 
    gather(everything(), key = "species", value = "weight", na.rm = TRUE) %>% 
    mutate(cut = cut(weight, breaks = params@w, right = FALSE)) %>% 
    group_by(species, cut) %>% 
    summarise(catch = n()) %>% 
    separate(cut, into = c(NA, "weight", "right", NA), 
             sep = "[\\[,\\)]", convert = TRUE) %>% 
    mutate(dw = right - weight) %>% 
    select(species, weight, dw, catch)

saveRDS(catch, "inst/algarve/catch.rds")
