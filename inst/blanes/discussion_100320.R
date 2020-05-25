devtools::install_github("gustavdelius/mizer", ref = "12c6e94426cd1225b99434c370eea90ede90f98e")
library(mizer)
library(tidyverse)
params <- readRDS("inst/blanes/params090320.rds")
catch <- readRDS("inst/blanes/catch.rds")

no_sp <- nrow(params@species_params)

spc <- params@species_params %>% 
    select(species, 
           catch_observed, 
           biomass_observed, 
           abundance_observed, 
           catchability) %>% 
    filter(!is.na(catch_observed)) %>% 
    mutate(ratio = catch_observed / biomass_observed)

params <- tuneParams(params, catch)


saveRDS(params, file="inst/blanes/params160320.rds", version = 2)
