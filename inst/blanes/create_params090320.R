devtools::install_github("gustavdelius/mizer", ref = "12c6e94426cd1225b99434c370eea90ede90f98e")
library(mizer)
library(tidyverse)
params <- readRDS("inst/blanes/params090220.rds")
catch <- readRDS("inst/blanes/catch.rds")

no_sp <- nrow(params@species_params)

## Predator-Prey interaction:  

# - All the fish (starting on hake down the list) should have 1 interaction with 
#   plankton
hake_idx <- which(params@species_params$species == "Hake")
params@species_params$interaction_p[hake_idx:no_sp] <- 1
params <- steady(params)

# - The benthic fish spotted flounder, black goby, gurnard, mullet, 
#   stripped red mullet should have 0.5 of interaction with plankton
sel <- c("Spotted flounder", "Black goby", "Gurnards", "Red mullet", "Striped red mullet")
id <- which(params@species_params$species %in% sel)
params@species_params$interaction_p[id] <- 0.5
params <- steady(params)

# - All the invertebrates (large crustacea, suprabenthic crustacea) should 
#   have 0 interaction
params@species_params$interaction_p[1:7] <- 0
params@species_params$interaction_p[[9]] <- 0
params <- steady(params)

# - (starfish, murex, angular crab and harbour crab a minimum interaction 0.05)
sel <- c("Starfish", "Murex", "Angular crab", "Harbour crab")
id <- which(params@species_params$species %in% sel)
params@species_params$interaction_p[id] <- 0.05
params <- steady(params)

# - Red Snapping shrimp should have 0 interaction with the fish 
#   (starting on spotted flounder down the list)
snap_idx <- which(params@species_params$species == "Red snapping shrimp")
spot_idx <- which(params@species_params$species == "Spotted flounder")
params@interaction[snap_idx, spot_idx:no_sp] <- 0
params <- steady(params)

# - Harbour crab, angular crab, murex and starfish, interaction with fish 
#   should be minimal, 0.05 starting on spotted flounder down the list)
sel <- c("Harbour crab", "Angular crab", "Murex", "Starfish")
id <- which(params@species_params$species %in% sel)
params@interaction[id, spot_idx:no_sp] <- 0.05
params <- steady(params)

## Hake maturity size of 25cm
hake_idx <- which(params@species_params$species == "Hake")
l_mat <- 25
w_mat <- params@species_params$a[[hake_idx]] * l_mat ^ params@species_params$b[[hake_idx]]
params@species_params$w_mat[[hake_idx]] <- w_mat
params@species_params$w_mat25[[hake_idx]] <- NA
params <- setReproduction(params) %>% steady()

## Shortfin squid
squid_idx <- which(params@species_params$species == "Shortfin squid")
l_mat <- 15
w_mat <- params@species_params$a[[squid_idx]] * l_mat ^ params@species_params$b[[squid_idx]]
l_inf <- 25
w_inf <- params@species_params$a[[squid_idx]] * l_inf ^ params@species_params$b[[squid_idx]]
params@species_params$w_mat[[squid_idx]] <- w_mat
params@species_params$w_mat25[[squid_idx]] <- NA
params@species_params$w_inf[[squid_idx]] <- w_inf
params <- setReproduction(params)


params <- tuneParams(params, catch)
# In the app I tuned a lot of the abundances and growth rates but only very
# roughly

saveRDS(params, file="inst/blanes/params090320.rds", version = 2)
