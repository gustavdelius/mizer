
devtools::install_github("sizespectrum/mizer")
library(mizerExperimental)
library(tidyverse)

params160320 <- readRDS("inst/blanes/params160320.rds")
params <- upgradeParams(params160320)

# reconstruct resource params
maxidx <- max(which(params@cc_pp > 0))  # The largest resource index
r_pp <- params160320@rr_pp[[maxidx]] / params@w_full[[maxidx]] ^ (params160320@n - 1)
params@resource_params <- list(
    r_pp = r_pp,
    lambda = params160320@lambda,
    kappa = params160320@kappa,
    n = params160320@n,
    w_pp_cutoff = max(params@w_full[params@cc_pp > 0]))

params <- params %>% 
    setComponent("detritus",
                 initial_value = params160320@initial_B[["detritus"]],
                 dynamics_fun = "constant_dynamics",
                 encounter_fun = "encounter_contribution",
                 component_params = 
                     list(rho = outer(params@species_params$rho_detritus, 
                                      params@w^params@resource_params$n))
                 ) %>% 
    setComponent("carrion",
                 initial_value = params160320@initial_B[["carrion"]],
                 dynamics_fun = "constant_dynamics",
                 encounter_fun = "encounter_contribution",
                 component_params = 
                     list(rho = outer(params@species_params$rho_carrion, 
                                      params@w^params@resource_params$n))
                 )


saveRDS(params, file = "inst/blanes/params250520.rds", version = 2)
