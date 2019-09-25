library(tidyverse)
library(googlesheets)

species_params <- "species_params" %>% 
    gs_title() %>% 
    gs_read()

sp <- species_params[species_params$`latin name` == "Lesueurigobius sp", ]

params_old <- readRDS("inst/blanes/params.rds")
no_sp_old <- length(params_old@species_params$species)
params <- addSpecies(params_old, sp)
params@rr_pp <- params_old@rr_pp
params@cc_pp <- params_old@cc_pp
params@mu_b[1:no_sp_old, ] <- params_old@mu_b
params@resource_dynamics <- params_old@resource_dynamics
params@resource_params <- params_old@resource_params
params@initial_B <- params_old@initial_B
for (res in names(params@resource_dynamics)) {
    res_var <- paste0("rho_", res)
    params <- set_species_param_default(params, res_var, 0)
}
params@rho <- array(0,
                    dim = dim(params_old@rho) + c(1, 0, 0),
                    dimnames = list(
                        "sp" = as.character(params@species_params$species),
                        "res" = dimnames(params_old@rho)$res,
                        "w" = dimnames(params_old@rho)$w))
params@rho[1:no_sp_old, , ] <- params_old@rho

params <- tuneParams(params)
saveRDS(params, "inst/blanes/params_added_species.rds")
