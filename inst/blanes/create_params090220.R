params <- readRDS("inst/blanes/params_modifiedSdJ2.rds")

# Add small pelagics
meso_sel <- params@w_full >= 1 & params@w_full <= 100
params@kappa <- 1000
params@cc_pp[meso_sel] <- params@kappa * params@w_full[meso_sel]^-params@lambda
comment(params@cc_pp) <- "Modelling small pelagic fish that have a benthic phase"
params@initial_n_pp[] <- params@cc_pp
params@rr_pp <- params@rr_pp * 100
comment(params@rr_pp) <- "High replenishment rate to keep close to carrying capacity"
params@species_params$interaction_p <- params@interaction[, 7]

params@linecolour <- c(params@linecolour, carrion = "lightgrey", detritus = "orange")

params <- tuneParams(params)
# In the app I tuned a lot of the abundances and growth rates but only very
# roughly

saveRDS(params, file="params090220.rds", version = 2)
