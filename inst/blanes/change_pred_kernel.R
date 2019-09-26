params <- readRDS("inst/blanes/params_added_species.rds")

ihake <- which(params@species_params$species == "Hake")
params@species_params$pred_kernel_type[ihake] <- "power_law"
params@species_params$kernel_exp <- NA
params@species_params$kernel_exp[ihake] <- -1.524
params@species_params$kernel_l_l <- NA
params@species_params$kernel_l_l[ihake] <- 1.283
params@species_params$kernel_u_l <- NA
params@species_params$kernel_u_l[ihake] <- 5.64
params@species_params$kernel_l_r <- NA
params@species_params$kernel_l_r[ihake] <- 6.6765
params@species_params$kernel_u_r <- NA
params@species_params$kernel_u_r[ihake] <- 5.054
params <- setPredKernel(params)
paams <- tuneParams(params)
save(params, file = "inst/blanes/params_hake_kernel.rds")
