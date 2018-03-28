 
######### returne abundance test
params <- set_scaling_model()
params@A[] <- NA
params@A[length(params@A)] <- 1
retune_abundance(params)

######### add_species test

params <- set_scaling_model()
 params@species_params$r_max <- params@species_params$w_mat
 params@species_params$r_max[] <- 10^50
 params@A[] <- NA
 species_params <- data.frame(
   species = "mullet",
   w_min = 0.001,
   w_inf = 251.94,
   w_mat = 16.48,
   h = NA, # will compute this later
   ks = 4,
   beta = 283,
   sigma = 1.8,
   z0 = 0,
   alpha = 0.4,
   erepro = 0.1,
   sel_func = "knife_edge", # not used but required
   knife_edge_size = 100,
   gear = "knife_edge_gear",
   k = 0,
   gamma = NA,
   w_min_idx = NA,
   r_max = 10^50
 )
 k_vb <- 0.6
 species_params$h <- 3*k_vb*(species_params$w_inf^(1/3))/(species_params$alpha*params@f0)
 ae <- sqrt(2*pi) * species_params$sigma * species_params$beta^(params@lambda-2) * exp((params@lambda-2)^2 * species_params$sigma^2 / 2)
 species_params$gamma <- (species_params$h / (params@kappa * ae)) * (params@f0 / (1 - params@f0))
 species_params$w_min_idx <- sum(params@w<=species_params$w_min)
 params_out <- add_species(params, species_params, mult = 5.5 * 10 ^ (8))
 sim <- project(params_out, t_max = 5, effort = 0)
 plot(sim)
 species_params$species <- "mullet2"
 species_params$w_inf <- 50
 params_out_2 <- add_species(params_out, species_params, mult = 5.5 * 10 ^ (8))
 sim <- project(params_out_2, t_max = 5, effort = 0)
 plot(sim)