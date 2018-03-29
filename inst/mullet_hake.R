 
######### returne abundance test
params <- set_scaling_model()
params@A[] <- NA
params@A[length(params@A)] <- 1
retune_abundance(params)

######### get scaling model

params <- set_scaling_model(max_w_inf = 5*10^3)
 params@species_params$r_max <- params@species_params$w_mat
 params@species_params$r_max[] <- 10^50
 params@A[] <- NA
 
######### add mullet
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

 ############# add hake 
 species_params <- data.frame(
   species = "hake",
   w_min = 0.001,
   w_inf = 3964,
   w_mat = 173.9, # use better value here later
   h = NA, # will compute this later
   #ks = 4, 
   ks = 1/2,
   beta = 11.02,
   sigma = 1.1,
   z0 = 0,
   alpha = 0.4,
   erepro = 0.1,
   sel_func = "knife_edge", # not used but required
   knife_edge_size = 100,
   gear = "knife_edge_gear",
   k = 0,
   gamma = NA,
   w_min_idx = NA,
   r_max = 10^50 #why do I need r_max after combining before
 )
 k_vb <- 0.1
 species_params$h <- 3*k_vb*(species_params$w_inf^(1/3))/(species_params$alpha*params@f0)
 ae <- sqrt(2*pi) * species_params$sigma * species_params$beta^(params@lambda-2) * exp((params@lambda-2)^2 * species_params$sigma^2 / 2)
 species_params$gamma <- (species_params$h / (params@kappa * ae)) * (params@f0 / (1 - params@f0))
 species_params$w_min_idx <- sum(params@w<=species_params$w_min)
 
 params_out_2 <- add_species(params_out, species_params, mult = 5.5 * 10 ^ (8))
 sim <- project(params_out_2, t_max = 5, effort = 0)
 plot(sim)
 
 # #18 #24 #29 Have got code that holds hake and mullet (pushed to inst/mullet_hake.R in adsp branch). 
 # Next I want to retune the abundance multipliers to be more reasonable, and do some experiments with fishing gears.
 