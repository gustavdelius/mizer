##set_trait_model <- function(no_sp = 10,
set_scaling_model <- function(no_sp = 11,
                            min_w_inf = 10,
                            ##max_w_inf = 1e5,
                            max_w_inf = 1000,
                            ## no_w = 100,
                            ## no_w = 701,
                            min_w = 0.001,
                            ## max_w = max_w_inf * 1.1,
                            max_w = max_w_inf,
                            ##min_w_pp = 1e-10,
                            min_w_pp = 1e-8,
                            no_w_pp = NA,
                            ##w_pp_cutoff = 1,
                            ## dont know where above is used
                            ##k0 = 50, # recruitment adjustment parameter
                            ## dont want to use above
                            n = 2/3,
                            ##p = 0.75,
                            # p=n
                            ##q = 0.9,
                            q = 3/4,
                            ## eta = 0.25,
                            eta = 0.2511886,
                            min_w_egg = 0.0001, ## new for this model type
                            ## r_pp = 4,
                            r_pp = 1e-1,
                            ##kappa = 0.005,
                            kappa = 7e10,
                            ##lambda = 2+q-n, ## we will force this later
                            ##alpha = 0.6,
                            alpha = 0.4,
                            ks = 4,
                           ## z0pre = 0.6, we dont wanr constant death term
                            h = 30,
                            beta = 100,
                            sigma = 1.3,
                           ## f0 = 0.5,
                            f0 = 0.6,
                            gamma = NA,
                           ## knife_edge_size = 1000,
                            knife_edge_size = 100,
                            gear_names = "knife_edge_gear",
                            chi = 0,
                            ##erepro = 0.1,
                           no_w = log10(max_w/min_w_egg)*100+1,
                            ...){
  
  #' ### Set parameters 
  #' Global Parameters
  erepro <- 0.1
  p <- n
  lambda <- 2+q-n
  
  # ----
  #' ### Set grid points and characteristic sizes 
  # we need to do coercion of such quantities to gridpoints
  
  
  # #  # #
  # max_w already set
  min_w <- min_w_egg
  minimum_egg <- log10(min_w_egg)
  maximum_egg <- log10(max_w_inf*min_w_egg/min_w_inf)
  dist_sp <- (maximum_egg-minimum_egg)/(no_sp-1)
  ## this goes wrong when no_sp=1, do we need if statement for this case ? 
  species <- 1:no_sp
  x_min <- seq(minimum_egg, by = dist_sp, length.out = no_sp)
  w_min <- 10^x_min
  w_inf <- min_w_inf*w_min/min_w_egg
  w_mat <- eta*w_inf
  
  # ----
  #' ### Build Params Object 
  
  species_params <- data.frame(
    species = 1:no_sp,
    w_min = w_min,
    w_inf = w_inf,
    w_mat = w_mat,
    h = h,
    ks = ks,
    beta = beta,
    sigma = sigma,
    z0 = 0,
    alpha = alpha,
    erepro = erepro,
    sel_func = "knife_edge", # not used but required
    knife_edge_size = knife_edge_size
  )
  
  params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
                        kappa = kappa, min_w = min_w, max_w = max_w, no_w = no_w, 
                        min_w_pp = min_w_pp, w_pp_cutoff = max_w, r_pp = r_pp,
                        chi = chi)
  
  gamma <- params@species_params$gamma[1]
  w <- params@w
  
  
  # ----
  
  # we need to calculate the steady state next, and we determine the solution, and output 
  # it as a new slot, and we need to determine the values of the erepro, cc_pp and mu_b, 
  # we have to determine the form of the analytic solution. So we should output that too, 
  # as part of the output of set_scaling_model.R
  
}