##set_trait_model <- function(no_sp = 10,
set_scaling_model <- function(no_sp = 11,
                            min_w_inf = 10,
                            ##max_w_inf = 1e5,
                            max_w_inf = 100,
                            ## no_w = 100,
                            no_w = 701,
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
                            eta_egg = 1e-05, ## new for this model type
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
                            erepro = 0.1,
                            ...){
  #' ### Set parameters 
  #' Global Parameters
  
  p <- n
  lambda <- 2+q-n
  
  # ----
  #' ### Set grid points and characteristic sizes 
  
  
  # max_w already set
  min_w <- eta_egg*min_w_inf
  # we need to do coercion of such quantities to gridpoints
  
  
  # #  # #
  
  dist_sp <- 0.2
  minimum_egg <- -4
  maximum_egg <- -1.9
  no_sp <- (maximum_egg-minimum_egg)/dist_sp + 1
  species <- 1:no_sp
  x_min <- seq(minimum_egg, by = dist_sp, length.out = no_sp)
  
  w_min <- 10^x_min
  w_inf <- 10^(x_min+5)
  w_mat <- 10^(x_min+4.4)  # This is about a quarter of w_inf
  min_w <- min(w_min)
  max_w <- max(w_inf)
  no_w <- log10(max_w/min_w)*100+1
  min_w_pp <- 1e-8
}