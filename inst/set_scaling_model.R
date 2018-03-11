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
  
############### next part just copied from stable_community.R, needs checking #######
  
  #' gamma is determined by mizerparams. Note that density dependence is currently off
  gamma <- params@species_params$gamma[1]
  w <- params@w
  
  
  # ----
  #' ### Determine analytic solution
  
  
  mu0 <- (1-f0) * sqrt(2*pi) * kappa * gamma * sigma *
    (beta^(n-1)) * exp(sigma^2 * (n-1)^2 / 2)
  hbar <- alpha * h * f0 - ks
  pow <- mu0/hbar/(1-n)
  n_mult <- (1 - (w/w_inf[1])^(1-n))^(pow-1) * (1 - (w_mat[1]/w_inf[1])^(1-n))^(-pow)
  n_mult[w < w_mat[1]] <- 1
  n_mult[w >= w_inf[1]] <- 0
  n_exact <- w  # Just to get array with correct dimensions and names
  n_exact <- ((w_min[1]/w)^(mu0/hbar) / (hbar * w^n)) * n_mult
  n_exact[w < w_min[1]] <- 0
  
  initial_n <- params@psi
  initial_n[,] <- 0
  w_inf_idx <- w_inf
  for (i in 1:no_sp) {
    w_inf_idx[i] <- length(w[w<=w_inf[i]])
    initial_n[i, params@species_params$w_min_idx[i]:
                (params@species_params$w_min_idx[i]+
                   (w_inf_idx[1]-params@species_params$w_min_idx[1]))] <-
      n_exact[params@species_params$w_min_idx[1]:
                (params@species_params$w_min_idx[1]+
                   (w_inf_idx[1]-params@species_params$w_min_idx[1]))] *
      (w_min[1]/w_min[i])^lambda
  }
  
  v <- sqrt(min(w_mat)*max(w_mat))
  v_idx <- length(w[w<v])
  #n_output <- initial_n*(kappa*w[v_idx]^(-lambda))/sum(sqrt(initial_n[,v_idx]*initial_n[,v_idx+19]))
  n_output <- initial_n*(kappa*w[v_idx]^(-lambda))/sum(initial_n[,v_idx])
  
  ######################################
  
  # ----
  #' ### Setup plankton
  
  plankton_vec <- (kappa*w^(-lambda))-colSums(n_output)
  plankton_vec[plankton_vec<0] <- 0
  ## need to put a waning here about -ve entries and maybe cutoff
  plankton_vec[min(which(plankton_vec==0)):length(plankton_vec)] <- 0
  params@cc_pp[sum(params@w_full<=w[1]):length(params@cc_pp)] <- plankton_vec
  initial_n_pp <- params@cc_pp
  # The cc_pp factor needs to be higher than the desired steady state in
  # order to compensate for predation mortality
  m2_background <- getM2Background(params, n_output, initial_n_pp)
  params@cc_pp <- (1+m2_background/params@rr_pp) * initial_n_pp
  
  # ----
  #' ### Setup background death and steplike psi
  
  m2 <- getM2(params, n_output, initial_n_pp)
  
  for (i in 1:no_sp) {
    params@psi[i, ] <- (w/w_inf[i])^(1-n)
    params@psi[i, w < (w_mat[i]-1e-10)] <- 0
    params@psi[i, w > (w_inf[i]-1e-10)] <- 1
    params@mu_b[i, ] <- mu0 * w^(n-1) - m2[i,]
    ## need to put a waning here about -ve entries and maybe cutoff
    params@mu_b[i,params@mu_b[i, ]<0] <- 0
  }
  
  # ----
  #' ### Set erepro to meet boundary condition
  
  rdi <- getRDI(params, n_output, initial_n_pp)
  gg <- getEGrowth(params, n_output, initial_n_pp)
  effort <- 0
  mumu <- getZ(params, n_output, initial_n_pp, effort = effort)
  erepro_final <- rdi
  for (i in (1:no_sp)){
    #  erepro_final[i] <- erepro*(gg[i,params@species_params$w_min_idx[i]]*n_output[i,params@species_params$w_min_idx[i]])/
    #    rdi[i]
    gg0 <- gg[i,params@species_params$w_min_idx[i]]
    mumu0 <- mumu[i,params@species_params$w_min_idx[i]]
    DW <- params@dw[params@species_params$w_min_idx[i]]
    erepro_final[i] <- erepro*(n_output[i,params@species_params$w_min_idx[i]]*(gg0+DW*mumu0))/rdi[i]
  }
  
  params@species_params$erepro <- erepro_final
  
  
  ## params@srr <- function(rdi, species_params) {return(rdi)}
  ## need to create a new params slots to hold n_output, and initial_n_pp
  return(params)
}

#19 have got a basic layout for the set_scaling_model function, but still need to coerce the
# grid points, include the n_output and initial_n_pp slots, include -ve warnings, add 
# fishing, and rmax. Next I will do add_species, then work on the system with background + 
# red mullet and make
