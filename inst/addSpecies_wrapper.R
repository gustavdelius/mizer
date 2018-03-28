#' Retunes abundance multipliers of background species so aggregate abundance is a 
#' power law.
#'
#' If N_i(w) is a steady state of the McKendrik-von Foerster equation with fixed growth 
#' and death rates, then A_i*N_i(w) is also a steady state, where A_i is an abundance multiplier.  
#' When we add a foreground species to our model, we want to choose new abundance multipliers of the 
#' background species so that we the abundance, summed over all species and background 
#' resources, is close to the power law kappa*w^(-lambda). 
#' 
#' retune_abundance operates of a params object, with a slot A. If i is a background 
#' species, then A[i]=NA, indicating we are allowed to retune the abundance 
#' multiplier  
#' 
#'
#' @param params A mizer params object with an A slot with 1's for species 
#' we wish to hold fixed the abundance multiplier of, and NA's for species that 
#' we shall vary the abundance multiplier of.

#' @export
#' @return An object of type \code{MizerParams}
#' @seealso \linkS4class{MizerParams}
#' @examples
#' \dontrun{
#' plotSpectra(sim)
#' }

retune_abundance <- function(params) {
  no_sp <- length(params@species_params$w_inf)
  # get a list of the species that we can tune abundance mulpliers of
  all_background <- is.na(params@A)
  largest_background <-
    which.max(params@species_params$w_inf[all_background])
  # we are assuming that the the abundance multiplier of the largest backgroud
  # species should be held fixed at 1, if though it was initially an NA.
  A2 <- params@A
  A2[largest_background] <- 1
  # we make a list L of species we will vary the abundance parameters of
  # (everything but largest background)
  L <- (1:no_sp)[is.na(A2)]
  idx_start <- sum(params@w <= min(params@species_params$w_mat))
  idx_stop <- sum(params@w < max(params@species_params$w_inf))
  RR <- matrix(0, nrow = length(L), ncol = length(L))
  QQ <- (1:length(L))
  Lcomp <- (1:no_sp)[!is.na(A2)]
  old_n <- params@initial_n
  cc <- colSums(params@initial_n[all_background, ])
  rho <- colSums(params@initial_n[Lcomp, ])
  den <- cc ^ 2
  for (i in (1:length(L))) {
    QQ[i] <-
      sum((params@initial_n[L[i], ] * (cc - rho) * params@dw / (den))[idx_start:idx_stop])
    for (j in (1:length(L))) {
      RR[i, j] <-
        sum((
          params@initial_n[L[i], ] * params@initial_n[L[j], ] * params@dw / (den)
        )[idx_start:idx_stop])
    }
  }
  A2[L] <- solve(RR, QQ)
  if (sum(A2 < 0) > 0) {
    stop('Abundance multipliers generated with negative entries')
  }
  return(A2)
}


################################################################################

add_species <- function(params, species_params, mult = 1.5 * 10 ^ (11)) {
  
  combi_species_params <- rbind(params@species_params, species_params)
  if (length(combi_species_params$r_max) == 0) {
    combi_species_params$r_max <- 10 ^ (50)
  }
  combi_params <-
    MizerParams(
      combi_species_params,
      p = params@p,
      n = params@n,
      q = params@q,
      lambda = params@lambda,
      f0 = params@f0,
      kappa = params@kappa,
      min_w = min(params@w),
      max_w = max(params@w),
      no_w = length(params@w),
      min_w_pp = min(params@w_full),
      w_pp_cutoff = max(params@w_full),
      r_pp = (params@rr_pp / (params@w_full ^ (params@p - 1)))[1]
    )
  
  
  combi_params@initial_n_pp <- params@initial_n_pp
  combi_params@cc_pp <- params@cc_pp
  
  new_sp <- length(params@species_params$species) + 1
  
  combi_params@initial_n[1:(new_sp - 1), ] <- params@initial_n
  combi_params@species_params$erepro[1:(new_sp - 1)] <-
    params@species_params$erepro
  combi_params@psi[1:(new_sp - 1), ] <- params@psi
  combi_params@mu_b[1:(new_sp - 1), ] <- params@mu_b
  
  # other important info to pass through correspond to
  # the parts of params that got modified after set_scaling made it initially.
  combi_params@species_params$erepro[new_sp] <- 0.1
  #! maybe we do not have to change psi[new_sp,]
  combi_params@psi[new_sp, ] <-
    (combi_params@w / combi_params@species_params$w_inf[new_sp]) ^ (1 - combi_params@n)
  combi_params@psi[new_sp, combi_params@w < (combi_params@species_params$w_mat[new_sp] - 1e-10)] <-
    0
  combi_params@psi[new_sp, combi_params@w > (combi_params@species_params$w_inf[new_sp] - 1e-10)] <-
    1
  combi_params@mu_b[new_sp, ] <- params@mu_b[(new_sp - 1), ]
  #! what about params@srr ? do I have to pass this through when rmax is off ?
  #! do I have to set rmax off if it is off in two inputs ?
  combi_params@srr <- params@srr
  # use rest of info to fill out n_new_sp properly
  # combi_params@initial_n[new_sp,] <- params@initial_n[(new_sp-1),]
  ################################
  combi_params@interaction[new_sp, new_sp] <- 0
  mumu <-
    getZ(combi_params,
         combi_params@initial_n,
         combi_params@initial_n_pp,
         effort = 0)[new_sp, ]
  gg <-
    getEGrowth(combi_params,
               combi_params@initial_n,
               combi_params@initial_n_pp)[new_sp, ]
  
  w_inf_idx <-
    sum(combi_params@w < combi_params@species_params$w_inf[new_sp])
  #! Alter this code here, so it avoids division by zero, in stunted growth case.
  integrand <-
    params@dw[combi_params@species_params$w_min_idx[new_sp]:w_inf_idx] * mumu[combi_params@species_params$w_min_idx[new_sp]:w_inf_idx] /
    gg[combi_params@species_params$w_min_idx[new_sp]:w_inf_idx]
  combi_params@initial_n[new_sp, ] <- 0
  
  combi_params@initial_n[new_sp, combi_params@species_params$w_min_idx[new_sp]:w_inf_idx] <-
    mult * exp(-cumsum(integrand)) / gg[combi_params@species_params$w_min_idx[new_sp]:w_inf_idx]
  
  combi_params@interaction[new_sp, new_sp] <- 1
  
  ##################
  #sim <- project(combi_params, t_max=5, effort = 0)
  #plot(sim)
  
  ###########
  #A <- rep(NA,length(combi_params@species_params$w_inf))
  #A <- c(params@A, 1)
  combi_params@A <- c(params@A, 1)
  A <- combi_params@A
  AA <- retune_abundance(combi_params)
  print(AA)
  new_n <- combi_params@initial_n
  for (i in 1:length(combi_params@species_params$w_inf)) {
    new_n[i, ] <- AA[i] * combi_params@initial_n[i, ]
  }
  combi_params@initial_n <- new_n
  
  
  
  
  for (i in 1:new_sp){
    mumu <-
      getZ(combi_params,
           combi_params@initial_n,
           combi_params@initial_n_pp,
           effort = 0)[i, ]
    gg <-
      getEGrowth(combi_params,
                 combi_params@initial_n,
                 combi_params@initial_n_pp)[i, ]
    
    gg0 <- gg[combi_params@species_params$w_min_idx[i]]
    mumu0 <- mumu[combi_params@species_params$w_min_idx[i]]
    rdi <-
      getRDI(combi_params,
             combi_params@initial_n,
             combi_params@initial_n_pp)[i]
    DW <-
      combi_params@dw[combi_params@species_params$w_min_idx[i]]
    #combi_params@species_params$erepro[i] <- combi_params@species_params$erepro[i]*(
    #  combi_params@initial_n[i,combi_params@species_params$w_min_idx[i]]*(gg0+DW*mumu0))/rdi
    
    H <- rdi / combi_params@species_params$erepro[i]
    X <-
      combi_params@initial_n[i, combi_params@species_params$w_min_idx[i]] *
      (gg0 + DW * mumu0)
    
    combi_params@species_params$erepro[i] <-
      combi_params@species_params$r_max[i] * X / (H * combi_params@species_params$r_max[i] +
                                                    H * X)
    
  }
  
  return(combi_params)
}

############### add mullet to set scaling ###############

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

# make gamma slot, then the combination should work. Note that f0 is not setup properly.

#h_gurn <- 3*params_data_NS$k_vb[gurn_sp]*((params_data_NS$w_inf[gurn_sp])^(1/3))/(alpha[rep_idx]*f0)

params_out <- add_species(params, species_params, mult = 5.5 * 10 ^ (8))


sim <- project(params_out, t_max = 5, effort = 0)
plot(sim)
species_params$species <- "mullet2"
species_params$w_inf <- 50

params_out_2 <- add_species(params_out, species_params, mult = 5.5 * 10 ^ (8))

sim <- project(params_out_2, t_max = 5, effort = 0)
plot(sim)

# #20 #42 Cleaning up the code for wrapper. Currently adding two mullet like species

# #20 #42 Removed A from retune_abundance