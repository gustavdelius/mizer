#' Retunes abundance multipliers of background species so aggregate abundance is a 
#' power law.
#'
#' If N_i(w) is a steady state of the McKendrik-von Foerster equation with fixed growth 
#' and death rates, then A_i*N_i(w) is also a steady state, where A_i is an abundance multiplier.  
#' When we add a foreground species to our model, we want to choose new abundance multipliers of the 
#' background species so that we the abundance, summed over all species and background 
#' resources, is close to cc(w), which is the aggregate abundance of all but 
#' the last species. We are assuming this last species is newly added, with A_i=1.  
#' 
#' retune_abundance operates of a params object, with a slot A. If i is a background 
#' species, then A_i=NA, indicating we are allowed to retune the abundance 
#' multiplier.
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
  # we are assuming that the abundance multiplier of the largest backgroud
  # species should be held fixed at 1, if though it was initially an NA.
  # to represent this we make a new background indicator A2
  A2 <- params@A
  A2[largest_background] <- 1
  # we make a list L of species we will vary the abundance parameters of
  # (everything but largest background)
  #! currently the code relies on L being a list, but we could switch 
  # to the TRUE/FALSE convention.
  L <- (1:no_sp)[is.na(A2)]
  # Determine the indices of the limits we shall integrate between 
  idx_start <- sum(params@w <= min(params@species_params$w_mat))
  idx_stop <- sum(params@w < max(params@species_params$w_inf))
  # The problem is to vary the abundance multupliers of species in L
  # so that, between the limits, to sum of the abundances 
  # of the species is "close" to the power law. 
  # More precisely, we find the abundance multipliers A_i so that 
  # the integral of the square of the relative distance (sum_{i not in L} A_i*N_i(w) + sum_{i not in L} N_i(w) - c(w))/c(w)
  # over w, between our limits, is minimized. Here c(w) is the sum of the abundances 
  # of all but the last (newly added) species, and L is the set of all background 
  # species, except the largest (that we keep the abundance multipliers of fixed).
  #
  #! how to define cc in general ? should it be smoothed ?
  # cc used to be defined as
  # cc <- colSums(params@initial_n[all_background, ])
  # but now we are assuming that the newly added species 
  # is the last one, and we are retunning to 
  # get similar to the aggregate abundance of the 
  # species (1:(no_sp-1))
  cc <- colSums(params@initial_n[1:(no_sp-1), ])
  # rho is the total abundance of all the species that have their abundance multipliers
  # held fixed.
  Lcomp <- (1:no_sp)[!is.na(A2)]
  rho <- colSums(params@initial_n[Lcomp, ])
  # We solve a linear system to find the abundance multipliers, first we initialize 
  # the matrix RR and vector QQ
  RR <- matrix(0, nrow = length(L), ncol = length(L))
  QQ <- (1:length(L))
  # Next we fill out the values of QQ and RR
  for (i in (1:length(L))) {
    QQ[i] <-
      sum((params@initial_n[L[i], ] * (cc - rho) * params@dw / (cc ^ 2))[idx_start:idx_stop])
    for (j in (1:length(L))) {
      RR[i, j] <-
        sum((
          params@initial_n[L[i], ] * params@initial_n[L[j], ] * params@dw / (cc ^ 2)
        )[idx_start:idx_stop])
    }
  }
  # Now we solve the linear system to find the abundance multipliers that 
  # yield our power law
  A2[L] <- solve(RR, QQ)
  if (sum(A2 < 0) > 0) {
    stop('Abundance multipliers generated with negative entries')
    #! we should add an extra iteration to solve this issue of -ve 
    # abundance multipliers by holding certain species off.
  }
  return(A2)
}


################################################################################

#' Adds a new species into the system, and sets its abundance to the steady state 
#' in the system where the new species does not self interact. Then the abundance 
#' multipliers of the background species are retuned to retain the old aggregate 
#' abundance curve, using retune_abundance().
#' 
#' Note that we assuming that the first species is a background species, and the
#' last species is a foreground species, with abundance multiplier mult. 
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

#! need to add in example code in help

#! assumption of first species being background.

add_species <- function(params, species_params, mult = 1.5 * 10 ^ (11)) {
  # create large r_max's if such slots are absent
  #! is it correct to make the rmax's like this
  if (is.null(params@species_params$r_max)){
    params@species_params$r_max <- params@species_params$w_inf
    params@species_params$r_max[] <- 10^50
  }
  if (is.null(species_params$r_max)){
    species_params$r_max <- 10^50
  }
  # add the new species (with parameters described by species_params), 
  # to make a larger species_params dataframe.
  combi_species_params <- rbind(params@species_params, species_params)
  # use dataframe and global settings from params to make a new MizerParams 
  # object.
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
  # Use the same resource specrum as params
  combi_params@initial_n_pp <- params@initial_n_pp
  combi_params@cc_pp <- params@cc_pp
  new_sp <- length(params@species_params$species) + 1
  # Initially use abundance curves for pre-existing species 
  # (we shall retune the abundance multipliers of such 
  # species from the background later)
  combi_params@initial_n[1:(new_sp - 1), ] <- params@initial_n
  # Use the same psi and mu_b as before for old species
  combi_params@psi[1:(new_sp - 1), ] <- params@psi
  combi_params@mu_b[1:(new_sp - 1), ] <- params@mu_b
  #! maybe we do not have to change psi[new_sp,]
  combi_params@psi[new_sp, ] <-
    (combi_params@w / combi_params@species_params$w_inf[new_sp]) ^ (1 - combi_params@n)
  combi_params@psi[new_sp, combi_params@w < (combi_params@species_params$w_mat[new_sp] - 1e-10)] <-
    0
  combi_params@psi[new_sp, combi_params@w > (combi_params@species_params$w_inf[new_sp] - 1e-10)] <-
    1
  # combi_params@mu_b[new_sp, ] <- params@mu_b[(new_sp - 1), ]
  combi_params@mu_b[new_sp, ] <- params@mu_b[1, ]
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

# #20 #42 cc is now cc <- colSums(params@initial_n[1:(no_sp-1), ])
# cc used to be defined as
# cc <- colSums(params@initial_n[all_background, ])
# but now we are assuming that the newly added foreground species 
# is the last one, and we are retuneing to 
# get similar to the aggregate abundance of the 
# species (1:(no_sp-1))

# #20 #42 Finished cleaning up and commenting retune_species, except for three issues
# (marked with #! in code): 1 currently the code relies on L being a list, but we could switch 
# to the TRUE/FALSE convention. 2   How to define cc in general ? should it be smoothed ? 
# should it relate to the original power law system, or just the last constructed system.
# 3 we should add an extra iteration to solve this issue of -ve 
# abundance multipliers by holding certain species off.

# #20 #42 Started writing help for add_species

# #20 #42 Progress on add_species clean up. Discussing erepro less
