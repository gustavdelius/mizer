
n2RECvec <- function(params,n_pp, n){
  no_sp <- dim(n)[1]
  indx <- (params@species_params$w_min_idx-1) * no_sp + (1:no_sp)
  c(n_pp,n[indx],colSums(n))
}

RECvec2n <- function(params,X){
  n_pp <- X[1:length(params@cc_pp)]
  no_sp <- dim(params@initial_n)[1]
  n_e <- X[(length(params@cc_pp)+1):(length(params@cc_pp)+no_sp)]
  n_c <- X[(length(params@cc_pp)+no_sp+1):length(X)]
  return(list(n_pp=n_pp,n_e=n_e,n_c=n_c))
}

myF <- function(params,X){
  inputs <- RECvec2n(params,X)
  dummy_n <- params@initial_n
  dummy_n[] <- 0
  dummy_n[dim(dummy_n)[1],] <- inputs$n_c
  gg <- getEGrowth(params,dummy_n,inputs$n_pp)
  mumu <- getZ(params,dummy_n,inputs$n_pp,effort = 0)
  p <- params
  nres <- p@initial_n
  for (new_sp in (1:no_sp)){
    w_inf_idx <- sum(p@w < p@species_params$w_inf[new_sp])
    idx <- p@species_params$w_min_idx[new_sp]:(w_inf_idx-1)
    if (any(gg[idx]==0)) {
      stop("Can not compute steady state due to zero growth rates")
    }
    nres[new_sp, ] <- 0
    nres[new_sp, p@species_params$w_min_idx[new_sp]:w_inf_idx] <- 
      c(1, cumprod(gg[idx] / ((gg + mumu * p@dw)[idx+1])))
    
    if (any(is.infinite(nres))) {
      stop("Candidate steady state holds infinities")
    }
    if (any(is.na(nres) || is.nan(nres))) {
      stop("Candidate steady state holds none numeric values")
    }
    nres[new_sp,] <- nres[new_sp,]*inputs$n_e[new_sp]/(nres[new_sp,p@species_params$w_min_idx[new_sp]]) 
  }
  Y <- n2RECvec(params,inputs$n_pp, nres)
  return(X-Y)
}


#################################################################
#################################################################

# Compute solution for new species
w_inf_idx <- sum(p@w < p@species_params$w_inf[new_sp])
idx <- p@species_params$w_min_idx[new_sp]:(w_inf_idx-1)
if (any(gg[idx]==0)) {
  stop("Can not compute steady state due to zero growth rates")
}
p@initial_n[new_sp, ] <- 0
p@initial_n[new_sp, p@species_params$w_min_idx[new_sp]:w_inf_idx] <- 
  c(1, cumprod(gg[idx] / ((gg + mumu * p@dw)[idx+1])))
if (any(is.infinite(p@initial_n))) {
  stop("Candidate steady state holds infinities")
}
if (any(is.na(p@initial_n) || is.nan(p@initial_n))) {
  stop("Candidate steady state holds none numeric values")
}

# Normalise solution
if (is.na(SSB)) {
  # If spawning stock biomass of new species is not supplied, 
  # normalise solution so that at its maximum it lies at half the 
  # power law, and then calculate its SSB.
  # We choose the maximum of the biomass density in log space
  # because that is always an increasing function at small size.
  idx <- which.max(p@initial_n[new_sp, ] * p@w^p@lambda)
  p@initial_n[new_sp, ] <- p@initial_n[new_sp, ] *
    p@kappa * p@w[idx]^(-p@lambda) / p@initial_n[new_sp, idx] / 2
  SSB <- sum(p@initial_n[new_sp, ] * p@w * p@dw * p@psi[new_sp, ])
} else {
  unnormalised_SSB <- sum(p@initial_n[new_sp,] * p@w * p@dw * 
                            p@psi[new_sp, ])
  p@initial_n[new_sp, ] <- p@initial_n[new_sp, ] * SSB / unnormalised_SSB
}
p@A <- c(params@A, SSB)
