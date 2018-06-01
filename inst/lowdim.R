
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
  no_sp <- dim(params@initial_n)[1]
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

get_wts <- function(params){
  return(c(params@w_full,params@species_params$w_min,params@w)^params@lambda)
}

myF_rescaled <- function(params,XX){
  wts <- get_wts(params)
  return(wts*myF(params,XX/wts))
}


params <- set_scaling_model(no_sp = 2, no_w=100)

sim <- project(params, t_max=50, effort = 0, initial_n = params@initial_n, t_save = 1, initial_n_pp = params@initial_n_pp)
plot(sim)

FF <- function(XX){
  return(myF_rescaled(params,XX))
}


wts <- get_wts(params)
trial_data_in <- wts*n2RECvec(params, n_pp = params@initial_n_pp, n = params@initial_n)

ss <- multiroot(f = FF, start = trial_data_in, rtol = 0, atol = 1e-15, ctol = 0)

outputs <- RECvec2n(params, ss$root/wts)


RECvec2sols <- function(params,X){
  inputs <- RECvec2n(params,X)
  dummy_n <- params@initial_n
  dummy_n[] <- 0
  dummy_n[dim(dummy_n)[1],] <- inputs$n_c
  gg <- getEGrowth(params,dummy_n,inputs$n_pp)
  mumu <- getZ(params,dummy_n,inputs$n_pp,effort = 0)
  no_sp <- dim(params@initial_n)[1]
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
  return(nres)
}


sim <- project(params, t_max=15, effort = 0, initial_n = RECvec2sols(params, ss$root/wts), t_save = 1, initial_n_pp = outputs$n_pp)
plot(sim)

