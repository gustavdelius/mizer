


#@ retune_abundance <- function(params,A){
  params <- set_scaling_model()
  A <- rep(1,length(params@species_params$w_inf))
  A[4:7] <- NA

  
  
  # get a list of the species that we can tune abundance mulpliers of
  all_background <- (1:length(A))[is.na(A)]
  largest_background <- which.max(params@species_params$w_inf[all_background])
  A2 <- A
  # we are assuming that the the abundance multiplier of the largest backgroud 
  # species should be held fixed at 1, if though it was initially an NA
  A2[largest_background] <- 1
  
  
  # we make a list L of species we will vary the abundance parameters of
  # (everything but largest background)
  #L <- all_background[all_background!=largest_background]
  L <- (1:length(A))[is.na(A2)]
  idx_start <- sum(params@w<=min(params@species_params$w_mat[L]))
  idx_stop <- sum(params@w<=max(params@species_params$w_inf[L]))
  RR <- matrix(0,nrow = length(L), ncol = length(L))
  QQ <- (1:length(L))
  cc <- params@kappa*params@w^(-params@lambda)
  Lcomp <- (1:length(A))[!is.na(A2)]
  
  old_n <- params@initial_n
  no_sp <- length(params@species_params$w_inf)
  
  for (i in 1:no_sp){
    old_n[i,] <- A2[i]*params@initial_n[i,]
  }
  #rho <- colSums(params@initial_n[Lcomp,])
  rho <- colSums(old_n[Lcomp,])
  
  
  for (i in (1:length(L))){
    QQ[i] <- sum((params@initial_n[L[i],]*(cc-rho)*params@dw/(cc^2))[idx_start:(idx_stop-1)])
    for (j in (1:length(L))){
      RR[i,j] <- sum((params@initial_n[L[i],]*params@initial_n[L[j],]*params@dw/(cc^2))[idx_start:(idx_stop-1)])
    }
  }
  
  A2[L] <- solve(RR,QQ)

  #@  return(A2)}
  new_n <- params@initial_n
  for (i in 1:no_sp){
    new_n[i,] <- A2[i]*params@initial_n[i,]
  }
  A2
  plot(params@w,new_n[1,],log="xy")
  for (i in (2:no_sp)){
    lines(params@w,new_n[i,])
  }
  
  #  #20 and #42 Made steps to make sure abundance multipliers are used. 
  # Finishing building retune_abundance in a seperate file. 
  # 
  
  #  #20 and #42 Have debugged a few little things. The resulting code (retune_abundance.R
  # in the adsp branch)
  # is tested on the scale free case, with some multipliers as 1, and 
  # some multipliers as NA. It does not work correctly,  A[5:7] = c(1.4,0.27,2.28)
  # I expected instead that the code would set these multipliers to unity.