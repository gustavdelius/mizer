library(mizer)
library(plyr)
# load parameter values
paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10)

# set q and input kappas
##kappaRval <- 2.184406e-01

kappaRval <- 0.4


kappaRstarval <- kappaRval

#kappaRstarval <- 1.141175e-01
qval <- 0.9

  kappaRstarval <- kappaRval
  
  # focus on third species
  myspno <- 3
  # choose n
  nval <- 2/3
  # determine prey community exponent
  lambdastar <- 2+qval-nval
  # choose prey switching exponent
  chi <- 0.5
  chi <- 0
  # setup preference function
  beta <- paramsConst@species_params$beta[1]
  sigmaval <- paramsConst@species_params$sigma[1] 
  
  # get lists of appropriate weights
  wvec <- paramsConst@w
  wvecfull <- paramsConst@w_full
  # determine phi=available energy
  determine_phi <- function(w){
    s <- exp(-(log(beta*wvecfull/w))^2/(2*sigmaval*sigmaval))
    integrand <- s*wvecfull*kappaRstarval*wvecfull^(-lambdastar)
    LL <- length(wvecfull)-1
    riemann <- sum(((wvecfull[2:(LL+1)]-wvecfull[1:LL])*integrand[1:LL]))
    return(riemann)
  }
  # compute phi for different weights
  phivals <- sapply(wvec,determine_phi)
  gamma <- paramsConst@species_params$gamma[1]
  # determine energy encountered 
  Ee <- gamma*phivals*wvec^qval
  # set constant feeding level = 1/2
  fbar <- 0.5
  #fbar <- 0.8
  
  # determine the max intake prefactor h value that returns this feeding level
  h <- (((Ee/fbar)-Ee)/(wvec^nval))[1]
  # define fraction of energy diverted to reproduction
  psi <- paramsConst@psi[myspno,]
  # define fraction of energy assimilated by prey
  alpha <- paramsConst@species_params$alpha[myspno]
  alpha <- 0.19
  # define growth rate with k=0 metabolic loss
  growth_rate <- (alpha*fbar*h*wvec^nval)*(1-psi)
  # get exponent of predator spectrum
  lambda <- (lambdastar-chi)/(1+chi)
  
  # here we have no starvation mortality or background mortality
  # determine mortality rates by integrating
  determine_mort <- function(w){
    s <- exp(-(log(beta*w/wvecfull))^2/(2*sigmaval*sigmaval))
    integrand <- (1-fbar)*s*kappaRval*(wvecfull^(-lambda))*gamma*wvecfull^qval
    LL <- length(wvecfull)-1
    riemann <- sum(((wvecfull[2:(LL+1)]-wvecfull[1:LL])*integrand[1:LL]))
    return(riemann)
  }
  mortvals <- sapply(wvec,determine_mort)
  A<-1
  # get asyptotic size of our species
  maxsize <- paramsConst@species_params$w_inf[myspno]
  maxsizeindex <- length(maxsize[wvec<maxsize])
  # divide out the w dependence from mortalirt rate, to get mortality rate prefactor
  realalpha_p <- max((mortvals/(wvec^(A+1+qval-((lambdastar+A)/(chi+1)))))[1:maxsizeindex])
  # check resulting function has constant value = mortality rate prefactor 
  XX <- (mortvals/(wvec^(A+1+qval-((lambdastar+A)/(chi+1)))))[1:maxsizeindex]
  #plot(XX)
  # get powerlaw form of mortality rate
  realmu_pp <- realalpha_p *wvec^(A+1+qval-((lambdastar+A)/(chi+1)))
  # check powerlaw matches reality up to maxsize
  #plot(wvec,mortvals,log="xy")
  #lines(wvec,realmu_pp)
  #abline(v=maxsize)
  # prepare integrand from pde steady state expression
  integrand <- realmu_pp/(growth_rate^(chi+1))
  # evaluate integral for w from eggsize to W
  determine_int <- function(w){
    L <- length(wvec[wvec<w])
    within <- sum(((wvec[2:(L+1)]-wvec[1:L])*integrand[1:L]))
    return(within)
  }
  intvals <- sapply(wvec[1:maxsizeindex],determine_int)
  # select a value of the integration constant C, and use the result of this integral
  # together with the growth rate to form a steady state of MvF eqn.
  result_n <- function(C){
    npart <- ((chi*(intvals+C))^(-1/chi))/growth_rate[1:maxsizeindex]
    nstart <- rep(0,length(wvec))
    nstart[1:maxsizeindex] <- npart
    return(nstart)
  }
  
  integrand <- realmu_pp/(growth_rate)
  # evaluate integral for w from eggsize to W
  determine_int <- function(w){
    L <- length(wvec[wvec<w])
    within <- sum(((wvec[2:(L+1)]-wvec[1:L])*integrand[1:L]))
    return(within)
  }
  intvals <- sapply(wvec[1:maxsizeindex],determine_int)
  
  result_chi0_n <- function(H){
    npart <- exp(-intvals)*H/growth_rate[1:maxsizeindex]
    nstart <- rep(0,length(wvec))
    nstart[1:maxsizeindex] <- npart
    return(nstart)
  }
  ngood <- result_chi0_n(1)
  
  
  
  
  
  
  #plot(wvec,result_chi0_n(1),log="xy")
  
  # plot an example solution
  #plot(wvec,result_n(11),log="xy")
  
  repro_eff <- 0.1
  
  # suppose the egg size for this species equals mizer's usual egg size.
  egg_size <- wvec[1]
  # abreviation for growth rate prefactor
  hbar <- alpha*fbar*h
  

  Lambda <- (lambdastar+A)/(chi+1)
  # get weight at maturity
  matsize <- paramsConst@species_params$w_mat[myspno]
  # determine kappastar by 
  # integrating `prey-switching modified` abundance over
  # w* using scale tranformed N solution
  y <- wvec
  chii <- chi
  integrand <- matsize*y^(Lambda*(chii+1))*ngood^(chii+1)*y^(-2)
  LL <- length(wvec)-1
  # kappa star is the w independent part of said abundance integral
  kappastaroutput <- sum(integrand[1:LL]*(wvec[2:(LL+1)]-wvec[1:LL]))
  
  # calculating kappa is similiar to above, but with chii=0
  y <- wvec
  chii <- 0
  integrand <- matsize*y^(Lambda*(chii+1))*ngood^(chii+1)*y^(-2)
  LL <- length(wvec)-1
  kappaoutput <- sum(integrand[1:LL]*(wvec[2:(LL+1)]-wvec[1:LL]))
  
  HH <- kappaRval/kappaoutput
  nvgood <- result_chi0_n(HH)
  
  # compare kappa ins with kappa outs
  c(kappaRval,kappaRstarval)
  c(kappaoutput,kappastaroutput)
  # determine reproductive efficiency required to make boundary conditions work
  # this may overwrite failure of previous netwton raphson to tune Cval
  #result_nnn <- result_n(Cval)
  #result_nnn <- result_chi0_n(1)
  result_nnn <- ngood
  result_nnn <- nvgood 
  LHS <- result_nnn[1]*growth_rate[1]
  RHSS <- sum(((result_nnn*psi*hbar*wvec^nval)[1:(length(wvec)-1)])*(wvec[2:length(wvec)]-wvec[1:(length(wvec)-1)]))/(2*egg_size)
  repro_efff <- LHS/RHSS
  repro_efff
  
  plot(wvec,ngood,log="xy")
  
  ngoodsmallkappa <- ngood
  
  LHSsmall <- result_nnn[1]*growth_rate[1]
  RHSSsmall <- sum(((result_nnn*psi*hbar*wvec^nval)[1:(length(wvec)-1)])*(wvec[2:length(wvec)]-wvec[1:(length(wvec)-1)]))/(2*egg_size)
  LHSsmall/RHSSsmall
  