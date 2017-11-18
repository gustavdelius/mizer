library(mizer)
library(plyr)
# load parameter values
paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10)

# set q and input kappas
kappaRval <- 5.273115e-11
kappaRstarval <- 2.697401e-18
qval <- 2.3
go2 <- function(kappaRval,kappaRstarval,qval){
# focus on third species
myspno <- 3
# choose n
nval <- 2/3
# determine prey community exponent
lambdastar <- 2+qval-nval
# choose prey switching exponent
chi <- 0.5
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
# determine the max intake prefactor h value that returns this feeding level
h <- (((Ee/fbar)-Ee)/(wvec^nval))[1]
# define fraction of energy diverted to reproduction
psi <- paramsConst@psi[myspno,]
# define fraction of energy assimilated by prey
alpha <- paramsConst@species_params$alpha[myspno]
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
# plot an example solution
#plot(wvec,result_n(11),log="xy")

repro_eff <- 0.1

# suppose the egg size for this species equals mizer's usual egg size.
egg_size <- wvec[1]
# abreviation for growth rate prefactor
hbar <- alpha*fbar*h

# Boundary conditions RHS(C)-LHS(C)
is_zero <- function(C){
  result_nnn <- result_n(C)
  LHS <- result_nnn[1]*growth_rate[1]
  RHS <- sum(((result_nnn*psi*hbar*wvec^nval)[1:(length(wvec)-1)])*(wvec[2:length(wvec)]-wvec[1:(length(wvec)-1)]))*repro_eff/(2*egg_size)
  return(RHS-LHS)
}
library(pracma)
# try and find C to solve BC
##Cval <- newtonRaphson(is_zero, 1)$root
##Cval <- 1.00875e+11 
Cval <- 1
#plot(wvec,result_n(Cval),log="xy")

# use this to define solution
ngood <- result_n(Cval)
# get abundance scaling exponent
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
# compare kappa ins with kappa outs
c(kappaRval,kappaRstarval)
c(kappaoutput,kappastaroutput)
# determine reproductive efficiency required to make boundary conditions work
# this may overwrite failure of previous netwton raphson to tune Cval
result_nnn <- result_n(Cval)
LHS <- result_nnn[1]*growth_rate[1]
RHSS <- sum(((result_nnn*psi*hbar*wvec^nval)[1:(length(wvec)-1)])*(wvec[2:length(wvec)]-wvec[1:(length(wvec)-1)]))/(2*egg_size)
repro_efff <- LHS/RHSS
repro_efff
return(c(kappaoutput,kappastaroutput,repro_efff))
}

itme <- function(Z,qval){
  return(go2(Z[1],Z[2],qval)[1:2])
}

Z0 <- c(5.273115e-11,2.697401e-18)
Z0 <- c(10^1,10^1)

QQval <- 2.1
T <- 500
kapparesults <- 1:T
Z0 <- itme(Z0,QQval)
for (t in (1:T)){
  kapparesults[t] <- Z0[1]
  Z0 <- itme(Z0,QQval)
}
plot(kapparesults,log="y")
# note the Z0 from iterating this is exactly self returning, although the epsilon
# is > 1
Z0
go2(Z0[1], Z0[2] ,QQval)

#go2(3.125824e-02, 3.896207e-05 ,1.5)
#go2(2.184797e-01,1.142112e-01,2.5)

go2(2.184406e-01,1.141175e-01,2.5)
itme(Z0,QQval)
go2(Z0[1],Z0[2],QQval)[1:2]

go2(0.2184406 ,0.1141175,1.9)

go2(0.2664,0.3990,2.1)
library("rootSolve")
multiroot()

QQ <- 2.1
zerme <- function(Z0){
return(itme(Z0,QQ)-Z0)
}
multiroot(zerme,c(1,1),maxiter = 10000)$root

QQvec <- seq(0.5,3,0.1)
kapstorevec <- QQvec
kapstarstorevec <- QQvec
epsilonvec <- QQvec

L <- length(QQvec)
for (i in (1:L)){
  zerme2 <- function(Z0){
    return(itme(Z0,QQvec[i])-Z0)
  }
  res <- multiroot(zerme2,c(1,1),maxiter = 10000)$root
  kapstorevec[i] <- res[1]
  kapstarstorevec[i] <- res[2]
  epsilonvec[i] <- go2(res[1],res[2],QQvec[i])[3]
  
}
itme(c(1,1),QQQ)
plot()
QQQ <- 0.5
zerme2 <- function(Z0){
  return(itme(Z0,QQQ)-Z0)
}
res <- multiroot(zerme2,itme(c(1,1),QQQ),maxiter = 10000)$root
go2(res[1],res[2],QQQ)
