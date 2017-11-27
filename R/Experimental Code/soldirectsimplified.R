library(mizer)
library(plyr)
paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10)
kappaRval <- 0.4
qval <- 0.8
myspno <- 3
nval <- 2/3
lambda <- 2+qval-nval
beta <- paramsConst@species_params$beta[1]
sigmaval <- paramsConst@species_params$sigma[1] 
wvec <- paramsConst@w
wvecfull <- paramsConst@w_full
determine_phi <- function(w){
  s <- exp(-(log(beta*wvecfull/w))^2/(2*sigmaval*sigmaval))
  integrand <- s*wvecfull*kappaRval*wvecfull^(-lambda)
  LL <- length(wvecfull)-1
  riemann <- sum(((wvecfull[2:(LL+1)]-wvecfull[1:LL])*integrand[1:LL]))
  return(riemann)
}
phivals <- sapply(wvec,determine_phi)
gamma <- paramsConst@species_params$gamma[1]
Ee <- gamma*phivals*wvec^qval
fbar <- 0.5
h <- (((Ee/fbar)-Ee)/(wvec^nval))[1]
psi <- paramsConst@psi[myspno,]
# check psi has appropriate form
#alpha <- paramsConst@species_params$alpha[myspno]
alpha <- 0.17
k <- 0
hbar <- alpha*fbar*h -k
growth_rate <- (hbar*wvec^nval)*(1-psi)

############
determine_mort <- function(w){
  s <- exp(-(log(beta*w/wvecfull))^2/(2*sigmaval*sigmaval))
  integrand <- (1-fbar)*s*kappaRval*(wvecfull^(-lambda))*gamma*wvecfull^qval
  LL <- length(wvecfull)-1
  riemann <- sum(((wvecfull[2:(LL+1)]-wvecfull[1:LL])*integrand[1:LL]))
  return(riemann)
}
mortvals <- sapply(wvec,determine_mort)
# work out mortvals alternative way, then compare
maxsize <- paramsConst@species_params$w_inf[myspno]
maxsizeindex <- length(maxsize[wvec<maxsize])
realalpha_p <- max((mortvals/(wvec^(1+1+qval-((lambda+1)/(1)))))[1:maxsizeindex])
realmu_pp <- realalpha_p *wvec^(1+1+qval-((lambda+1)/(1)))
plot(wvec,mortvals,log="xy")
lines(wvec,realmu_pp)
abline(v=maxsize)
# solve for steady state
solintegrand <- realmu_pp/(growth_rate)
determine_int <- function(w){
  L <- length(wvec[wvec<w])
  within <- sum(((wvec[2:(L+1)]-wvec[1:L])*solintegrand[1:L]))
  return(within)
}
solintvals <- sapply(wvec[1:maxsizeindex],determine_int)

result_chi0_n <- function(H){
  npart <- exp(-solintvals)*H/growth_rate[1:maxsizeindex]
  nstart <- rep(0,length(wvec))
  nstart[1:maxsizeindex] <- npart
  return(nstart)
}
n1 <- result_chi0_n(1)
y <- wvec
matsize <- paramsConst@species_params$w_mat[myspno]
egg_size <- wvec[1]
kappaintegrand <- matsize*y^((lambda+1)*(1))*n1^(1)*y^(-2)
LL <- length(wvec)-1
kappaoutput <- sum(kappaintegrand[1:LL]*(wvec[2:(LL+1)]-wvec[1:LL]))
HH <- kappaRval/kappaoutput
nvgood <- result_chi0_n(HH)
LHS <- nvgood[1]*growth_rate[1]
RHSS <- sum(((nvgood*psi*hbar*wvec^nval)[1:(length(wvec)-1)])*(wvec[2:length(wvec)]-wvec[1:(length(wvec)-1)]))/(2*egg_size)
repro_efff <- LHS/RHSS
repro_efff