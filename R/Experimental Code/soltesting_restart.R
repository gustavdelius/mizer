library(mizer)
library(plyr)
paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10)
simConst <- project(paramsConst, t_max=150, effort = 0)
plot(simConst)

kappaRval <- 0.005
kappaRstarval <- 0.004


myspno <- 5
nval <- 2/3
qval <- 0.9
lambdastar <- 2+qval-nval
chi <- 0.05
beta <- paramsConst@species_params$beta[1]
sigmaval <- paramsConst@species_params$sigma[1] 

# determine phi
wvec <- paramsConst@w

determine_phi <- function(w){
  s <- exp(-(log(beta*wvec/w))^2/(2*sigmaval*sigmaval))
  integrand <- s*wvec*kappaRstarval*wvec^(-lambdastar)
  LL <- length(wvec)-1
  riemann <- sum(((wvec[2:(LL+1)]-wvec[1:LL])*integrand[1:LL]))
  return(riemann)
}
#determine_phi(1)

phivals <- sapply(wvec,determine_phi)
gamma <- paramsConst@species_params$gamma[1]
h <- paramsConst@species_params$h[1]
Ee <- gamma*phivals*wvec^qval
feeding_lvl <- Ee/(Ee+h*wvec^nval)
plot(wvec,feeding_lvl,log="x")
# hopefully dip at LHS is caused because finite w resolution cannot 
# account for beta
fbar <- feeding_lvl[length(feeding_lvl)-2]

psi <- paramsConst@psi[myspno,]
alpha <- paramsConst@species_params$alpha

growth_rate <- (alpha*fbar*h*wvec^nval)*(1-psi)

################

alpha_p <- (1-fbar)*sqrt(2*pi)*kappaRval*gamma*sigmaval*(beta^(1+qval-lambda))*
  exp(((sigmaval^2)*(1+qval-lambda)^2)/2)
alpha_p
mu_p <- alpha_p*wvec^(nval-1)

# check this vs numeric integral

determine_mort <- function(w){
  s <- exp(-(log(beta*w/wvec))^2/(2*sigmaval*sigmaval))
  integrand <- (1-fbar)*s*kappaRval*(wvec^(-lambda))*gamma*wvec^qval
  #integrand <- (1-feeding_lvl)*s*kappaRval*(wvec^(-lambda))*gamma*wvec^qval
  
  LL <- length(wvec)-1
  riemann <- sum(((wvec[2:(LL+1)]-wvec[1:LL])*integrand[1:LL]))
  return(riemann)
}

mortvals <- sapply(wvec,determine_mort)

A<-1
mu_pp <- alpha_p*wvec^(A+1+qval-((lambdastar+A)/(chi+1)))


plot(wvec,mortvals,log="xy")
lines(wvec,mu_pp)
# work out mortality integral again, 
#and check whether an adjustment should be added