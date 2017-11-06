library(mizer)
library(plyr)
paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10)
simConst <- project(paramsConst, t_max=150, effort = 0)
plot(simConst)

paramsConst@species_params$h

fbar <- 0.4
hbar <- (paramsConst@species_params$alpha*paramsConst@species_params$h*fbar)[1]
# perhaps better to let hbar=10 ?
#paramsConst@species_params$k
myspno <- 5
psi <- paramsConst@psi[myspno,]
nval <- 2/3
wvec <- paramsConst@w
growth_rate <- (1-psi)*hbar*(wvec^nval)  

qval <- 0.9
lambdastar <- 2+qval-nval
chi <- 0.05
lambda <- ((lambdastar+1)/(chi+1))-1
thetabar <- 1
beta <- paramsConst@species_params$beta[1]
sigmaval <- paramsConst@species_params$sigma[1] 
kappaRval <- 0.005
gammaval <- paramsConst@species_params$gamma[1]
alpha_p <- thetabar*(1-fbar)*sqrt(2*pi)*kappaRval*gammaval*sigmaval*(beta^(1+qval-lambda))*
  exp(((sigmaval^2)*(1+qval-lambda)^2)/2)
alpha_p
mu_p <- alpha_p*wvec^(nval-1)

integrand <- mu_p/(growth_rate^(chi+1))

determine_int <- function(w){
  L <- length(wvec[wvec<w])
  within <- sum(((wvec[2:L+1]-wvec[1:L])*integrand[1:L]))
  return(within)
}

intvals <- sapply(wvec,determine_int)
intvals2 <- rep(2,length(wvec))
intvals2[1:72] <- intvals[1:72]

result_n <- function(C){
  return(((chi*(intvals2+C))^(-1/chi))/growth_rate)
}
plot(wvec[1:71],result_n(1)[1:71],log="xy")

repro_eff <- 0.1
egg_size <- wvec[1]

is_zero <- function(C){
  result_nn <- result_n(C)
  result_nnn <- rep(0,length(wvec))
  result_nnn[1:71] <- result_nn[1:71]
  LHS <- result_nnn[1]*growth_rate[1]
  RHS <- sum(((result_nnn*psi*hbar*wvec^nval)[1:(length(wvec)-1)])*(wvec[2:length(wvec)]-wvec[1:(length(wvec)-1)]))*repro_eff/(2*egg_size)
  return(RHS-LHS)
  
}

plot(10:20,sapply(10:20,is_zero))

plot(wvec[1:71],result_n(11.5)[1:71]*wvec[1:71],log="xy")
abline(v=149)
