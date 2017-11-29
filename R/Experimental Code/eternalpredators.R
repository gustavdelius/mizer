

library(pracma)
source("R/Experimental Code/soldirectsimplified.R")
logdiff <- log(wvec)[2] -log(wvec)[1]
bigenough <- wvec[length(wvec)]*beta*exp(5*sigmaval)
wvecright <- exp(seq(log(wvec[length(wvec)]),log(bigenough),by=logdiff))
Lr <- length(wvecright)-1

mort_eter <- function(w){
  return(sum(((1-fbar)*exp(-(((log(wvecright[1:Lr])-log(w)-log(beta))^2)/(2*sigmaval*sigmaval))
  )*gamma*kappaRval*wvecright[1:Lr]^(qval-lambda))*(wvecright[2:(Lr+1)]-wvecright[1:Lr])))
}
plot(wvec,sapply(wvec,mort_eter),log="x")
etmort <- realmu_pp-mortvals
lines(wvec,etmort)

WW <- 10^5

direct_et_mort <- function(w){
const <- (-beta^(1 - lambda + qval))*exp((1/2)*((1 - lambda + qval)^2)*sigmaval^2)  
const*(
  -1 + fbar)*gamma*kappaRval*sqrt(pi/2)*sigmaval*w^(1 - lambda + qval)*(sqrt(1/sigmaval^2)*sigmaval + erf(((1 - lambda + qval)*sigmaval^2 + log(beta) + log(w) - 
                                                                                             log(WW))/(sqrt(2)*sigmaval)))
}

plot(wvec,sapply(wvec,direct_et_mort),log="x")

const <- (-beta^(1 - lambda + qval))*exp((1/2)*((1 - lambda + qval)^2)*sigmaval^2)  


(-beta^(1 - lambda + qval))*exp((1/2)*(1 - lambda + qval)^2*sigmaval^2)
#etmortgood[w_] := (-beta^(1 - lambda + q))*
#  E^((1/2)*(1 - lambda + q)^2*sigma^2)*(-1 + fbar)*gamma*kappa*
#  Sqrt[Pi/2]*sigma*w^(1 - lambda + q)*
#  (Sqrt[1/sigma^2]*sigma + 
#     Erf[((1 - lambda + q)*sigma^2 + Log[beta] + Log[w] - 
#            Log[W])/(Sqrt[2]*sigma)])

  
  (Sqrt[1/sigma^2]*sigma + Erf(((1 - lambda + q)*sigma^2 + log(beta) + log(w) - 
                 log(W))/(sqrt(2)*sigma)))
  
  
# implement and plot the above
# where variables go through mizer
# grid point maker
# increasingly accurate 
