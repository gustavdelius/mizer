


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
