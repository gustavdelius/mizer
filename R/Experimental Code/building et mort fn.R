library(mizer)
library(plyr)


source("R/Experimental Code/soldirectsimplified.R")




#betaval <- res@species_params$beta[1]
#lambdaval <- lambda
#sigmaval <- sigma[1]
#sigmaval <- res@species_params$sigma[1]
#wvec <- res@w
#gammaval <- res@species_params$gamma[1]
#WW <- max(object$w_inf)
#fbar <- f0
#kappaval <- kappa
#qval <- q

prm <- list("w" = wvec,"beta"=beta,"lambda"=lambda,
            "sigma"=sigmaval,"gamma"=gamma, "WW"=max(wvec),"fbar"=fbar,"kappa"=kappaRval,
            "q"=qval,"n"=nval)
wvec <- prm$w
betaval <- prm$beta
#lambdaval <- prm$lambda
sigmaval <- prm$sigma
gammaval <- prm$gamma
WW <- prm$WW
fbar <- prm$fbar
kappaval <- prm$kappa
qval <- prm$q
nval <- prm$n
lambdaval <- 2+qval-nval
ans <- (-betaval^(1 - lambdaval + qval))*exp((1/2)*((1 - lambdaval + qval)^2)*sigmaval^2)*(-1 + fbar)*gammaval*kappaval*sqrt(pi/2)*sigmaval*wvec^(1 - lambdaval + qval)*(sqrt(1/sigmaval^2)*sigmaval + erf(((1 - lambdaval + qval)*sigmaval^2 + log(betaval) + log(wvec) -log(WW))/(sqrt(2)*sigmaval)))

plot(wvec,ans,log="x")

et_mort_any <- function(prm){
  wvec <- prm$w
  betaval <- prm$beta
  #lambdaval <- prm$lambda
  sigmaval <- prm$sigma
  gammaval <- prm$gamma
  WW <- prm$WW
  fbar <- prm$fbar
  kappaval <- prm$kappa
  qval <- prm$q
  nval <- prm$n
  lambdaval <- 2+qval-nval
  ans <- (-betaval^(1 - lambdaval + qval))*exp((1/2)*((1 - lambdaval + qval)^2)*sigmaval^2)*(-1 + fbar)*gammaval*kappaval*sqrt(pi/2)*sigmaval*wvec^(1 - lambdaval + qval)*(sqrt(1/sigmaval^2)*sigmaval + erf(((1 - lambdaval + qval)*sigmaval^2 + log(betaval) + log(wvec) -log(WW))/(sqrt(2)*sigmaval)))
  return(ans)
}

plot(wvec,et_mort_any(prm),log="x")

et_mort_num <- function(prm){
  wvec <- prm$w
  betaval <- prm$beta
  #lambdaval <- prm$lambda
  sigmaval <- prm$sigma
  gammaval <- prm$gamma
  WW <- prm$WW
  fbar <- prm$fbar
  kappaval <- prm$kappa
  qval <- prm$q
  nval <- prm$n
  lambdaval <- 2+qval-nval
  logdiff <- log(wvec)[2] -log(wvec)[1]
  bigenough <- wvec[length(wvec)]*betaval*exp(5*sigmaval)
  wvecright <- exp(seq(log(wvec[length(wvec)]),log(bigenough),by=logdiff))
  Lr <- length(wvecright)-1
  
  mort_eter <- function(w){
    return(sum(((1-fbar)*exp(-(((log(wvecright[1:Lr])-log(w)-log(betaval))^2)/(2*sigmaval*sigmaval))
    )*gamma*kappaval*wvecright[1:Lr]^(qval-lambda))*(wvecright[2:(Lr+1)]-wvecright[1:Lr])))
  }
  ans <- sapply(wvec,mort_eter)
 
  return(ans)
}

lines(wvec,et_mort_num(prm))

source("R/MizerParams-classETRN.R")
source("R/wrapper_functionsETRN.R")

paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10,n=2/3,q=0.8,eta=0.25,
                               k0=10^(50),kappa=0.4, alpha=0.17,
                               h=1614.363,f0=0.5,ks=0,z0=0,gamma=660.2633
)
lines(paramsConst@w,paramsConst@mu_et) 

plot(paramsConst@w,paramsConst@mu_et,log="x") 

plot(wvec,et_mort_any(prm),log="x")
lines(wvec,et_mort_num(prm))
lines(paramsConst@w,paramsConst@mu_et) 

plot(wvec,et_mort_num(prm),log="x")
lines(wvec,et_mort_num(prm))
lines(paramsConst@w,paramsConst@mu_et) 

lines(wvec,etmort)

