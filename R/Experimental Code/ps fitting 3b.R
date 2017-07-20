# I need to modify this code so it has time varying growth rates, and so that 
# it is run dynamically, with changing fishing effort.

library(devtools)
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library(mizer)
library(deSolve)
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
mytimes <- seq(0,100,0.1)
param1 <- MizerParams(params_data, interaction = inter, no_w = 100)
test<-project(param1,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)

hdef <- param1@species_params[,14]
gammadef <- param1@species_params[,15]
w_infdef <- params_data[["w_inf"]]

sizeHistory <- function(param_data, times=seq(0,1,0.1), h, gamma, w_inf=w_infdef ){
  params_dataB <- param_data
  
  params_dataB[["h"]] <- h
  params_dataB[["gamma"]] <- gamma
  params_dataB[["w_inf"]] <- w_inf
  
  param1b <- MizerParams(params_dataB, interaction = inter, no_w = 100)
  
  test<-project(param1b,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)
  ns <- dim(test@n)[2]
  gg <- getEGrowth(param1b, test@n[dim(test@n)[1],,], test@n_pp[dim(test@n_pp)[1],])
  sols <- as.list(1:ns)
  MM <- matrix(0,nrow=length(times),ncol=dim(test@n)[2])
  for (i in (1:ns)){
    gini <- approxfun(param1b@w, gg[i,])
    
    myodefun <- function(t, state, parameters){
      return(list(gini(state)))
    }
    out <- ode(y = param1@w[1], times = times, func = myodefun, parms = 1)
    MM[,i] <- out[,2]
  }
  return(MM)
}

plot(mytimes,sizeHistory(params_data,mytimes, hdef, gammadef )[,4])
points(mytimes,sizeHistory(params_data,mytimes, hdef+50, gammadef )[,4],col="Red")
points(mytimes,sizeHistory(params_data,mytimes, hdef, gammadef*10^10 )[,4],col="Green")

a <- 0.006
b <- 3.05
exp(log(MM[dim(MM)[1],4]/a)/b)

# use size history and a and b to get mizer simulated time vs length curve

herring_length <- function(params_data,mytimes, hdef, gammadef,w_infdef){
  return(sapply(sizeHistory(params_data,mytimes, hdef, gammadef,w_infdef)[,4],function(x) exp(log(x/a)/b)))
}

plot(mytimes,herring_length(params_data,mytimes, hdef, gammadef,w_infdef ))

### abundance based catches, with selectivity should be coded later

load("Herring.Rdata")
#out is name of this data = k, Linf, 
library(mvtnorm)

out

out2 <- out
out2[,2] <- out[,2]/10

# out2 is mikes samples from p(k|x), where k=(VB parameter, L_inf)
# we should fit this distribution using kernel density estimation 
# (approximating dsbn by sum of gaussians)

SIGMA<-cov(out2[,1:2])
MU<-colMeans(out2[,1:2])

dmvnorm(c(.58, 31),MU,SIGMA,log=T)


lengths <- herring_length(params_data,mytimes, hdef+5, gammadef,w_infdef )

# here data is set us like fishing catch data, but from mizer, with X=time, Y=weight
dats <- data.frame(X=mytimes, Y= lengths)
vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*X)))}

vbTyp

obs_Linf <- 29.3
obs_k <- 0.6
# Fit a von Bertalanffy growth curve to mizer data X, where dats is the mizer-generated 
# data
fitTyp<-nls(Y~vbTyp(X,Linf,k),data=dats,start=list(Linf=obs_Linf,k=obs_k))
coef(fitTyp)

MU

dmvnorm(coef(fitTyp),MU,SIGMA,log=T)

#############
#############

loglike <- function(h, gamma,w_inf){
  lengths <- herring_length(params_data,mytimes, h, gamma, w_inf )
  dats <- data.frame(X=mytimes, Y= lengths)
  vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*X)))}
  obs_Linf <- 29.3
  obs_k <- 0.6
  fitTyp<-nls(Y~vbTyp(X,Linf,k),data=dats,start=list(Linf=obs_Linf,k=obs_k))
  return(dmvnorm(coef(fitTyp),MU,SIGMA,log=T))
}

loglike(hdef,gammadef,w_infdef)

loglike(hdef,gammadef,10*w_infdef)


w2 <- w_infdef
LIKE <- 1:10

for (i in 1:10){
  w2[4] <- 0.01*i
  LIKE[i] <- loglike(hdef,gammadef,w2)
}
plot(LIKE)

############
############

### need to multiply up to go from mm to cm, in his data.

# This is the likelyhood p(k|theta)
exp(dmvnorm(coef(fitTyp),MU,SIGMA,log=T))

exp(dmvnorm(MU,MU,SIGMA,log=T))

MU

