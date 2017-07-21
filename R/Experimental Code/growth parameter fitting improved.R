library(devtools)
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library(mizer)
library(deSolve)
load("Herring.Rdata")
#out is name of this data = k, Linf, 
library(mvtnorm)

a <- 0.006
b <- 3.05

out2 <- out
out2[,2] <- out[,2]/10
SIGMA<-cov(out2[,1:2])
MU<-colMeans(out2[,1:2])



params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
mytimes <- seq(0,100,0.1)
param1 <- MizerParams(params_data, interaction = inter, no_w = 100)
test<-project(param1,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)

loglikk <- function(winfher){
  wB <- params_data[["w_inf"]]
  wB[4] <- winfher
  
  params_dataB <- params_data
  params_dataB[["h"]] <- param1@species_params[,14]
  params_dataB[["gamma"]] <- param1@species_params[,15]
  params_dataB[["w_inf"]] <- wB
  param1B <- MizerParams(params_dataB, interaction = inter, no_w = 100)
  testB<-project(param1B,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)
  
  ns <- dim(testB@n)[2]
  nbb <- testB@n[dim(testB@n)[1],,]
  ggB <- getEGrowth(param1B, testB@n[dim(testB@n)[1],,], testB@n_pp[dim(testB@n_pp)[1],])
  sols <- as.list(1:ns)
  weightsB <- matrix(0,nrow=length(mytimes),ncol=dim(testB@n)[2])
  for (i in (1:ns)){
    gini <- approxfun(param1B@w, ggB[i,])
    
    myodefun <- function(t, state, parameters){
      return(list(gini(state)))
    }
    weightsB[,i] <- ode(y = param1B@w[1], times = mytimes, func = myodefun, parms = 1)[,2]
  }
  
  lengthsB <- sapply(weightsB[,4],function(x) exp(log(x/a)/b))
  
  ############
  
  datsB <- data.frame(X=mytimes, Y= lengthsB)
  vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*X)))}
  obs_Linf <- 29.3
  obs_k <- 0.6
  # Fit a von Bertalanffy growth curve to mizer data X, where dats is the mizer-generated 
  # data
  fitTypB<-nls(Y~vbTyp(X,Linf,k),data=datsB,start=list(Linf=obs_Linf,k=obs_k))
  
  vbfitB <- coef(fitTypB)
  
  MU
  
  loglikeB <- dmvnorm(c(vbfitB[["k"]],vbfitB[["Linf"]]),MU,SIGMA,log=T)
  return(loglikeB)

}


loglikk(params_data[["w_inf"]][4])

loglikk(10)
loglikk(20)



loglikevec <- 1:20

for (i in (1:20)){
  loglikevec[i] <- loglikk(20*i)
}
w_inf_herring <- 10*(1:20)
plot(10*(1:20),loglikevec)

############################### specific run ###############
# generate likelihood for default settings, first by making params and sim object, test.
params_dataA <- params_data
params_dataA[["h"]] <- param1@species_params[,14]
params_dataA[["gamma"]] <- param1@species_params[,15]
params_dataA[["w_inf"]] <- params_data[["w_inf"]]
param1A <- MizerParams(params_dataA, interaction = inter, no_w = 100)
testA<-project(param1A,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)

ns <- dim(testA@n)[2]
naa <- testA@n[dim(testA@n)[1],,]
gg <- getEGrowth(param1A, testA@n[dim(testA@n)[1],,], testA@n_pp[dim(testA@n_pp)[1],])
sols <- as.list(1:ns)
weightsA <- matrix(0,nrow=length(mytimes),ncol=dim(testA@n)[2])
for (i in (1:ns)){
  gini <- approxfun(param1A@w, gg[i,])
  
  myodefun <- function(t, state, parameters){
    return(list(gini(state)))
  }
  weightsA[,i] <- ode(y = param1@w[1], times = mytimes, func = myodefun, parms = 1)[,2]
}

lengthsA <- sapply(weightsA[,4],function(x) exp(log(x/a)/b))

############

datsA <- data.frame(X=mytimes, Y= lengthsA)
vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*X)))}
obs_Linf <- 29.3
obs_k <- 0.6
# Fit a von Bertalanffy growth curve to mizer data X, where dats is the mizer-generated 
# data
fitTypA<-nls(Y~vbTyp(X,Linf,k),data=datsA,start=list(Linf=obs_Linf,k=obs_k))

vbfitA <- coef(fitTypA)

loglikeA <- dmvnorm(c(vbfitA[["k"]],vbfitA[["Linf"]]),MU,SIGMA,log=T)


#####################################################


# run ps debug first to arm the rest

vbfitD <- vbfitA

