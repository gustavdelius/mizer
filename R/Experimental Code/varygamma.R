
set.seed(5)

library(devtools)
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library(mizer)

library(mizer)
library(FME)
library(ggplot2)
library(reshape2)
library(deSolve)
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library("plot3D")
library(rgl)
library("plot3Drgl")
library(optimx)
library(deSolve)

set.seed(123)

source("R/project.R")
source("R/MizerSim-class.R")
load("Landings.RData")
load("Fmat.RData")

params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")

Fmat <- Fmat[,11:(dim(Fmat)[2]-1)]

#### rest of primining
load("Herring.Rdata")
#out is name of this data = k, Linf, 
library(mvtnorm)

set.seed(5)

a <- 0.006
b <- 3.05
load("Herring.Rdata")

out2 <- out
out2[,2] <- out[,2]/10
SIGMA<-cov(out2[,1:2])
MU<-colMeans(out2[,1:2])

####

# $$$$$ $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#uncomment out fn part later

#loggam <- -25
loggam <- -24
loglikk <- function(loggam){
#wB <- params_data[["w_inf"]]
#wB[4] <- winfher

#param1 <- params

param1 <- MizerParams(params_data, interaction = inter, no_w = 100)

wG <- log(param1@species_params[,15])
# log(param1@species_params[4,15]) = -25

wG[4] <- loggam

params_dataB <- params_data
# params_dataB[["h"]] <- param1@species_params[,14]
params_dataB[["gamma"]] <- exp(wG)
params_dataB[["w_inf"]] <- params_data[["w_inf"]]

params <- MizerParams(params_dataB, interaction = inter, no_w = 100)
primer <- project(params, effort = t(Fmat), dt = 0.1, t_save =.1)
sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =.1,initial_n= primer@n[nrow(primer@n),,], initial_n_pp=primer@n_pp[nrow(primer@n_pp),])


ExtractGrowthRates <- function(sim){
  HH <- sim@growth
  HH[dim(HH)[1],,] <- getEGrowth(params,sim@n[nrow(sim@n),,],sim@n_pp[nrow(sim@n_pp),])
  return(HH)
}
TimesWeKnowG <- as.numeric(rownames(sim@n_pp))
w_pts <- params@w
get_w_index <- function(W){
  return(length(w_pts[w_pts <=W]))
}
get_t_index <- function(T){
  return(length(TimesWeKnowG[TimesWeKnowG<=T]))
}
eg <- ExtractGrowthRates(sim)
#!!!!!!!!!!!!!!!!!!!
#Ginter <- function(W,T,eg){
#  return(eg[get_t_index(T),,get_w_index(W)])
#}
#myodefun <- function(t, state, parameters){
#  return(list(Ginter(state,t,eg)))
#}
#ageWeightGenBirth <- function(t0) {
#  return(ode(y = rep(w_pts[1],12),
#             times = TimesWeKnowG[TimesWeKnowG>=t0], func = myodefun, parms = 1))
#  
#}
#EarlyBirthData <- ageWeightGenBirth(TimesWeKnowG[1])
#weightsB <- EarlyBirthData[,2:13]
#!!!!!!!!!!!!!!
Ginter <- function(W,T,eg){
  return(eg[get_t_index(T),,get_w_index(W)])
}
EarlyBirthData <- matrix(0,nrow=length(TimesWeKnowG),ncol=12)
for (i in (1:12)){
  myodefun <- function(t, state, parameters){
    return(list(Ginter(state,t,eg)[i]))
  }
  
  ageWeightGenBirth <- function(t0) {
    # for some reason lsoda gives errors, so we use euler for now
    # return(ode(y = w_pts[1],
    #           times = TimesWeKnowG[TimesWeKnowG>=t0], func = myodefun, parms = 1))
    return(ode(y = w_pts[1],
               times = TimesWeKnowG[TimesWeKnowG>=t0], func = myodefun, parms = 1, method = "euler"))    
  }
  EarlyBirthData[,i] <- ageWeightGenBirth(TimesWeKnowG[1])[,2]
}
weightsB <- EarlyBirthData
#!!!!!!!!!!!!!!!!!
lengthsB <- sapply(weightsB[,4],function(x) exp(log(x/a)/b))
# # # # # # mini weld
# should edit list before use, to avoid using data 
#with extreemal ages
#datsB <- data.frame(X=TimesWeKnowG, Y= lengthsB)
#datsB <- data.frame(X=TimesWeKnowG[lengthsB>1], Y= lengthsB[lengthsB>1])
datsB <- data.frame(X=TimesWeKnowG-TimesWeKnowG[1], Y= lengthsB)

param1 <- params

W_egg <- param1@w[1]
L_egg <- (W_egg/a)^(1/b)
gett0 <- function(Linf,k){
  return(log(1-L_egg/Linf)/k)
}
vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*(X-gett0(Linf,k)))))}
obs_k <- MU[[1]]
obs_Linf <- MU[[2]]
fitTypB<-nls(Y~vbTyp(X,Linf,k),data=datsB,start=list(Linf=obs_Linf,k=obs_k))
vbfitB <- coef(fitTypB)
loglikeB <- dmvnorm(c(vbfitB[["k"]],vbfitB[["Linf"]]),MU,SIGMA,log=T)
loglikeB
return(loglikeB)}


loglikk(-24)

ln_gam <- seq( -26, -24, .5)
ln_like <- sapply(ln_gam, loglikk)
plot(ln_gam,ln_like)

loglikk(-28)

loglikk(-26)

ln_gam <- seq( -26.3, -24, .1)
ln_like <- sapply(ln_gam, loglikk)
plot(ln_gam,ln_like)


# some type of error from line 149 of debug, about
#   return(log(1-L_egg/Linf)/k)
