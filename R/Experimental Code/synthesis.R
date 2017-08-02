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

params <- MizerParams(params_data, interaction = inter)

# run through all the fishing effort simulation once, with default 
# initial conditions, to make sensible initial conditions
# It would probably be smarter to cycle it round the 
# fishing effort years until it becomes periodic, we can 
# put this in later, but it looks like 1 cycle is enough

primer <- project(params, effort = t(Fmat), dt = 0.1, t_save =.1)




sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =.1,initial_n= primer@n[nrow(primer@n),,], initial_n_pp=primer@n_pp[nrow(primer@n_pp),])

#sim <- project(params, t_save =1)


# we should get the growth rates at the tsave points

plot(sim)

plot(sim@growth[,1,1])

# need to make initial condition properly
head(getEGrowth(params,sim@n[nrow(sim@n),,],sim@n_pp[nrow(sim@n_pp),])[1,])

head(getEGrowth(params,sim@n[nrow(sim@n)-1,,],sim@n_pp[nrow(sim@n_pp)-1,])[1,])
head(sim@growth[432,1,])

#####

# The sim@growth slot that we made mizer produce during sim
# is such that 
# sim@growth[t,,] == getEGrowth(params,sim@n[t,,],sim@n_pp[t,])
# holds for each t except the last, in that the 
# final growth rates are just copied from the penultimate ones in 
# this is fixed by the following command:
sim@growth[dim(sim@growth)[1],,] <- getEGrowth(params,sim@n[nrow(sim@n),,],sim@n_pp[nrow(sim@n_pp),])
plot(sim@growth[,1,1])


# Here ExtractGrowthRates(sim)[k,i,t] equals
# the growth rate g(W,T) where W=params@w[k] and 
# T = as.numeric(rownames(sim@n_pp))[t]
# is the t th time point at which mizer saves output


ExtractGrowthRates <- function(sim){
  HH <- sim@growth
  HH[dim(HH)[1],,] <- getEGrowth(params,sim@n[nrow(sim@n),,],sim@n_pp[nrow(sim@n_pp),])
  return(HH)
}

plot(ExtractGrowthRates(sim)[,1,1])

###################################
######################################

#param1 could be any params

#test<-project(param1,effort=0.5,t_max = 3, dt = 0.25, t_save = 1) ## run 3 years
#test1<-project(param1,effort=0.5,t_max = 1, dt = 0.25, t_save = 1, initial_n= test@n[nrow(test@n),,], initial_n_pp=test@n_pp[nrow(test@n_pp),]) # continue from previous run for 1 more year
#test2<-project(param1,effort=0.5,t_max = 4, dt = 0.25, t_save = 1) # run 4 years in one go.
#getYield(test1)-getYield(test2)[4,] # check to see that the results are the same

############

# put fixed point initial condition in
# check I understand what times I have the growth rates for
# extract growth rates
# backup on github
# investigate alpha, make markdown


#############
#############

#plank <- sim@n_pp
#as.numeric(rownames(plank))

eg <- ExtractGrowthRates(sim)
eg

TimesWeKnowG <- as.numeric(rownames(sim@n_pp))

# Here ExtractGrowthRates(sim)[k,i,t] equals
# the growth rate g(W,T) where W=params@w[k] and 
# T = as.numeric(rownames(sim@n_pp))[t]
# is the t th time point at which mizer saves output

##############

# 
w_pts <- params@w

#get_indexes <- function(W,T){
#  return(c(length(w_pts[w_pts <= W]),length(TimesWeKnowG[TimesWeKnowG <- T])))
#   
#  }


#inter_round_down <- function(W,T,eg){
#  return(eg[length(w_pts[w_pts <= W]), ,
#     length(TimesWeKnowG[TimesWeKnowG <- T])])
#}

w_pts[length(w_pts[w_pts <=1000])]

get_w_index <- function(W){
  return(length(w_pts[w_pts <=W]))
}

w_pts[get_w_index(1000)]

#length(TimesWeKnowG[TimesWeKnowG<=1999.7])

get_t_index <- function(T){
  return(length(TimesWeKnowG[TimesWeKnowG<=T]))
}
TimesWeKnowG[get_t_index(1999.734)]

eg <- ExtractGrowthRates(sim)

# Ginter approximates the growth rate @ weight W and time T, 
# by the growth rate known in eg at W* and T* where W* is the max val 
# less than or equal to W where the growth rate is known, and 
# T* is the max val 
# less than or equal to T where the growth rate is known

# Here Ginter(W,T,egg)[i] gives resulting growth rate for species i.

Ginter <- function(W,T,eg){
  return(eg[get_t_index(T),,get_w_index(W)])
}

Ginter(1000,1997.734,eg)

myw <- (1:400)
plot(myw,sapply(myw, function(x) Ginter(x,1997.734,eg)[4]))

params_data$w_inf[4]

Ginter(370,1997.734,eg)[4]

# Here ExtractGrowthRates(sim)[k,i,t] equals
# the growth rate g(W,T) where W=params@w[k] and 
# T = as.numeric(rownames(sim@n_pp))[t]
# is the t th time point at which mizer saves output

################

# get (age, length) points for a fish born at t0
# (I should work on all species in parralel, but lets start with 
# herring)

## prev approach
#weightsB <- matrix(0,nrow=length(mytimes),ncol=dim(testB@n)[2])
#for (i in (1:ns)){
#  gini <- approxfun(param1B@w, ggB[i,])
#  
#  myodefun <- function(t, state, parameters){
#    return(list(gini(state)))
#  }
#  weightsB[,i] <- ode(y = param1B@w[1], times = mytimes, func = myodefun, parms = 1)[,2]
#}
#
#lengthsB <- sapply(weightsB[,4],function(x) exp(log(x/a)/b))
##



myodefun <- function(t, state, parameters){
  return(list(Ginter(state,t,eg)))
}

ageWeight <- ode(y = rep(w_pts[1],12), times = TimesWeKnowG, func = myodefun, parms = 1)

head(ageWeight)

# time vs weight curve for a herring born at the start, in 1967

# we want a more general function where we can vary the birth time

plot(TimesWeKnowG,ageWeight[,4+1])

ageWeightGenBirth <- function(t0) {
  return(ode(y = rep(w_pts[1],12),
             times = TimesWeKnowG[TimesWeKnowG>=t0], func = myodefun, parms = 1))
  
}

# weight for a herring born in 1980

plot(ageWeightGenBirth(1980)[,1],ageWeightGenBirth(1980)[,4+1])

EarlyBirthData <- ageWeightGenBirth(TimesWeKnowG[1])
plot(EarlyBirthData[,1],EarlyBirthData[,4+1])

## connect this with prevnomod..

####################

###################################
################################### weld
###################################

a <- 0.006
b <- 3.05

out2 <- out
out2[,2] <- out[,2]/10
SIGMA<-cov(out2[,1:2])
MU<-colMeans(out2[,1:2])



params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
mytimes <- seq(0.5,10,0.1)
param1 <- MizerParams(params_data, interaction = inter, no_w = 100)
test<-project(param1,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)

loglikk <- function(loggam){
  #wB <- params_data[["w_inf"]]
  #wB[4] <- winfher
  wG <- log(param1@species_params[,15])
  # log(param1@species_params[4,15]) = -25
  
  wG[4] <- loggam
  
  params_dataB <- params_data
  params_dataB[["h"]] <- param1@species_params[,14]
  #params_dataB[["gamma"]] <- param1@species_params[,15]
  params_dataB[["gamma"]] <- exp(wG)
  #params_dataB[["w_inf"]] <- wB 
  params_dataB[["w_inf"]] <- params_data[["w_inf"]]
  param1B <- MizerParams(params_dataB, interaction = inter, no_w = 100)
  testB<-project(param1B,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)

  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
  params <- MizerParams(params_data, interaction = inter)
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
  
  
  Ginter <- function(W,T,eg){
    return(eg[get_t_index(T),,get_w_index(W)])
  }
  
  myodefun <- function(t, state, parameters){
    return(list(Ginter(state,t,eg)))
  }
  
  
  
  
  #EarlyBirthData <- ageWeightGenBirth(TimesWeKnowG[1])
  #plot(EarlyBirthData[,1],EarlyBirthData[,4+1])
  
  
  #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  
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
  #vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*X)))}
  
  W_egg <- param1@w[1]
  L_egg <- (W_egg/a)^(1/b)
  gett0 <- function(Linf,k){
    return(log(1-L_egg/Linf)/k)
  }
  
  #-k*(X-gett0(Linf,k))
  
  #vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*X)))}
  
  vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*(X-gett0(Linf,k)))))}
  
  
  
  # use MU to get observations
  
  #obs_Linf <- 29.3
  #obs_k <- 0.6
  
  # MU[[1]]
  
  obs_k <- MU[[1]]
  obs_Linf <- MU[[2]]
  
  
  # Fit a von Bertalanffy growth curve to mizer data X, where dats is the mizer-generated 
  # data
  fitTypB<-nls(Y~vbTyp(X,Linf,k),data=datsB,start=list(Linf=obs_Linf,k=obs_k))
  
  vbfitB <- coef(fitTypB)
  
  MU
  
  loglikeB <- dmvnorm(c(vbfitB[["k"]],vbfitB[["Linf"]]),MU,SIGMA,log=T)
  return(loglikeB)
  
}

loggamma <- seq(-26.3,-25,0.01)
ln_likelihood <- sapply(loggamma,loglikk)

plot(loggamma,ln_likelihood)

# most likely gamma
exp(loggamma[which(ln_likelihood==max(ln_likelihood))])

# actual gamma used by mizer
param1@species_params[4,15]



#plot(((-27):(-10)),sapply(seq(-27,-10,1),loglikk))

#loglikk(params_data[["w_inf"]][4])

#loglikk(10)
#loglikk(20)



#loglikevec <- 1:20

#for (i in (1:20)){
#  loglikevec[i] <- loglikk(20*i)
#}
#w_inf_herring <- 10*(1:20)
#plot(w_inf_herring,loglikevec)

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
W_egg <- param1@w[1]
L_egg <- (W_egg/a)^(1/b)
gett0 <- function(Linf,k){
  return(log(1-L_egg/Linf)/k)
}

#-k*(X-gett0(Linf,k))

#vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*X)))}

vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*(X-gett0(Linf,k)))))}

# should these ve set equal to MU ?

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
####################

param1@w[1]
a
b

params_data[["w_inf"]][4]
a*MU[2]^b

(params_data[["w_inf"]][4]/a)^(1/b)
MU[2]

W_egg <- param1@w[1]
L_egg <- (W_egg/a)^(1/b)
gett0

# are a and b for cm or mm, or does it matter

################ MCMC ####################

T <- 10000
theta <- -25
outputt <- 1:T
outputt[1] <- theta
for (t in (1:(T-1))){
  theta <- outputt[t]
  thetap <- rnorm(1,theta,0.001)
  r <- runif(1,0,1)
  outputt[t+1] <- theta
  if (r<exp(loglikk(thetap)-loglikk(theta))){
    outputt[t+1] <- thetap
  }
}
# mizer gamma
param1@species_params[4,15]

# MCMC gamma
exp(mean(outputt[200:1000]))

# Best fit gamma
exp(loggamma[which(ln_likelihood==max(ln_likelihood))])
