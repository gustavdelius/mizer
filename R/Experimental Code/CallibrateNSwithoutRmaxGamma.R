library(mizer)
library(plyr)

params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
#load("Fmat.RData")
#Fmat <- Fmat[,11:(dim(Fmat)[2]-1)]
load("Landings.RData")
#params_data$catchability <- as.numeric(colMeans((f_history)[19:29,]))
landings <- t(landings)
landings[is.na(landings)] <- 0
landings <- landings[18:(dim(landings)[1]-1),]
params_data$sel_func <- "sigmoid_length"
params_data$l25 <- c(7.6, 9.8, 8.7, 10.1, 11.5, 19.8, 16.4, 19.8, 11.5,
                     19.1, 13.2, 35.3)
params_data$l50 <- c(8.1, 11.8, 12.2, 20.8, 17.0, 29.0, 25.8, 29.0, 17.0,
                     24.3, 22.9, 43.6)
params_data$a <- c(0.007, 0.001, 0.009, 0.002, 0.010, 0.006, 0.008, 0.004,
                   0.007, 0.005, 0.005, 0.007)
params_data$b <- c(3.014, 3.320, 2.941, 3.429, 2.986, 3.080, 3.019, 3.198,
                   3.101, 3.160, 3.173, 3.075)
f_history <- read.csv("./vignettes/NS_f_history.csv", row.names=1)
f_history <- as(f_history, "matrix")
#params_data$catchability <- as.numeric(f_history["1990",])
params_data$catchability <- as.numeric(colMeans((f_history)[19:29,]))
params <- MizerParams(params_data, inter, kappa = 9.27e10)
#relative_effort <- sweep(f_history,2,f_history["1990",],"/")
relative_effort <- sweep(f_history,2,colMeans((f_history)[19:29,]),"/")
relative_effort[as.character(1988:1992),]
initial_effort <- matrix(relative_effort[1,],byrow=TRUE, nrow=100,
                         ncol=ncol(relative_effort), dimnames = list(1867:1966))
##relative_effort <- rbind(initial_effort,relative_effort)

load("for_richard.RData")
# find entry with Haddock in it, and write in the Rmax manually, while checking against for_richard.Rdata

### more accurate part specification

getnsindex <- function(s){
  matcher <- match(params_data$species,s)==1
  matcher[is.na(matcher)] <- FALSE
  return((1:length(matcher))[matcher])
}
getnsindex("Haddock")
getnsindex("Cod")

capacity <- exp(25.210)
rmax <- exp(c(26.7, 26, 30.67,26.56,23.1,26.03, 22.95, 25.38,30.56, 28.375, 22.77,26.92))
rmax[getnsindex("Sprat")] <- exp(26.659)
rmax[getnsindex("Sandeel")] <- exp(26.008)
rmax[getnsindex("Norway pout")] <- exp(30.684)
rmax[getnsindex("Dab")] <- exp(23.108)
rmax[getnsindex("Herring")] <- exp(26.556)
rmax[getnsindex("Sole")] <- exp(22.948)
rmax[getnsindex("Whiting")] <- exp(26.034)
rmax[getnsindex("Plaice")] <- exp(30.562)
rmax[getnsindex("Haddock")] <- exp(28.375)
rmax[getnsindex("Saithe")] <- exp(26.920)
rmax[getnsindex("Cod")] <- exp(22.767)
dd <- params_data
dd$r_max <- rmax[1:(length(dd$r_max))]
dd$r_max <- rep(10^49,12)
dd$erepro <- rep(0.9,12)


params <- MizerParams(dd, interaction = inter, kappa=capacity)
#params <- MizerParams(params_data, interaction = inter, kappa=9.27e10)

epsi <- 0.5
kappaR2 <- capacity
lambda2 <- 2+0.8-(2/3)
source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
source("./R/Experimental Code/projectmodPREYSWITCH2.R")


simini <- project(params, effort = relative_effort, dt = 0.1, t_save =1)
sim <- project(params, effort = relative_effort, dt = 0.1, t_save =1, initial_n=simini@n[dim(simini@n)[1],,],initial_n_pp=simini@n_pp[dim(simini@n_pp)[1],])

vv <- log(getYield(sim)*10^(-6))


j <- 3
params_data$species[j]
plot(log(landings[,j]))
lines(vv[(dim(vv)[1]+1-dim(landings)[1]):dim(vv)[1],j])


vv <- log((getYield(sim)+10^(-10))*10^(-6))
# sum of squares
sum((vv[(dim(vv)[1]+1-dim(landings)[1]):dim(vv)[1],]-log(10^(-10)+landings[,]))^2)

####################### make function to minimize #############

#from 1:13 from params@species_params$gamma

#mypar <- c(log(88830864522),rep(0.9,12))
mypar <- c(log(88830864522),log(c(3.591300e-11, 1.692620e-11, 9.573345e-11, 1.264439e-11, 5.229280e-11, 8.319884e-11,
                                  3.506987e-11, 3.319288e-11, 3.204404e-11, 4.544627e-11, 1.798108e-10, 1.385290e-10)))

runnit <- function(par=mypar){
  dd <- params_data
  dd$r_max <- rep(10^49,12)
  #dd$erepro <- par[2:13]
  dd$gamma <- exp(par[2:13])
  
  params <- MizerParams(dd, interaction = inter, kappa=exp(par[1]))
  #params <- MizerParams(params_data, interaction = inter, kappa=9.27e10)
  
  epsi <- 0.5
  kappaR2 <- exp(par[1])
  lambda2 <- 2+0.8-(2/3)
  source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
  source("./R/Experimental Code/projectmodPREYSWITCH2.R")
  
  
  simini <- project(params, effort = relative_effort, dt = 0.1, t_save =1)
  sim <- project(params, effort = relative_effort, dt = 0.1, t_save =1, initial_n=simini@n[dim(simini@n)[1],,],initial_n_pp=simini@n_pp[dim(simini@n_pp)[1],])
  return(sim)
}

simm <- runnit()
vv <- log((getYield(simm)+10^(-10))*10^(-6))
j <- 3
params_data$species[j]
plot(log(landings[,j]))
lines(vv[(dim(vv)[1]+1-dim(landings)[1]):dim(vv)[1],j])





minme <- function(par=mypar){
  sim <- runnit(par)
  vv <- log((getYield(sim)+10^(-10))*10^(-6))
  # sum of squares
  #addon <- 10^49
  #if (all(par[2:13]>=rep(0,12))) {
  #  if (all(par[2:13]<=rep(1,12))) {
  #    addon <- 0
  #  }
  #}
  #return(addon+sum((vv[(dim(vv)[1]+1-dim(landings)[1]):dim(vv)[1],]-log(10^(-10)+landings[,]))^2))
  gb <- getBiomass(sim)
  extinctionPunishment <- 10^(49)
  if (min(gb[dim(gb)[1],])>0) {
    if ((max(gb[dim(gb)[1],])/min(gb[dim(gb)[1],]))<=10^2) {
      extinctionPunishment <- 0
    }
  }
  return(extinctionPunishment+sum((vv[(dim(vv)[1]+1-dim(landings)[1]):dim(vv)[1],]-log(10^(-10)+landings[,]))^2))
}
minme()

opout <- c(25.01638, -23.74906, -25.36743, -22.75435, -25.31926, -23.82399, -23.69747,
   -24.45443, -23.59808, -22.60934, -23.74995, -22.87058, -22.00435)
  

op <- optim(par=opout, fn=minme, method = "SANN", control = list(maxit = 1000))
mypar
op$par

#[1]  25.32816 -24.21185 -25.00448 -22.68365 -25.08402 -21.67672 -24.53056
#[8] -24.51131 -24.55053 -24.01618 -25.77193 -21.16473 -21.80833

#newer
#[1]  25.01638 -23.74906 -25.36743 -22.75435 -25.31926 -23.82399 -23.69747
#[8] -24.45443 -23.59808 -22.60934 -23.74995 -22.87058 -22.00435

simm <- runnit(op$par)
#simm <- runnit(mypar)

vv <- log((getYield(simm)+10^(-10))*10^(-6))
for (j in (1:12)){
  params_data$species[j]
  plot(log(landings[,j]))
  lines(vv[(dim(vv)[1]+1-dim(landings)[1]):dim(vv)[1],j])
}
  
# add condition to promote coexistence
# 
plot(simm)
gb <- getBiomass(simm)
extinctionPunishment <- 0
if (max(gb[dim(gb)[1],])/min(gb[dim(gb)[1],])>10^3){
  extinctionPunishment <- 10^(49)
}

