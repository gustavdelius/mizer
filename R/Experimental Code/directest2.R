set.seed(5)
library(devtools)
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library(mizer)
library(FME)
library(reshape2)
library(deSolve)
library("plot3D")
library(rgl)
library("plot3Drgl")
library(optimx)
library(mvtnorm)


# Run appropriate core mizer files again, that have been 
# rewritten to include sim@growth slots
source("R/project.R")
source("R/MizerSim-class.R")

# Fmat is fishing effort data, we use it from 1967 until 2010
# for our mizer runs
load("Fmat.RData")

# Load basic North sea parameters
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")

### extract real herring data 
# save and load to external file for speed and low mem

# data originates from the large file
# north_sea_data.csv
# at
# https://www.dropbox.com/s/078bwgewimwvnbs/north_sea_data.csv?dl=0

#data1 <- read.csv("north_sea_data.csv")
#data2 <- data1[which(data1$latin_name=="Clupea harengus"),]
#data3 <- data2[data2$Year==2010,]
#write.csv(data3, "herringage2010.csv")

# we just load up the subset of the data required here
Hdata <- read.csv("herringage2010.csv")

#############

# Set Fmat be be over the proper years
Fmat <- Fmat[,11:(dim(Fmat)[2]-1)]

# a and b are such that W=a*L^b for Herring (species # 4)
# TO DO: get these a's and b's for other NS species, and generalize 
a <- 0.006
b <- 3.05

Fmat[1,]

dim(Fmat)

# fishing effort data for 2010
effort2010 <- Fmat[,44]


params_init <- MizerParams(params_data, interaction = inter, no_w = 100)
primer <- project(params_init, effort = effort2010, t_max=100 , dt = 0.1, t_save =.1)


# gamma is the volumetric search rate constant (for Herring)
# set the log gamma value to work with.
loggam <- -24
# mizer-encoded `empirical` value is -25.21224
# $$$$$$$$ # loglikk <- function(loggam){

# Initialize mizer with altered gamma
param1 <- MizerParams(params_data, interaction = inter, no_w = 100)
wG <- log(param1@species_params[,15])
wG[4] <- loggam
params_dataB <- params_data
## params_dataB[["h"]] <- param1@species_params[,14]
params_dataB[["gamma"]] <- exp(wG)
params_dataB[["w_inf"]] <- params_data[["w_inf"]]
params <- MizerParams(params_dataB, interaction = inter, no_w = 100)

###primer <- project(params, effort = t(Fmat), dt = 0.1, t_save =.1)
###sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =.1,initial_n= primer@n[nrow(primer@n),,], initial_n_pp=primer@n_pp[nrow(primer@n_pp),])


# The fishing efforts are held fixed at those of 2010, 
# we hope the system can settle to fixture in 20 years for 
# different values of gamma
sim <- project(params, effort = effort2010, t_max=20 , dt = 0.1, t_save =.1,initial_n= primer@n[nrow(primer@n),,], initial_n_pp=primer@n_pp[nrow(primer@n_pp),])

plot(sim)

GR <- getEGrowth(params, sim@n[nrow(sim@n),,], sim@n_pp[nrow(sim@n_pp),])

dim(GR)

the_times <- seq(0,40,0.5)

plot(params@w,GR[4,],log="xy")

GRint <- approxfun(params@w,GR[4,])

some_weights <- (1:1000)
plot(some_weights, sapply(some_weights,GRint))

plot(params@w[params@w<1000],GR[4,params@w<1000])

#########################

# maximum age empirically recorded
ma <- max(Hdata$Age)
the_times <- seq(0,1+ma,0.5)

  myodefun <- function(t, state, parameters){
    return(list(GRint(state)))
  }
  weightsA <- ode(y = param1@w[1], times = the_times, func = myodefun, parms = 1)[,2]
plot(the_times,weightsA )

weightsA[length(weightsA)]

lengthsA <- sapply(weightsA,function(x) exp(log(x/a)/b))

plot(the_times,lengthsA )

params_data[["w_inf"]][4]

###############


# for now we just add 0.5 to empirically recorded age


emp_ages <- Hdata$Age
emp_lengths <- Hdata$Lngt
points(emp_ages,emp_lengths)
dim(data3)

plot(emp_ages,emp_lengths)
lines(the_times,lengthsA )

#ma <- max(data3$Age)
#AA <- (0:ma)
#AL <- (0:ma)
#for (i in (0:ma)){
#  AL[i+1] <- mean(data3[data3$Age==i,]$Lngt)
#}
#plot(AA,AL)


## compute likelihood
# need to estimate sigma noise

# have fractional part of age as described by 
# distribution in terms of extra parameters

# how to test log normallcy of data ?

emp_ages[1]
emp_lengths[1]

mizer_interp <- approxfun(the_times,lengthsA)

mizer_interp(emp_ages[1]+0.5)

#ss <- mizer_interp(emp_ages[1]+0.5)-emp_lengths[1]

#ss<- sum((sapply(emp_ages+0.5,mizer_interp) - emp_lengths)^2)
#ss

ss <- sum((emp_lengths-sapply(emp_ages+0.5,mizer_interp))^2)

# We should convert to logs in order to properly evaluate likelihood
# and maybe we should use
#ss <- sum((log(emp_lengths)-log(sapply(emp_ages+0.5,mizer_interp)))^2)



#mizer_interp(1.5)
#a1a <- sapply(emp_ages+0.5,mizer_interp) 
#a1b <- emp_lengths
#sum((a1b-a1a)^2)

#emp_ages[2109]
