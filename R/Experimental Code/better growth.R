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