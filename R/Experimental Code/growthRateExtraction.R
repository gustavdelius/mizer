library(devtools)
#install_github("gustavdelius/mizer", ref = "master")
#install_github("drfinlayscott/mizer", ref = "master")


library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library(mizer)
#setwd("C:/Users/richa/Dropbox/Job/YorkJob/JulyMizerDB/mizer-fftclean2")


source("R/project.R")
source("R/MizerSim-class.R")


params_data <- read.csv("./vignettes/NS_species_params.csv")

inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")

####
params <- MizerParams(params_data, interaction = inter, no_w = 100)

ptm <- proc.time() 
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
Tsim <- (proc.time() - ptm)[[3]]

plot(sim)

param1 <- params

test<-project(param1,effort=0.5,t_max = 3, dt = 0.25, t_save = 1) ## run 3 years

dim(test@growth)

plot(test@growth[,1,1])
plot(test@growth[1,1,])
