source('./R/MizerParams-class.r')
source('./R/MizerSim-class.R')
source('./R/project_methods.R')
source('./R/selectivity_funcs.R')
source('./R/summary_methods.R')
source('./R/wrapper_functions.R')
source('./R/project_methods.R')
source('./R/plots.R')
source('./R/project.R')
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)

    data(NS_species_params_gears)
    data(inter)
    params <- MizerParams(NS_species_params_gears, inter)
    sim <- project(params, effort=1, t_max=20, dt = 0.5, t_save = 0.5)
    time_range <- 15:20
    expect_that(length(dim(getM2(sim, time_range=time_range))), equals(3))
    time_range <- 20
  ##  expect_that(length(dim(getM2(sim, time_range=time_range))), equals(2))
    ### expect_that(getM2(sim, time_range=time_range), equals(getM2(sim@params, sim@n[as.character(time_range),,], sim@n_pp[as.character(time_range),])))
    
aq1 <- getM2(sim, time_range=time_range)
aq2 <- getM2(sim@params, sim@n[as.character(time_range),,], sim@n_pp[as.character(time_range),])

dim(aq1)

plot(aq1[1,])
lines(aq2[1,])

ttot <- 0
for (i in (1:dim(aq1)[1])){
ttot <- ttot + sum(aq1[i,]!=aq2[i,])    
}
ttot==0

dim(aq1)
dim(aq2)
