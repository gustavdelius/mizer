
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


params <- set_community_model(z0 = 0.1, f0 = 0.7, alpha = 0.2, recruitment = 4e7)

slotNames(params)

sim <- project(params, t_max=50, effort = 0)

plot(sim)

A <- sim@n_d

B <- sim@n_pp

C <- sim@n

dim(C)

plot(C[51,1,])

plot(A[51,])
