library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
#library(mizer)
params_data <- read.csv("./vignettes/NS_species_params.csv")

inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")

####
params <- MizerParams(params_data, interaction = inter, no_w = 100,n=0.7)

ptm <- proc.time() 
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
plot(sim)

dim(large_predation(params))

plot(params@w,large_predation(params)[1,],log="x")

params@mu_ext <- large_predation(params)

plot(params@w,getM2(params,sim@n[dim(sim@n)[1],,],sim@n_pp[dim(sim@n_pp)[1],])[1,],log="x")
lines(params@w,large_predation(params)[1,])

params@mu_ext <- large_predation(params)
plot(params@w,getM2(params,sim@n[dim(sim@n)[1],,],sim@n_pp[dim(sim@n_pp)[1],])[1,],log="x")
