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

data(inter)
data(NS_species_params_gears)
# Multiple species params
params <- MizerParams(NS_species_params_gears, inter=inter)
single_params <- MizerParams(NS_species_params_gears[1,])
# Multiple species, single effort sim
sim1 <- project(params, effort = 1, t_max = 10)
n_mult <- sim1@n[11,,]
n_single <- matrix(sim1@n[11,1,],nrow=1)
dimnames(n_single) <- list(sp = "Sprat", w = dimnames(sim1@n)$w)
npp <- sim1@n_pp[11,]
# Multiple species, mult effort sim
effort <- array(abs(rnorm(40)),dim=c(10,4))
single_effort <- array(abs(rnorm(10)),dim=c(10,1))

length(dim(getPhiPrey(params,n_mult,npp)))

length(dim(getPhiPrey(single_params,n_single,npp)))

length(npp)
length(params@w_full)

getPhiPrey(single_params,n_single,npp)

length(single_params)

dim(n_single)

dim(npp)

length(npp)


getPhiPrey(single_params,n_single,npp)

dim(matrix(npp,nrow=1))


n_single

expect_that(length(dim(getPhiPrey(params,n_mult,npp))), is_identical_to(length(dim(getPhiPrey(single_params,n_single,npp)))))


getPhiPrey(single_params,n_single,npp)

######### problem is discrepency between the following:

length(single_params@w_full)
length(npp)


