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
inter[,"Dab"] <- 0 # Dab not eaten by anything
params <- MizerParams(NS_species_params_gears, inter)
m2 <- getM2(params,get_initial_n(params),params@cc_pp)
##expect_that(all(m2["Dab",] == 0), is_true())

m2["Dab",]

m2
dim(m2)

#### did we lose dimension names ?

dim(NS_species_params_gears)

NS_species_params_gears

get_initial_n(params)

params@n

rownames(params@interaction)
