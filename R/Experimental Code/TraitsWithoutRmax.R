library(mizer)
library(plyr)

epsi <- 0.25

kappaR2 <- 0.005
lambda2 <- 2+0.9-(2/3)
source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
source("./R/Experimental Code/projectmodPREYSWITCH2.R")
paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5,k0=50*10^(50))
simConst <- project(paramsConst, t_max=50, effort = 0)
plot(simConst)


##################################
epsi <- .1
kappaR2 <- 10^(11)
lambda2 <- 2+0.9-(2/3)
source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
source("./R/Experimental Code/projectmodPREYSWITCH2.R")

no_sp <- 12
params_data <- read.csv("./vignettes/NS_species_params.csv")
reduced_data <- params_data[1:no_sp,]
reduced_data$r_max <- rep(10^(70),length(reduced_data$r_max))
params <- MizerParams(reduced_data)
sim <- project(params, effort = 1, t_max = 200, dt = 0.01, t_save = 1)
ini_n <- sim@n[dim(sim@n)[1],,]
ini_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]

sim2 <- project(params, effort = 1, t_max = 50, dt = 0.01, t_save = 1,initial_n=ini_n,initial_n_pp=ini_n_pp)
plot(sim2)

sim3 <- project(params, effort = 1, t_max = 50, dt = 0.05, t_save = 1,initial_n=ini_n,initial_n_pp=ini_n_pp)
plot(sim3)

sim4 <- project(params, effort = 1, t_max = 50, dt = 0.001, t_save = 1,initial_n=ini_n,initial_n_pp=ini_n_pp)
plot(sim4)

#################################

no_sp <- 11
params_data <- read.csv("./vignettes/NS_species_params.csv")
reduced_data <- params_data[1:no_sp,]
reduced_data$r_max <- rep(10^(70),length(reduced_data$r_max))
params <- MizerParams(reduced_data)
sim <- project(params, effort = 1, t_max = 200, dt = 0.1, t_save = 1)
ini_n <- sim@n[dim(sim@n)[1],,]
ini_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]

sim2 <- project(params, effort = 1, t_max = 50, dt = 0.01, t_save = 1,initial_n=ini_n,initial_n_pp=ini_n_pp)
plot(sim2)

# # # 
no_sp <- 11
params_data <- read.csv("./vignettes/NS_species_params.csv")
reduced_data <- params_data[-11,]
reduced_data$r_max <- rep(10^(70),length(reduced_data$r_max))
params <- MizerParams(reduced_data)
sim <- project(params, effort = 1, t_max = 200, dt = 0.1, t_save = 1)
ini_n <- sim@n[dim(sim@n)[1],,]
ini_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]

sim2 <- project(params, effort = 1, t_max = 50, dt = 0.01, t_save = 1,initial_n=ini_n,initial_n_pp=ini_n_pp)
plot(sim2)

##############################
##############################

library(mizer)
library(plyr)

epsi <- 0.1

kappaR2 <- 0.005
lambda2 <- 2+0.9-(2/3)
source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
source("./R/Experimental Code/projectmodPREYSWITCH2.R")

params_data <- read.csv("./vignettes/NS_species_params.csv")
reduced_data <- params_data
reduced_data$r_max <- rep(10^(70),length(reduced_data$r_max))
params <- MizerParams(reduced_data)
sim <- project(params, effort = 1, t_max = 200, dt = 0.1, t_save = 1)
ini_n <- sim@n[dim(sim@n)[1],,]
ini_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]

simA <- project(params, effort = 1, t_max = 50, dt = 0.01, t_save = 1,initial_n=ini_n,initial_n_pp=ini_n_pp)
plot(simA)
