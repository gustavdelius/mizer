library(mizer)
library(plyr)

epsi <- 0.05

kappaR2 <- 0.005
lambda2 <- 2+0.9-(2/3)
#source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
#source("./R/Experimental Code/projectmodPREYSWITCH2.R")
source("./R/Experimental Code/project_methods_medium_change.R")

paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5,k0=50*10^(50),lambda = ((3+0.9-(2/3))/(epsi+1))-1)
simConst <- project(paramsConst, t_max=250, effort = 0)
plot(simConst)


##################################
epsi <- .1
kappaR2 <- 10^(11)
lambda2 <- 2+0.9-(2/3)
#source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
#source("./R/Experimental Code/projectmodPREYSWITCH2.R")
source("./R/Experimental Code/project_methods_medium_change.R")
no_sp <- 12
params_data <- read.csv("./vignettes/NS_species_params.csv")
reduced_data <- params_data[1:no_sp,]
reduced_data$r_max <- rep(10^(70),length(reduced_data$r_max))
params <- MizerParams(reduced_data)
sim <- project(params, effort = 1, t_max = 100, dt = 0.1, t_save = 1)
plot(sim)
