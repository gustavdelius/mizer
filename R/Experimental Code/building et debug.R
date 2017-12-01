rm(list=ls())
.rs.restartR()


library(mizer)
library(plyr)

mmu_et <- readRDS(file = "my_data.rds")

epsi <- 0.01
et_on <- 1
source("R/Experimental Code/project_methods_much_switch2.R")

source("R/MizerParams-classETRN.R")
source("R/wrapper_functionsETRN.R")
# for some reason the code breaks when project_methods is executed after mizer_prms
paramsConst2 <- set_trait_model(no_sp = 10, min_w_inf = 10,n=2/3,q=0.8,eta=0.25,
                                k0=10^(50),kappa=10^(50), alpha=0.17,
                                h=1614.363,f0=0.5,ks=0,z0=0,gamma=660.2633
)

sim <- project(paramsConst2, t_max=150, effort = 0)
plot(sim)
