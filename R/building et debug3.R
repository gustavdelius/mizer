rm(list=ls())
.rs.restartR()
library(mizer)
library(plyr)
paramsConst2 <- set_trait_model(no_sp = 10, min_w_inf = 10,n=2/3,q=0.8,eta=0.25,
                                k0=10^(50),kappa=10^(50), alpha=0.17,
                                h=1614.363,f0=0.5,ks=0,z0=0,gamma=660.2633
)
sim <- project(paramsConst2, t_max=15, effort = 0)
plot(sim)
source("R/MizerParams-classETRNnew2.R")
# must be some new addition to MizerParams-classETRNnew2.R that causes this bug
# however it seems to get triggered only when project methods is called afterwards
# and then one tries to use getPhiPrey. It is strange since 
# MizerParams-classETRNnew2.R is identical to MizerParams-class appart
# from it has a new clone slot psi2
source("R/project_methods.R")
n <- sim@n[dim(sim@n)[1],,]
n_pp <- sim@n_pp[dim(sim@n_pp)[1],]
gtphi <- getPhiPrey(sim@params, n=n, n_pp=n_pp)

