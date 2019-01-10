library(progress)
library(mizer)

##########################
# load data
data(inter)
params_data_NS <- read.csv("./vignettes/NS_species_params.csv")
params <- MizerParams(params_data_NS, interaction = inter,r_pp=5)
# default r_pp = 10
# run to a steady state of the original system
sim <- project(params, effort = 1, t_max = 50, dt = 0.1, t_save = 1)
plot(sim)
# use steady state as new initial condition
params@initial_n <- sim@n[dim(sim@n)[1],,]
params@initial_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]
# set erepro to compensate for turning off r_max
params@species_params$erepro <- params@species_params$erepro*
    getRDD(params,params@initial_n,params@initial_n_pp)/getRDI(params,params@initial_n,params@initial_n_pp)
#params@srr <- function(rdi, species_params) {
#    return(rdi)
#}
params@species_params$r_max  <- Inf

sim <- project(params, effort = 1, t_max = 50, dt = 0.1, t_save = 1)
plot(sim)

# this code is unstable when r_pp=5 on line 7, and goes into a small oscillation when r_pp=10 
# (default) and seems more stable when r_pp = 100. However, it is difficult to tell, because 
# in this code, the distance from the steady state one starts at also depends on r_pp. It would 
# be code to do a new version of this code, where for both r_pp values one sets up a steady state 
# (by what procedure - we maybe want Newton-Raphson to be exact), and then we can run our system at 
# a fixed pertubation away from the steady state, and gather evidence about whether it goes 
# oscillatory like a hopf bifurcation as r_pp changes.