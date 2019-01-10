library(progress)
library(mizer)

##########################
# load data
data(inter)
params_data_NS <- read.csv("./vignettes/NS_species_params.csv")
params <- MizerParams(params_data_NS, interaction = inter,r_pp=10)
# default r_pp = 10
# run to a steady state of the original system
sim <- project(params, effort = 1, t_max = 500, dt = 0.1, t_save = 1)
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

#params@rr_pp[] <- 

#res@rr_pp[] <- r_pp * res@w_full^(n - 1)
params@species_params$r_max  <- Inf

########## run simulation without using r_max (note we get a better steady state if we repeat 
# the above code with r_pp = 15)

sim <- project(params, effort = 1, t_max = 50, dt = 0.1, t_save = 1)
plot(sim)

######### run experiment again with higher r_pp

rr_pp2 <- 15
params2 <- params
params2@rr_pp[] <- rr_pp2*params2@w_full^(params2@n-1)
sim2 <- project(params2, effort = 1, t_max = 50, dt = 0.1, t_save = 1)
plot(sim2)

######### run experiment again with lower r_pp

rr_pp3 <- 5
params3 <- params
params3@rr_pp[] <- rr_pp3*params3@w_full^(params3@n-1)
sim3 <- project(params3, effort = 1, t_max = 50, dt = 0.1, t_save = 1)
plot(sim3)

# It seems decreasing r_pp to 5 hurts co-existence more than increasing r_pp from 10 to 15 does.

# wrote code that sets up NS steady state withour rmax, and investigates how 
# altering the value of r_pp effects stability in a couple of instances.