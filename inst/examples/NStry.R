library(mizer)
data(NS_species_params)
data(inter)
params <- MizerParams(params_data, interaction = inter)
#params <- MizerParams(params_data)
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
plot(sim)
plot(project(params))

library(progress)
library(mizer)

##########################

params_data_NS <- read.csv("./vignettes/NS_species_params.csv")
params <- MizerParams(params_data_NS, interaction = inter)
sim <- project(params, effort = 1, t_max = 50, dt = 0.1, t_save = 1)
plot(sim)
params@initial_n <- sim@n[dim(sim@n)[1],,]
params@initial_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]
params@species_params$erepro <- params@species_params$erepro*
    getRDD(params,params@initial_n,params@initial_n_pp)/getRDI(params,params@initial_n,params@initial_n_pp)
#params@srr <- function(rdi, species_params) {
#    return(rdi)
#}
params@species_params$r_max  <- Inf
sim <- project(params, effort = 1, t_max = 50, dt = 0.1, t_save = 1)
plot(sim)
##################


params_data_NS_no_rmax <- read.csv("./vignettes/NS_species_params.csv")

params_data_NS_no_rmax$r_max[] <- Inf
params <- MizerParams(params_data_NS_no_rmax, interaction = inter)
params@initial_n <- sim@n[dim(sim@n)[1],,]
params@initial_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]

#repro <- erepro *RDD/RDI

params@species_params$erepro <- params@species_params$erepro*
    getRDD(params,params@initial_n,params@initial_n_pp)/getRDI(params,params@initial_n,params@initial_n_pp)



sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
plot(sim)
params <- steady(params, effort = 1)

rdd_temp <- getRDD(params,params@initial_n,params@initial_n_pp)
gg <- getEGrowth(params,params@initial_n,params@initial_n_pp)

#params@species_params$w_min
#params@species_params$w_min[12]

# all w_min = w[1] in NS case, so below is simpler
rdd_desired <- gg[,1]*params@initial_n[,1]

# since r_max -> inf, rescaling the erepos is simpler
params@species_params$erepro <- params@species_params$erepro*rdd_desired/rdd_temp

sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
plot(sim)

params@species_params$erepro
# got the north sea model from the vignette, and turned off the r_max, and used steady so run 
# to a new steady state (when egg influx is held constant), and then the erepros are retunned 
# so that we find a new steady state where no stock recruitmemnt relationship is imposed.

######################################

params_data_NS <- read.csv("./vignettes/NS_species_params.csv")
params <- MizerParams(params_data_NS, interaction = inter)
sim <- project(params, effort = 1, t_max = 50, dt = 0.1, t_save = 1)
plot(sim)
params@initial_n <- sim@n[dim(sim@n)[1],,]
params@initial_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]
params@species_params$erepro <- params@species_params$erepro*
    getRDD(params,params@initial_n,params@initial_n_pp)/getRDI(params,params@initial_n,params@initial_n_pp)
#params@srr <- function(rdi, species_params) {
#    return(rdi)
#}
params@species_params$r_max  <- Inf
sim <- project(params, effort = 1, t_max = 50, dt = 0.1, t_save = 1)
plot(sim)
