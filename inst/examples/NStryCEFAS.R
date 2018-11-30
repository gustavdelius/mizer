
######### initiate old system

library(mizer)
data(NS_species_params)
data(inter)

params_data_NS <- read.csv("./vignettes/NS_species_params.csv")
params <- MizerParams(params_data_NS, interaction = inter)
simR <- project(params, effort = 1, t_max = 500, dt = 0.1, t_save = 1)
plot(simR)

########## run old system, with new effort levels


simRe <- project(params, effort = 1.1, t_max = 10, dt = 0.1, t_save = 1,
                 initial_n = simR@n[dim(simR@n)[1],,], initial_n_pp = 
                     simR@n_pp[dim(simR@n_pp)[1],])
plot(simRe)

########### initialize new system, without rmax

params@initial_n <- simR@n[dim(simR@n)[1],,]
params@initial_n_pp <- simR@n_pp[dim(simR@n_pp)[1],]
params@species_params$erepro <- params@species_params$erepro*
    getRDD(params,params@initial_n,params@initial_n_pp)/getRDI(params,params@initial_n,params@initial_n_pp)

params@species_params$r_max  <- Inf
simF <- project(params, effort = 1, t_max = 50, dt = 0.05, t_save = 1)
plot(simF)

############# run new system, with different fishing effort

simFe <- project(params, effort = 1.1, t_max = 50, dt = 0.1, t_save = 1,
                 initial_n = simF@n[dim(simF@n)[1],,], initial_n_pp = 
                     simF@n_pp[dim(simF@n_pp)[1],])

plot(simFe)

###################

plotBiomass(simRe)

plotBiomass(simFe)
