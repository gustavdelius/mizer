library(mizer)
library(FME)
library(ggplot2)
library(reshape2)
library(deSolve)
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library("plot3D")
library(rgl)
library("plot3Drgl")
library(optimx)
set.seed(123)
load("Fmat.RData")
Fmat <- Fmat[,11:(dim(Fmat)[2]-1)]
load("Landings.RData")
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
params_NS <- MizerParams(params_data, interaction = inter)
sim_ini <- project(params_NS, effort = t(Fmat), dt = 0.1, t_save =1)
#plot(sim)
sim_NS <- project(params_NS, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=sim_ini@n[dim(sim_ini@n)[1],,],initial_n_pp=sim_ini@n_pp[dim(sim_ini@n_pp)[1],])
plot(sim_NS)

whiting_no <- 6

biomass_NS <- 10*sum(getBiomass(sim_NS)[dim(getBiomass(sim_NS))[1],])/9
targettraitbiomass <- biomass_NS - getBiomass(sim_NS)[dim(getBiomass(sim_NS))[1],whiting_no]

no_traits <- 10
params_traits <- set_trait_model(no_sp = no_traits+1, min_w_pp=min(params_NS@w_full), min_w=min(params_NS@w),max_w=max(params_NS@w),no_w=length(params_NS@w))
sim_traits <- project(params_traits)
plot(sim_traits)
#all(params_traits@w_full==params_NS@w_full)

scalec <- targettraitbiomass/sum(getBiomass(sim_traits)[100,])
# try scalec = 10^13

kappa_def <- 0.005
k0_def <- 50
#scalec <- 10^13
paramsK <- set_trait_model(no_sp = no_traits+1, min_w_pp=min(params_NS@w_full), min_w=min(params_NS@w),max_w=max(params_NS@w),no_w=length(params_NS@w),kappa=scalec*kappa_def,
                           k0=scalec*k0_def,f0=0.6)
simK <- project(paramsK)
plot(simK)


kvb<-(paramsK@species_params$h)*(paramsK@species_params$alpha
)*ff0*((paramsK@species_params$w_inf
)^(-1/3))/3

prim <- paramsK@species_params[c("species","w_inf","w_mat","beta","sigma","r_max")] 
prim$k_vb <- kvb
#paramsK@species_params[c("w_inf","w_mat","beta","sigma","r_max","k_vb")] 
prim

params_prim <- MizerParams(prim)
sim_prim <- project(params_prim)
plot(sim_prim)

rp <- 11
prim_mod <- prim
prim_mod$species[rp] <- params_data$species[whiting_no]
prim_mod$w_inf[rp] <- params_data$w_inf[whiting_no]
prim_mod$w_mat[rp] <- params_data$w_mat[whiting_no]
prim_mod$beta[rp] <- params_data$beta[whiting_no]
prim_mod$sigma[rp] <- params_data$sigma[whiting_no]
prim_mod$r_max[rp] <- params_data$r_max[whiting_no]
prim_mod$k_vb[rp] <- params_data$k_vb[whiting_no]

params_prim_mod <- MizerParams(prim_mod)
sim_prim_mod <- project(params_prim_mod)
plot(sim_prim_mod)


# test this code from reset
# turn effort on in hypbrid model and examine default fishing gears
# change w_cut


