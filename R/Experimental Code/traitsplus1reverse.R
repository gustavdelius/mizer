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

ff0 <- 0.6
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

sim_prim_mod <- project(params_prim_mod,effort=1)
plot(sim_prim_mod)

#sim_ini <- project(params_NS, effort = t(Fmat), dt = 0.1, t_save =1)

Fmat
no_traits+1
eff <- 0
hybrid_Fmat <- matrix(eff, nrow=(no_traits+1), ncol=ncol(Fmat))
hybrid_Fmat[rp,] <- Fmat[whiting_no,]
rownames(hybrid_Fmat) <- params_prim_mod@species_params$species
colnames(hybrid_Fmat) <- colnames(Fmat)
hybrid_Fmat <- hybrid_Fmat[,24:44]
sim_prim_mod <- project(params_prim_mod,effort = t(hybrid_Fmat))
plot(sim_prim_mod)

# extract landings data
# write code where background fishing effort, Rmax, kappa_R and k0
# can be varied, and we can tune to match landings data

plot(rownames(getYield(sim_prim_mod))
,getYield(sim_prim_mod)[,rp],log="xy")
whiting_landings <- landings[whiting_no,18:(dim(landings)[2]-1)]
#length(getYield(sim_prim_mod)[,rp])
#length(whiting_landings)
plot(rownames(getYield(sim_prim_mod))
      ,whiting_landings,log="xy")

colnames(landings[,18:(dim(landings)[2]-1)])
colnames(hybrid_Fmat)
##########
whiting_landings <- landings[whiting_no,41:(dim(landings)[2]-1)]

plot(rownames(getYield(sim_prim_mod))
     ,getYield(sim_prim_mod)[,rp])
lines(rownames(getYield(sim_prim_mod))
     ,100000*whiting_landings)
#getYield(sim_prim_mod)[,rp]
#whiting_landings

##################### figure out how to change r_max etc.
#####################

params_data_altered <- prim_mod
rmax_trait <- 1.655195e+08
params_data_altered$r_max[1:(rp-1)] <- rmax_trait*prim_mod$r_max[1:(rp-1)]/min(prim_mod$r_max[1:(rp-1)])
rmax_whiting <- 5.480000e+11
params_data_altered$r_max[rp] <- rmax_whiting
capacity <- 10^(11)
params_prim_mod_altered <- MizerParams(params_data_altered,kappa=capacity)

sim_prim_mod_altered <- project(params_prim_mod_altered,effort = t(hybrid_Fmat))
plot(sim_prim_mod_altered)

plot(rownames(getYield(sim_prim_mod_altered))
     ,getYield(sim_prim_mod_altered)[,rp])
lines(rownames(getYield(sim_prim_mod_altered))
      ,100000*whiting_landings)

dosim <- function(capacity=10^(11),rmax_trait=1.655195e+08,rmax_whiting=5.480000e+11){
  params_data_altered <- prim_mod
  params_data_altered$r_max[1:(rp-1)] <- rmax_trait*prim_mod$r_max[1:(rp-1)]/min(prim_mod$r_max[1:(rp-1)])
  params_data_altered$r_max[rp] <- rmax_whiting
  params_prim_mod_altered <- MizerParams(params_data_altered,kappa=capacity)
  sim_prim_mod_altered_warmup <- project(params_prim_mod_altered,effort = t(hybrid_Fmat))
  sim_prim_mod_altered <- project(params_prim_mod_altered,effort = t(hybrid_Fmat),initial_n=sim_prim_mod_altered_warmup@n[dim(sim_prim_mod_altered_warmup@n)[1],,],initial_n_pp=sim_prim_mod_altered_warmup@n_pp[dim(sim_prim_mod_altered_warmup@n_pp)[1],])
return(sim_prim_mod_altered)
}

mydosim <- dosim()

plot(rownames(getYield(mydosim))
     ,getYield(mydosim)[,rp])
lines(rownames(getYield(mydosim))
      ,1000000*whiting_landings)
# ss(theta) is sum(ln(landings_sim)-ln(landings_emp))^2

SS <- function(capacity=10^(11),rmax_trait=1.655195e+08,rmax_whiting=5.480000e+11){
  mydosim <- dosim(capacity, rmax_trait,rmax_whiting  )
  return(sum((log(getYield(mydosim)[,rp])-log(whiting_landings))^2))
}
SS(10^(11),1.655195e+08,5.480000e+11)

myfun <- function(par){
  return(SS(par[1],par[2],par[3]))
}
op <- optim(par=c(10^(11),1.655195e+08,5.480000e+11),myfun)
guesssim <- dosim(op$par[1],op$par[2],op$par[3])

plot(guesssim)

plot(rownames(getYield(guesssim))
     ,getYield(guesssim)[,rp])
lines(rownames(getYield(guesssim))
      ,whiting_landings)


log10(10^(11))
twodf <- function(parr=c(11,log10(1.655195e+08))){
  return(SS(10^parr[1],10^parr[2],5.480000e+11))
}
twodf(c(11,log10(1.655195e+08)))
op2 <- optim(par=c(11,log10(1.655195e+08)),twodf,lower=c(5,5),upper = c(15,15))

guesssim2 <- dosim(10^(op2$par[1]),10^(op2$par[2]),5.480000e+11)
plot(guesssim2)
#twodf

plot(rownames(getYield(guesssim2))
     ,getYield(guesssim2)[,rp])
lines(rownames(getYield(guesssim2))
      ,whiting_landings)
