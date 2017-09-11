# initialization
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
sim_NS <- project(params_NS, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=sim_ini@n[dim(sim_ini@n)[1],,],initial_n_pp=sim_ini@n_pp[dim(sim_ini@n_pp)[1],])

# determine biomass of north sea under assumption that 12 NS species account for 90%
# of biomass
whiting_no <- 6
biomass_NS <- 10*sum(getBiomass(sim_NS)[dim(getBiomass(sim_NS))[1],])/9
# determine the biomass of the background species required to make this happen
targettraitbiomass <- biomass_NS - getBiomass(sim_NS)[dim(getBiomass(sim_NS))[1],whiting_no]
no_traits <- 10
# make a mizer params object for the trait based model with consistent grid points
params_traits <- set_trait_model(no_sp = no_traits+1, min_w_pp=min(params_NS@w_full), min_w=min(params_NS@w),max_w=max(params_NS@w),no_w=length(params_NS@w))
sim_traits <- project(params_traits)
# get rescaling factor for volume, to make trait biomass sensible
scalec <- targettraitbiomass/sum(getBiomass(sim_traits)[100,])
kappa_def <- 0.005
k0_def <- 50
#scalec <- 10^13
# rescale trait based model
paramsK <- set_trait_model(no_sp = no_traits+1, min_w_pp=min(params_NS@w_full), min_w=min(params_NS@w),max_w=max(params_NS@w),no_w=length(params_NS@w),kappa=scalec*kappa_def,
                           k0=scalec*k0_def,f0=0.6)
#simK <- project(paramsK)
#plot(simK)
# rearrange mizer vignette equation (8.1) to reverse engineer k_vb
ff0 <- 0.6
kvb<-(paramsK@species_params$h)*(paramsK@species_params$alpha
)*ff0*((paramsK@species_params$w_inf
)^(-1/3))/3

prim <- paramsK@species_params[c("species","w_inf","w_mat","beta","sigma","r_max")] 
prim$k_vb <- kvb

# add in whiting
rp <- 11
prim_mod <- prim
#prim_mod$species[rp] <- params_data$species[whiting_no]
prim_mod$species[rp] <- no_traits+1
prim_mod$w_inf[rp] <- params_data$w_inf[whiting_no]
prim_mod$w_mat[rp] <- params_data$w_mat[whiting_no]
prim_mod$beta[rp] <- params_data$beta[whiting_no]
prim_mod$sigma[rp] <- params_data$sigma[whiting_no]
prim_mod$r_max[rp] <- params_data$r_max[whiting_no]
prim_mod$k_vb[rp] <- params_data$k_vb[whiting_no]

# make fishing effort matrix
eff <- 0
hybrid_Fmat <- matrix(eff, nrow=(no_traits+1), ncol=ncol(Fmat))
hybrid_Fmat[rp,] <- Fmat[whiting_no,]
rownames(hybrid_Fmat) <- prim_mod$species
colnames(hybrid_Fmat) <- colnames(Fmat)
hybrid_Fmat <- hybrid_Fmat[,24:44]

#
whiting_landings <- landings[whiting_no,41:(dim(landings)[2]-1)]


dosim <- function(capacity=10^(11),rmax_trait=1.655195e+08,rmax_whiting=5.480000e+11){
  params_data_altered <- prim_mod
  params_data_altered$r_max[1:(rp-1)] <- rmax_trait*prim_mod$r_max[1:(rp-1)]/min(prim_mod$r_max[1:(rp-1)])
  params_data_altered$r_max[rp] <- rmax_whiting
  params_prim_mod_altered <- MizerParams(params_data_altered,kappa=capacity)
  hybrid_Fmat2 <- cbind(hybrid_Fmat,hybrid_Fmat,hybrid_Fmat,hybrid_Fmat)
  colnames(hybrid_Fmat2) <- 1990:2073
  sim_prim_mod_altered_warmup <- project(params_prim_mod_altered,effort = t(hybrid_Fmat2))
  sim_prim_mod_altered <- project(params_prim_mod_altered,effort = t(hybrid_Fmat),initial_n=sim_prim_mod_altered_warmup@n[dim(sim_prim_mod_altered_warmup@n)[1],,],initial_n_pp=sim_prim_mod_altered_warmup@n_pp[dim(sim_prim_mod_altered_warmup@n_pp)[1],])
  return(sim_prim_mod_altered)
}

SS <- function(capacity=10^(11),rmax_trait=1.655195e+08,rmax_whiting=5.480000e+11){
  mydosim <- dosim(capacity, rmax_trait,rmax_whiting  )
  return(sum((log(getYield(mydosim)[,rp])-log(whiting_landings))^2))
}

twodfB <- function(parr=c(11,log10(5.480000e+11))){
  return(SS(10^parr[1],1.655195e+08,10^parr[2]))
}
twodfB(c(11,log10(5.480000e+11)))
op2B <- optim(par=c(11,log10(5.480000e+11)),twodfB,lower=c(5,5),upper = c(15,15))

mgood_par <- c(op2B$par[1],log10(1.655195e+08),op2B$par[2])

guesssim4m <- dosim(10^(mgood_par[1]),10^(mgood_par[2]),10^(mgood_par[3]))
plot(guesssim4m)

plot(rownames(getYield(guesssim4m))
     ,getYield(guesssim4m)[,rp])
lines(rownames(getYield(guesssim4m))
      ,whiting_landings)

#mgood_par <- c(9.412267, 8.218849, 6.633932)
