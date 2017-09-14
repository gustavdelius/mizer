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
Fmat2 <- Fmat
Fmat2[1,] <- Fmat[1,]
Fmat2[2,] <- Fmat[2,]
Fmat2[3,] <- Fmat[3,]
Fmat2[4,] <- Fmat[5,]
Fmat2[5,] <- Fmat[4,]
Fmat2[6,] <- Fmat[8,]
Fmat2[7,] <- Fmat[7,]
Fmat2[8,] <- Fmat[6,]
Fmat2[9,] <- Fmat[9,]
Fmat2[10,] <- Fmat[10,]
Fmat2[11,] <- Fmat[12,]
Fmat2[12,] <- Fmat[11,]
Fmat <- Fmat2

load("Landings.RData")
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")

mypar <- c(11,log10(params_data$r_max))
mysd <- 1

# initialize by running from 1957 to 2010, to generate close to realistic initial conditions 
# myn to be used for later runs

parass <- mypar
capacity <- 10^(parass[1])
rmax <- 10^(parass[2:(1+length(parass))])
dd <- params_data
dd$r_max <- rmax[1:(length(dd$r_max))]
params <- MizerParams(dd, interaction = inter, kappa=capacity)
sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =1)
#plot(sim)
myn <- sim@n[dim(sim@n)[1],,]
sim_ini <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=myn)
sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=sim_ini@n[dim(sim_ini@n)[1],,],initial_n_pp=sim_ini@n_pp[dim(sim_ini@n_pp)[1],])



#plot(sim)

##### prepare landings data

L <- t(landings)
L[is.na(L)] <- 0
log_landings <- log10(10^(-10)+L[18:(dim(L)[1]-1),])
Y <- log_landings

#k <- 7
#plot(log(getYield(sim)[,k]))
#plot(log_landings[,k])
################
model <- function(paras)
{
  capacity <- 10^(paras[1])
  rmax <- 10^(paras[2:(1+length(paras))])
  dd <- params_data
  dd$r_max <- rmax[1:(length(dd$r_max))]
  params <- MizerParams(dd, interaction = inter, kappa=capacity)
  sim_ini <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=myn)
  sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=sim_ini@n[dim(sim_ini@n)[1],,],initial_n_pp=sim_ini@n_pp[dim(sim_ini@n_pp)[1],])
  return(log10((10^(-10))+getYield(sim)))
}


# make a test time series v using PDE model

v <- model(mypar)
names(dimnames(v))[names(dimnames(v))=="sp"] <- "Species"
ym <- melt(v)
ggplot(ym) + geom_line(aes(x=time,y=value, colour=
                             Species, linetype=Species)) + scale_y_continuous(
                               name="Log Yield") + scale_x_continuous(name="Time")
v-log_landings

# make function to evaluate model cost for a particular parameter choice

fish_model_cost <- function(par=mypar,YY=log_landings, sd=mysd){
  return(sum((model(par)-YY)^2/(sd^2)))
}

##################

AA <- load(file="MCMCrun11.RData")
AA

hist(MCMCrun11$pars[,1])

mean_vec <- (1:dim(MCMCrun11$pars)[2])
sd_vec <- (1:dim(MCMCrun11$pars)[2])
for (i in (1:dim(MCMCrun11$pars)[2])){
  mean_vec[i] <- mean(MCMCrun11$pars[,i])  
  sd_vec[i] <- sd(MCMCrun11$pars[,i])  
  
}
mean_vec
#[1]  6.519508  6.225212  7.892919  7.008220  7.533711  6.449963
#[7]  1.791425  3.408887 11.008537  8.710987 13.248105  6.428867
#[13] 12.705354

sd_vec
#[1] 0.4549491 0.8688397 0.5798118 1.2575846 0.2726693 0.7010928
#[7] 2.6602041 1.2918777 0.5207947 0.7173913 1.6309563 0.8548662
#[13] 1.7299909

paras <- mean_vec
capacity <- 10^(paras[1])
rmax <- 10^(paras[2:(1+length(paras))])
dd <- params_data
dd$r_max <- rmax[1:(length(dd$r_max))]
params <- MizerParams(dd, interaction = inter, kappa=capacity)
sim_ini <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=myn)
sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=sim_ini@n[dim(sim_ini@n)[1],,],initial_n_pp=sim_ini@n_pp[dim(sim_ini@n_pp)[1],])
plot(log_landings[,8])
lines(log10(getYield(sim)[,8]))

################ done prev fit

Fmat3 <- Fmat[,24:44]
# replace odd Northen pout 2005 entry
Fmat3[3,16] <- 0.2
FmatLong <- cbind(Fmat3,Fmat3,Fmat3,Fmat3)
colnames(FmatLong) <- 1990:2073

# run 12 species model again over this shorter period

landingsshort <- landings[,41:(dim(landings)[2]-1)]
paras <- mean_vec
capacity <- 10^(paras[1])
rmax <- 10^(paras[2:(1+length(paras))])
dd <- params_data
dd$r_max <- rmax[1:(length(dd$r_max))]
params_NS <- MizerParams(dd, interaction = inter, kappa=capacity)
sim_ini <- project(params_NS, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=myn)
sim_NS <- project(params_NS, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=sim_ini@n[dim(sim_ini@n)[1],,],initial_n_pp=sim_ini@n_pp[dim(sim_ini@n_pp)[1],])
plot(log10(landingsshort[,8]))
lines(log10(getYield(sim_NS)[,8]))

##############

params_NS <- MizerParams(params_data, interaction = inter)
sim_ini <- project(params_NS, effort = t(Fmat), dt = 0.1, t_save =1)
sim_NS <- project(params_NS, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=sim_ini@n[dim(sim_ini@n)[1],,],initial_n_pp=sim_ini@n_pp[dim(sim_ini@n_pp)[1],])

# determine biomass of north sea under assumption that 12 NS species account for 90%
# of biomass
whiting_no <- 6
biomass_NS <- 10*sum(getBiomass(sim_NS)[dim(getBiomass(sim_NS))[1],])/9
# determine the biomass of the background species required to make this happen
targettraitbiomass <- biomass_NS - getBiomass(sim_NS)[dim(getBiomass(sim_NS))[1],whiting_no]

# determine smallest thing fished
smallfish <- min(params_NS@species_params$knife_edge_size)
targettraitbiomassseq <- (10*rowSums(getBiomass(sim_NS,min_w=smallfish))/9) - getBiomass(
    sim_NS)[,whiting_no]

# compute biomass of plankton

# example of verified biomass computation
#sum(sim_NS@n[45,7,]*sim_NS@params@w*sim_NS@params@dw)
#getBiomass(sim_NS)[45,7]

cdn <- (0.0001<sim_NS@params@w_full)&sim_NS@params@w_full<0.001
planktonbiomasses <- 1:dim(sim_NS@n_pp)[1]
for (tt in (1:dim(sim_NS@n_pp)[1])){
  planktonbiomasses[tt] <- sum(sim_NS@n_pp[tt,cdn]*sim_NS@params@w_full[cdn]*sim_NS@params@dw_full[cdn])
}

############

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
hybrid_Fmat <- matrix(eff, nrow=(no_traits+1), ncol=ncol(Fmat3))
hybrid_Fmat[rp,] <- Fmat3[whiting_no,]
rownames(hybrid_Fmat) <- prim_mod$species
colnames(hybrid_Fmat) <- colnames(Fmat3)


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

SS()

twodfB <- function(parr=c(11,log10(5.480000e+11))){
  return(SS(10^parr[1],1.655195e+08,10^parr[2]))
}
twodfB(c(11,log10(5.480000e+11)))
op2B <- optim(par=c(11,log10(5.480000e+11)),twodfB,lower=c(5,5),upper = c(15,15))

mgood_par <- c(op2B$par[1],log10(1.655195e+08),op2B$par[2])
# found vals  8.795080 8.218849 6.742935

guesssim4m <- dosim(10^(mgood_par[1]),10^(mgood_par[2]),10^(mgood_par[3]))
plot(guesssim4m)

plot(rownames(getYield(guesssim4m))
     ,getYield(guesssim4m)[,rp])
lines(rownames(getYield(guesssim4m))
      ,whiting_landings)

sigwhiting <- 1
sigtraits <- 1
sigplankton <- 1
SSadvanced <- function(capacity=10^(11),rmax_trait=1.655195e+08,rmax_whiting=5.480000e+11){
  mydosim <- dosim(capacity, rmax_trait,rmax_whiting  )
  planktonbiomassesT <- 1:dim(mydosim@n_pp)[1]
  for (tt in (1:dim(mydosim@n_pp)[1])){
    planktonbiomassesT[tt] <- sum(mydosim@n_pp[tt,cdn]*mydosim@params@w_full[cdn]*mydosim@params@dw_full[cdn])
  }
  return(sum((log(getYield(mydosim)[,rp])-log(whiting_landings))^2)*sigwhiting^(-2)
         +sum((log(rowSums(getBiomass(mydosim,min_w=smallfish)[,1:no_traits]))-log(targettraitbiomassseq))^2)*sigtraits^(-2)
         +sum((log(planktonbiomasses)-log(planktonbiomassesT))^2)*sigplankton^(-2))
}



# sort plankton biomass
# correct w_cut thing
# compute actual probability of old & new fits
