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






#################

library(mizer)
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
params <- MizerParams(params_data)
params <- MizerParams(params_data, interaction = inter)
params_data_gears <- params_data
params_data_gears$gear <- c("Industrial","Industrial","Industrial",
                            "Pelagic","Beam","Otter",
                            "Beam","Otter","Beam",
                            "Otter","Otter","Otter")
params_gears <- MizerParams(params_data_gears, interaction = inter)
params <- MizerParams(params_data, interaction = inter)

effmatns <- matrix(0,12,100)
for (i in (1:100)){
  effmatns[,i] <- Fmat[,dim(Fmat)[2]]
}
rownames(effmatns) <- params@species_params$species
colnames(effmatns) <- 1:100

sim <- project(params, effort = t(effmatns), t_max = 10, dt = 0.1, t_save = 1)
plot(sim)

smallestcaught <- min(params@species_params$knife_edge_size)
gb <- getBiomass(sim,min_w=smallestcaught)
biomasstarget <- sum(gb[dim(gb)[1],])*10/9 - gb[dim(gb)[1],6]
############################
params_NS <- params

no_traits <- 10
params_traits <- set_trait_model(no_sp = no_traits+1, min_w_pp=min(params@w_full), min_w=min(params@w),max_w=max(params@w),no_w=length(params@w))
sim_traits <- project(params_traits)
# get rescaling factor for volume, to convert to above kappa_R
scalec <- 10^11/(5*10^(-3))
kappa_def <- 5*10^(-3)
k0_def <- 50
paramsK <- set_trait_model(no_sp = no_traits+1, min_w_pp=min(params_NS@w_full), min_w=min(params_NS@w),max_w=max(params_NS@w),no_w=length(params_NS@w),kappa=scalec*kappa_def,
                           k0=scalec*k0_def,f0=0.6)
wegg <- params_NS@w[1]
cut <- 0.1
paramsK <- set_trait_model(no_sp = no_traits+1, min_w_pp=min(params_NS@w_full), min_w=min(params_NS@w),max_w=max(params_NS@w),no_w=length(params_NS@w),kappa=scalec*kappa_def,
                           k0=scalec*k0_def,f0=0.6,w_pp_cutoff=cut,min_w_inf=cut,
                           max_w_inf=max(params_NS@species_params$w_inf))
simK <- project(paramsK)
plot(simK)

# rearrange mizer vignette equation (8.1) to reverse engineer k_vb
ff0 <- 0.6
kvb<-(paramsK@species_params$h)*(paramsK@species_params$alpha
)*ff0*((paramsK@species_params$w_inf
)^(-1/3))/3

prim <- paramsK@species_params[c("species","w_inf","w_mat","beta","sigma","r_max")] 
prim$k_vb <- kvb

# add in whiting
rp <- 11
whiting_no <- 6
prim_mod <- prim
#prim_mod$species[rp] <- params_data$species[whiting_no]
prim_mod$species[rp] <- no_traits+1
prim_mod$w_inf[rp] <- params_data$w_inf[whiting_no]
prim_mod$w_mat[rp] <- params_data$w_mat[whiting_no]
prim_mod$beta[rp] <- params_data$beta[whiting_no]
prim_mod$sigma[rp] <- params_data$sigma[whiting_no]
prim_mod$r_max[rp] <- params_data$r_max[whiting_no]
prim_mod$k_vb[rp] <- params_data$k_vb[whiting_no]

effmat <- matrix(0,no_traits+1,100)
effmat[rp,] <- rep(1,100)
rownames(effmat) <- 1:(no_traits+1)
colnames(effmat) <- 1:100
effmat[(no_traits+1),] <- rep(Fmat[6,dim(Fmat)[2]],dim(effmat)[2]) 


dosim <- function(capacity=10^(11),rmax_trait=1.655195e+08,rmax_whiting=5.480000e+11){
  params_data_altered <- prim_mod
  params_data_altered$r_max[1:(rp-1)] <- rmax_trait*prim_mod$r_max[1:(rp-1)]/min(prim_mod$r_max[1:(rp-1)])
  params_data_altered$r_max[rp] <- rmax_whiting
  params_prim_mod_altered <- MizerParams(params_data_altered,kappa=capacity,
                                         min_w_pp=min(params_NS@w_full),
                                         min_w=min(params_NS@w),max_w=max(params_NS@w),
                                         no_w=length(params_NS@w), w_pp_cutoff=cut)
 
  #sim_prim_mod_altered_warmup <- project(params_prim_mod_altered,effort = 1)
  sim_prim_mod_altered_warmup <- project(params_prim_mod_altered,effort = t(effmat))
  
  return(sim_prim_mod_altered_warmup)
}
plot(dosim())

traitsim <- dosim()



evalsim <- function(capacity=10^(11),rmax_trait=1.655195e+08,rmax_whiting=5.480000e+11){
  traitsim <- dosim(capacity,rmax_trait,rmax_whiting)
  gbtraits <- getBiomass(traitsim,min_w=smallestcaught)
  actualbiomass <- sum(gbtraits[dim(gbtraits)[1],1:no_traits])
  outputval <- (log(getYield(sim)[100,whiting_no]) - log(getYield(traitsim)[100,rp]))^2 + (
    (log(actualbiomass) - log(biomasstarget))^2 )
  return(outputval)
}
  
#gb <- getBiomass(sim,min_w=smallestcaught)
#biomasstarget <- sum(gb[dim(gb)[1],])*10/9 - gb[dim(gb)[1],6]



two <- function(parr=c(log10(1.655195e+08),log10(5.480000e+11))){
  return(evalsim(10^(11),10^parr[1],10^parr[2]))
}
two()


#opt <- optim(par=c(log10(1.655195e+08),log10(5.480000e+11)),two,
#             method = "SANN",
#             control = list(maxit = 100))

#opt$par
# 11.43791 12.52455

plot(dosim(capacity=10^(11), rmax_trait=10^opt$par[1], rmax_whiting=10^opt$par[2]))

 




#opt2 <- optim(par=c(opt$par[1],opt$par[2]),two,
#             method = "SANN",
#             control = list(maxit = 100))


evalsimsimple <- function(capacity=10^(11),rmax_trait=1.655195e+08,rmax_whiting=5.480000e+11){
  traitsim <- dosim(capacity,rmax_trait,rmax_whiting)
  gbtraits <- getBiomass(traitsim,min_w=smallestcaught)
  actualbiomass <- sum(gbtraits[dim(gbtraits)[1],1:no_traits])
  outputval <- (log(getYield(sim)[100,whiting_no]) - log(getYield(traitsim)[100,rp]))^2 
  return(outputval)
}

twosimple <- function(parr=c(log10(5.480000e+11))){
  return(evalsimsimple(10^(11),1.655195e+08,10^parr[1]))
}

rmaxwhitingvec <- seq(from=11,to=18,by=0.5)
ssvec <- sapply(rmaxwhitingvec,twosimple)
plot(ssvec)

rmaxwhitingvec[2]

#11.5
goodsim <- dosim(capacity=10^(11), rmax_trait=1.655195e+08, rmax_whiting=10^(11.5))
plot(goodsim)

bmg <- getBiomass(goodsim, min_w=smallestcaught)
actualtraitbm <- sum(bmg[101,1:10])
actualtraitbm
biomasstarget

# target 6.821727e+12
# actual 10^(11.39981)

# this seems ok, in the sense that it reproduces landings and plankton
# however it gets the biomass of the traits a bit wrong, 
# but if one wants to estimate this 

##################

#optsimple <- optim(par=c(log10(5.480000e+11)),twosimple,
#             method = "SANN",
#             control = list(maxit = 100))
