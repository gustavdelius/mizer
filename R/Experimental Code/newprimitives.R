# try and tune trait based model to match biomass etc. for simple vignette
# version of mizer

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
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
plot(sim)

# get trait based model, scale it up until its plankton kappa_R
# is the same as that of the above model

# add whiting

# we want to vary R_max of whiting, and k0 under traits, so we match
# the 90% biomass, and the biomass of the whiting

params_NS <- params

#################################

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

dosim <- function(capacity=10^(11),rmax_trait=1.655195e+08,rmax_whiting=5.480000e+11){
  params_data_altered <- prim_mod
  params_data_altered$r_max[1:(rp-1)] <- rmax_trait*prim_mod$r_max[1:(rp-1)]/min(prim_mod$r_max[1:(rp-1)])
  params_data_altered$r_max[rp] <- rmax_whiting
  params_prim_mod_altered <- MizerParams(params_data_altered,kappa=capacity,
                                         min_w_pp=min(params_NS@w_full),
                                         min_w=min(params_NS@w),max_w=max(params_NS@w),
                                         no_w=length(params_NS@w), w_pp_cutoff=cut)
#  sim_prim_mod_altered_warmup <- project(params_prim_mod_altered,effort = c(rep(0,no_traits),1))
  effmat <- matrix(0,no_traits+1,100)
  effmat[rp,] <- rep(1,100)
  rownames(effmat) <- 1:(no_traits+1)
  colnames(effmat) <- 1:100
  #sim_prim_mod_altered_warmup <- project(params_prim_mod_altered,effort = 1)
  sim_prim_mod_altered_warmup <- project(params_prim_mod_altered,effort = t(effmat))
  
    return(sim_prim_mod_altered_warmup)
}
plot(dosim())
