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
paramsK <- set_trait_model(no_sp = no_traits+12, min_w_pp=min(params_NS@w_full), min_w=min(params_NS@w),max_w=max(params_NS@w),no_w=length(params_NS@w),kappa=scalec*kappa_def,
                           k0=scalec*k0_def,f0=0.6,w_pp_cutoff=wegg,min_w_inf=4*wegg,
                           max_w_inf=500000)
# need a way to choose this max_w_inf value exactly, so 
# that we 

paramsK@species_params$w_inf
simK <- project(paramsK)
plot(simK)

ff0 <- 0.6
kvb<-(paramsK@species_params$h)*(paramsK@species_params$alpha
)*ff0*((paramsK@species_params$w_inf
)^(-1/3))/3

prim <- paramsK@species_params[c("species","w_inf","w_mat","beta","sigma","r_max")] 
prim$k_vb <- kvb

######################################

# add in whiting
prim_mod <- prim
for (i in (1:12)){
  rp <- no_traits+i
  prim_mod$species[rp] <- rp
  prim_mod$w_inf[rp] <- params_data$w_inf[i]
  prim_mod$w_mat[rp] <- params_data$w_mat[i]
  prim_mod$beta[rp] <- params_data$beta[i]
  prim_mod$sigma[rp] <- params_data$sigma[i]
  prim_mod$r_max[rp] <- params_data$r_max[i]
  prim_mod$k_vb[rp] <- params_data$k_vb[i]
}

params_prim_mod_check <- MizerParams(prim_mod)
sim_prim_mod_check <- project(params_prim_mod_check,effort = 1)
plot(sim_prim_mod_check)

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
  #hybrid_Fmat2 <- cbind(hybrid_Fmat,hybrid_Fmat,hybrid_Fmat,hybrid_Fmat)
  #colnames(hybrid_Fmat2) <- 1990:2073
  #sim_prim_mod_altered_warmup <- project(params_prim_mod_altered,effort = t(hybrid_Fmat2))
  #sim_prim_mod_altered <- project(params_prim_mod_altered,effort = t(hybrid_Fmat),initial_n=sim_prim_mod_altered_warmup@n[dim(sim_prim_mod_altered_warmup@n)[1],,],initial_n_pp=sim_prim_mod_altered_warmup@n_pp[dim(sim_prim_mod_altered_warmup@n_pp)[1],])
  sim_prim_mod_altered <- project(params_prim_mod_altered,effort = 1,t_max=200)
  return(sim_prim_mod_altered)
}


#plot(dosim(10^(11),1.655195e+08,5.480000e+11))

#plot(dosim(10^(11),1.655195e+07,5.480000e+11))

simplussmall <- dosim(10^(11),1.655195e+06,5.480000e+11)
plot(simplussmall)

moresimplussmall <- dosim(10^(11),1.655195e+07,5.480000e+11)
plot(moresimplussmall)

