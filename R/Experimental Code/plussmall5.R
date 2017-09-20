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
Fmat3 <- Fmat
# make fishing effort matrix
eff <- 0
hybrid_Fmat <- matrix(eff, nrow=(no_traits+12), ncol=ncol(Fmat3))
rownames(hybrid_Fmat) <- prim_mod$species
colnames(hybrid_Fmat) <- colnames(Fmat3)


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
  hybrid_Fmat[rp,] <- Fmat3[i,]
  
}

params_prim_mod_check <- MizerParams(prim_mod)
sim_prim_mod_check <- project(params_prim_mod_check,effort = 1)
plot(sim_prim_mod_check)

plotSpectra(sim_prim_mod_check)
nn <- sim_prim_mod_check@n
ww <- sim_prim_mod_check@params@w

plot(ww,nn[101,11,]*ww,log="xy")

whiting_landings <- landings[whiting_no,41:(dim(landings)[2]-1)]
#standard_interaction <- mean(inter)
standard_interaction <- 1
biginter <- matrix(standard_interaction,nrow=22,ncol=22)
for (i in (1:12)){
  for (j in (1:12)){
    biginter[no_traits+i,no_traits+j] <- inter[i,j]
  }
}
rownames(hybrid_Fmat) <- prim_mod$species
colnames(hybrid_Fmat) <- colnames(Fmat3)

dosim <- function(capacity=10^(11),rmax_trait=1.655195e+08,rmax_whiting=5.480000e+11){
  params_data_altered <- prim_mod
  params_data_altered$r_max[1:(rp-1)] <- rmax_trait*prim_mod$r_max[1:(rp-1)]/min(prim_mod$r_max[1:(rp-1)])
  params_data_altered$r_max[rp] <- rmax_whiting
  params_prim_mod_altered <- MizerParams(params_data_altered,
                                         kappa=capacity,w_pp_cutoff=1,interaction = biginter)
  #hybrid_Fmat2 <- cbind(hybrid_Fmat,hybrid_Fmat,hybrid_Fmat,hybrid_Fmat)
  #colnames(hybrid_Fmat2) <- 1990:2073
  #sim_prim_mod_altered_warmup <- project(params_prim_mod_altered,effort = t(hybrid_Fmat2))
  #sim_prim_mod_altered <- project(params_prim_mod_altered,effort = t(hybrid_Fmat),initial_n=sim_prim_mod_altered_warmup@n[dim(sim_prim_mod_altered_warmup@n)[1],,],initial_n_pp=sim_prim_mod_altered_warmup@n_pp[dim(sim_prim_mod_altered_warmup@n_pp)[1],])
  #sim_prim_mod_altered <- project(params_prim_mod_altered,effort = 1,t_max=200)
  sim_prim_mod_altered <- project(params_prim_mod_altered,effort = t(hybrid_Fmat))
  return(sim_prim_mod_altered)
}


moresimplussmall2 <- dosim(10^(11),1.655195e+08,5.480000e+11)
plot(moresimplussmall2)
ps <- plotSpectra(moresimplussmall2)

plot(moresimplussmall2@params@w_full,moresimplussmall2@params@w_full*moresimplussmall2@n_pp[201,],log="xy",ylim=c(10^(-5),10^20))
for (i in (1:22)){
  lines(moresimplussmall2@params@w,moresimplussmall2@params@w*moresimplussmall2@n[201,i,],log="xy",ylim=c(10^(-5),10^20))
}


#moresimplussmall2@params@species_params
#plot(getBiomass(moresimplussmall2)[201,],log="y")
#gb2 <- getBiomass(moresimplussmall2,max_w=1)
#gb2[201,]
#getRDI(moresimplussmall2@params, moresimplussmall2@n[201,,], moresimplussmall2@n_pp[201,])
#getSSB(moresimplussmall2)[201,]
#gg <- getEGrowth(moresimplussmall2@params, moresimplussmall2@n[201,,], moresimplussmall2@n_pp[201,])
#plot(moresimplussmall2@params@w,gg[13,],log="xy")
#plot(moresimplussmall2@params@w,getFeedingLevel(moresimplussmall2)[201,13,],log="xy")
#plot(moresimplussmall2@params@w,moresimplussmall2@n[201,12,],log="xy",ylim=c(10^5,10^14))

plot(moresimplussmall2@params@w_full,moresimplussmall2@n_pp[201,],log="xy",ylim=c(10^(-5),10^20))
for (i in (1:22)){
  lines(moresimplussmall2@params@w,moresimplussmall2@n[201,i,],log="xy",ylim=c(10^(-5),10^20))
}
gm2 <- getM2(moresimplussmall2@params,moresimplussmall2@n[201,,],moresimplussmall2@n_pp[201,])
plot(moresimplussmall2@params@w,gm2[13,],log="x")


# if one has a system where the plankton lines up, 
# and one scales the k0 by a particular factor then one should 
# change the kappa_R be the same factor in order to make things work.


#########

# cut back the plankton spectrum, and setup the new species so that it takes up the same biomass