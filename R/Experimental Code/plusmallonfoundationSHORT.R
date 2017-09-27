library(mizer)
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
load("Landings.RData")
landings <- t(landings)
landings[is.na(landings)] <- 0
landings <- landings[18:(dim(landings)[1]-1),]
params_data$sel_func <- "sigmoid_length"
params_data$l25 <- c(7.6, 9.8, 8.7, 10.1, 11.5, 19.8, 16.4, 19.8, 11.5,
                     19.1, 13.2, 35.3)
params_data$l50 <- c(8.1, 11.8, 12.2, 20.8, 17.0, 29.0, 25.8, 29.0, 17.0,
                     24.3, 22.9, 43.6)
params_data$a <- c(0.007, 0.001, 0.009, 0.002, 0.010, 0.006, 0.008, 0.004,
                   0.007, 0.005, 0.005, 0.007)
params_data$b <- c(3.014, 3.320, 2.941, 3.429, 2.986, 3.080, 3.019, 3.198,
                   3.101, 3.160, 3.173, 3.075)
f_history <- read.csv("./vignettes/NS_f_history.csv", row.names=1)
f_history <- as(f_history, "matrix")
params_data$catchability <- as.numeric(colMeans((f_history)[19:29,]))
params <- MizerParams(params_data, inter, kappa = 9.27e10)
relative_effort <- sweep(f_history,2,colMeans((f_history)[19:29,]),"/")
relative_effort[as.character(1988:1992),]
initial_effort <- matrix(relative_effort[1,],byrow=TRUE, nrow=100,
                         ncol=ncol(relative_effort), dimnames = list(1867:1966))
getnsindex <- function(s){
  matcher <- match(params_data$species,s)==1
  matcher[is.na(matcher)] <- FALSE
  return((1:length(matcher))[matcher])
}

capacity <- exp(25.210)
rmax <- exp(c(26.7, 26, 30.67,26.56,23.1,26.03, 22.95, 25.38,30.56, 28.375, 22.77,26.92))
rmax[getnsindex("Sprat")] <- exp(26.659)
rmax[getnsindex("Sandeel")] <- exp(26.008)
rmax[getnsindex("Norway pout")] <- exp(30.684)
rmax[getnsindex("Dab")] <- exp(23.108)
rmax[getnsindex("Herring")] <- exp(26.556)
rmax[getnsindex("Sole")] <- exp(22.948)
rmax[getnsindex("Whiting")] <- exp(26.034)
rmax[getnsindex("Plaice")] <- exp(30.562)
rmax[getnsindex("Haddock")] <- exp(28.375)
rmax[getnsindex("Saithe")] <- exp(26.920)
rmax[getnsindex("Cod")] <- exp(22.767)
dd <- params_data
dd$r_max <- rmax[1:(length(dd$r_max))]

params <- MizerParams(dd, interaction = inter, kappa=capacity)
simini <- project(params, effort = relative_effort, dt = 0.1, t_save =1)
sim <- project(params, effort = relative_effort, dt = 0.1, t_save =1, initial_n=simini@n[dim(simini@n)[1],,],initial_n_pp=simini@n_pp[dim(simini@n_pp)[1],])

vv <- log(getYield(sim)*10^(-6))
j <- 12
params_data$species[j]
plot(log(landings[,j]))
lines(vv[(dim(vv)[1]+1-dim(landings)[1]):dim(vv)[1],j])
vv <- log((getYield(sim)+10^(-10))*10^(-6))
# sum of squares
sum((vv[(dim(vv)[1]+1-dim(landings)[1]):dim(vv)[1],]-log(10^(-10)+landings[,]))^2)

################################

smallestcaught <- min(params_data$a*((params_data$l50)^params_data$b))
gb <- getBiomass(sim,min_w=smallestcaught)
biomasstarget <- sum(gb[dim(gb)[1],])*10/9 - gb[dim(gb)[1],6]

###########

params_NS <- params

no_traits <- 10
params_traits <- set_trait_model(no_sp = no_traits+1, min_w_pp=min(params@w_full), min_w=min(params@w),max_w=max(params@w),no_w=length(params@w))
#sim_traits <- project(params_traits)
# get rescaling factor for volume, to convert to above kappa_R
scalec <- capacity/(5*10^(-3))
kappa_def <- 5*10^(-3)
k0_def <- 50
paramsK <- set_trait_model(no_sp = no_traits+1, min_w_pp=min(params_NS@w_full), min_w=min(params_NS@w),max_w=max(params_NS@w),no_w=length(params_NS@w),kappa=scalec*kappa_def,
                           k0=scalec*k0_def,f0=0.6)
wegg <- params_NS@w[1]
cut <- 1
paramsK <- set_trait_model(no_sp = no_traits+12, min_w_pp=min(params_NS@w_full), min_w=min(params_NS@w),max_w=max(params_NS@w),no_w=length(params_NS@w),kappa=scalec*kappa_def,
                           k0=scalec*k0_def,f0=0.6,w_pp_cutoff=cut,min_w_inf=4*wegg,
                           max_w_inf=500000)
# need a way to choose this max_w_inf value exactly, so 
# that we 

paramsK@species_params$w_inf
simK <- project(paramsK)
plot(simK)

##################################################################
##################################################################

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

relative_effort_big <- matrix(0,dim(relative_effort)[1],dim(relative_effort)[2]+no_traits)
for (i in (1:12)){
    relative_effort_big[,no_traits+i] <- relative_effort[,i]
}
rownames(relative_effort_big) <- rownames(relative_effort)
colnames(relative_effort_big) <- 1:(no_traits+12)
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
  ##sim_prim_mod_altered <- project(params_prim_mod_altered,effort = t(hybrid_Fmat))
  sim_prim_mod_altered <- project(params_prim_mod_altered,effort = relative_effort_big)
  return(sim_prim_mod_altered)
}

plotSpectra(dosim(10^(11),1.655195e+08,5.480000e+11))

moresimplussmall2 <- dosim(10^(11),1.655195e+08,5.480000e+11)
plot(moresimplussmall2)
ps <- plotSpectra(moresimplussmall2)


moresimplussmall2 <- dosim(10^(11),1.655195e+9,5.480000e+11)

plot(moresimplussmall2@params@w_full,moresimplussmall2@params@w_full*moresimplussmall2@n_pp[dim(moresimplussmall2@n)[1],],log="xy",ylim=c(10^(-5),10^20))
for (i in (1:22)){
  lines(moresimplussmall2@params@w,moresimplussmall2@params@w*moresimplussmall2@n[dim(moresimplussmall2@n)[1],i,],log="xy",ylim=c(10^(-5),10^20))
}

###################

moresimplussmall2 <- dosim(10^(11),8.655195e+10,5.480000e+11)

plot(moresimplussmall2@params@w_full,moresimplussmall2@params@w_full*moresimplussmall2@n_pp[dim(moresimplussmall2@n)[1],],log="xy",ylim=c(10^(-5),10^20))
for (i in (1:22)){
  lines(moresimplussmall2@params@w,moresimplussmall2@params@w*moresimplussmall2@n[dim(moresimplussmall2@n)[1],i,],log="xy",ylim=c(10^(-5),10^20))
}
plot(moresimplussmall2)

# this is a good starting point for parameter fitting
# dosim(10^(10),7.655195e+9,5.480000e+11) 


###########

basicsmallfishsetup <- dosim(10^(11),8.655195e+10,5.480000e+11)
plot(basicsmallfishsetup)

################
#################
#capacity=10^(11),rmax_trait=1.655195e+08,rmax_whiting=5.480000e+11
# if needs be we can minimize the measure of how the community spectrum deviates 
# from the plankton community slopeif needs be, but we could alternatively just 
# set the rmaxes of the traits one at a time, so they line up, but the vignette
# suggests exact line ups dont always happen
capacity <- 10^(11)
rmax_trait <- 8.655195e+10
rmax_whiting <- 5.480000e+11
  params_data_altered <- prim_mod
  params_data_altered$r_max[1:(rp-1)] <- rmax_trait*prim_mod$r_max[1:(rp-1)]/min(prim_mod$r_max[1:(rp-1)])
  params_data_altered$r_max[rp] <- rmax_whiting
  params_prim_mod_altered <- MizerParams(params_data_altered,
                                         kappa=capacity,w_pp_cutoff=1,interaction = biginter)
  sim_prim_mod_altered <- project(params_prim_mod_altered,effort = relative_effort_big)
 plot(sim_prim_mod_altered) 
 # add facility to change the rmax of the other 11 NS species, and tune to fit landings data