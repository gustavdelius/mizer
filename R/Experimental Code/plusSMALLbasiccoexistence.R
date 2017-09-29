####### load external files

library(mizer)
f_history <- read.csv("./vignettes/NS_f_history.csv", row.names=1)
f_history <- as(f_history, "matrix")
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
load("Landings.RData")
landings <- t(landings)
landings[is.na(landings)] <- 0
landings <- landings[18:(dim(landings)[1]-1),]


####### setup `improved` NS model, using mike's parameters


params_data$sel_func <- "sigmoid_length"
params_data$l25 <- c(7.6, 9.8, 8.7, 10.1, 11.5, 19.8, 16.4, 19.8, 11.5,
                     19.1, 13.2, 35.3)
params_data$l50 <- c(8.1, 11.8, 12.2, 20.8, 17.0, 29.0, 25.8, 29.0, 17.0,
                     24.3, 22.9, 43.6)
params_data$a <- c(0.007, 0.001, 0.009, 0.002, 0.010, 0.006, 0.008, 0.004,
                   0.007, 0.005, 0.005, 0.007)
params_data$b <- c(3.014, 3.320, 2.941, 3.429, 2.986, 3.080, 3.019, 3.198,
                   3.101, 3.160, 3.173, 3.075)
params_data$catchability <- as.numeric(colMeans((f_history)[19:29,]))
params <- MizerParams(params_data, inter, kappa = 9.27e10)
effort_history_85_95 <- sweep(f_history,2,colMeans((f_history)[19:29,]),"/")
getnsindex <- function(s){
  matcher <- match(params_data$species,s)==1
  matcher[is.na(matcher)] <- FALSE
  return((1:length(matcher))[matcher])
}
capacity_imp <- exp(25.210)
rmax_imp <- exp(c(26.7, 26, 30.67,26.56,23.1,26.03, 22.95, 25.38,30.56, 28.375, 22.77,26.92))
rmax_imp[getnsindex("Sprat")] <- exp(26.659)
rmax_imp[getnsindex("Sandeel")] <- exp(26.008)
rmax_imp[getnsindex("Norway pout")] <- exp(30.684)
rmax_imp[getnsindex("Dab")] <- exp(23.108)
rmax_imp[getnsindex("Herring")] <- exp(26.556)
rmax_imp[getnsindex("Sole")] <- exp(22.948)
rmax_imp[getnsindex("Whiting")] <- exp(26.034)
rmax_imp[getnsindex("Plaice")] <- exp(30.562)
rmax_imp[getnsindex("Haddock")] <- exp(28.375)
rmax_imp[getnsindex("Saithe")] <- exp(26.920)
rmax_imp[getnsindex("Cod")] <- exp(22.767)
dd <- params_data
dd$r_max <- rmax_imp[1:(length(dd$r_max))]
params_data_imp <- dd
params_imp <- MizerParams(dd, interaction = inter, kappa=capacity_imp)
#simini <- project(params_imp, effort = effort_history_85_95, dt = 0.1, t_save =1)
#sim_imp <- project(params_imp, effort = effort_history_85_95, dt = 0.1, t_save =1, initial_n=simini@n[dim(simini@n)[1],,],initial_n_pp=simini@n_pp[dim(simini@n_pp)[1],])
#plot(sim_imp)
###########################
########################### make a scaled up trait based model
no_traits <- 10
scalec <- capacity_imp/(5*10^(-3))
kappa_def <- 5*10^(-3)
k0_def <- 50
paramsK <- set_trait_model(no_sp = no_traits+1, min_w_pp=min(params_imp@w_full), min_w=min(params_imp@w),max_w=max(params_imp@w),no_w=length(params_imp@w),kappa=scalec*kappa_def,
                           k0=scalec*k0_def,f0=0.6)
wegg <- params_imp@w[1]
cut <- 1
# set so max weight of trait is 10g
mintraightweight <- 4*wegg
maxtraightweight <- 5
paramsK <- set_trait_model(no_sp = no_traits+12, min_w_pp=min(params_imp@w_full), min_w=min(params_imp@w),max_w=max(params_imp@w),no_w=length(params_imp@w),kappa=scalec*kappa_def,
                           k0=scalec*k0_def,f0=0.6,w_pp_cutoff=cut,min_w_inf=mintraightweight,
                           max_w_inf=mintraightweight*(maxtraightweight/mintraightweight)^((no_traits+12-1)/(no_traits-1)))
#simK <- project(paramsK)
#plot(simK)
########################### reverse engineer vb parameters for scaled up trait model
ff0 <- 0.6
kvb<-(paramsK@species_params$h)*(paramsK@species_params$alpha
)*ff0*((paramsK@species_params$w_inf
)^(-1/3))/3
prim <- paramsK@species_params[c("species","w_inf","w_mat","beta","sigma","r_max")] 
prim$k_vb <- kvb
# note prim is the scaled up no_traits+12 trait based model, 
################### we shall now overwrite the last 12 species with those from the North sea
prim_mod <- prim
for (i in (1:12)){
  rp <- no_traits+i
  prim_mod$species[rp] <- rp
  prim_mod$w_inf[rp] <- params_data_imp$w_inf[i]
  prim_mod$w_mat[rp] <- params_data_imp$w_mat[i]
  prim_mod$beta[rp] <- params_data_imp$beta[i]
  prim_mod$sigma[rp] <- params_data_imp$sigma[i]
  prim_mod$r_max[rp] <- params_data_imp$r_max[i]
  prim_mod$k_vb[rp] <- params_data_imp$k_vb[i]
}
params_data$catchability <- as.numeric(colMeans((f_history)[19:29,]))
prim_mod$sel_func <- "sigmoid_length"
prim_mod$l25 <- rep(7.6,no_traits+12)
prim_mod$l25[(no_traits+1):(no_traits+12)] <- c(7.6, 9.8, 8.7, 10.1, 11.5, 19.8, 16.4, 19.8, 11.5,
                                                19.1, 13.2, 35.3)
prim_mod$l50 <- rep(8.1,no_traits+12)
prim_mod$l50[(no_traits+1):(no_traits+12)] <- c(8.1, 11.8, 12.2, 20.8, 17.0, 29.0, 25.8, 29.0, 17.0,
                                                24.3, 22.9, 43.6)
prim_mod$a <- rep(0.007,no_traits+12)
prim_mod$a[(no_traits+1):(no_traits+12)] <- c(0.007, 0.001, 0.009, 0.002, 0.010, 0.006, 0.008, 0.004,
                                              0.007, 0.005, 0.005, 0.007)
prim_mod$b <- rep(3.014,no_traits+12)
prim_mod$b[(no_traits+1):(no_traits+12)] <- c(3.014, 3.320, 2.941, 3.429, 2.986, 3.080, 3.019, 3.198,
                                              3.101, 3.160, 3.173, 3.075)
########## ? we still need to setup fishing gear and catchability here on prim_mod
# here prim_mod is the overwritten version
################### setup fishing efforts
effort_history_85_95_big <- matrix(0,dim(effort_history_85_95)[1],dim(effort_history_85_95)[2]+no_traits)
for (i in (1:12)){
  effort_history_85_95_big[,no_traits+i] <- effort_history_85_95[,i]
}
rownames(effort_history_85_95_big) <- rownames(effort_history_85_95)
colnames(effort_history_85_95_big) <- 1:(no_traits+12)
################ setup interaction matrix ##################
standard_interaction <- 1
biginter <- matrix(standard_interaction,nrow=22,ncol=22)
for (i in (1:12)){
  for (j in (1:12)){
    biginter[no_traits+i,no_traits+j] <- inter[i,j]
  }
}
################## run combined model
capacity <- 10*capacity_imp
rmax_trait <- 1*min(prim_mod$r_max[1:no_traits])
params_data_altered <- prim_mod
#k00 <- 10^(12)
#ratioo <- 2
#for (j in (1:no_traits)){
#  params_data_altered$r_max[j] <- k00*((ratioo)^(j-1))
#}
params_data_altered$r_max[1:no_traits] <- rmax_trait*prim_mod$r_max[1:no_traits]/min(prim_mod$r_max[1:no_traits])
#params_data_altered$r_max[1:no_traits] <- 1:(no_traits)


params_prim_mod_altered <- MizerParams(params_data_altered,
                                       kappa=capacity,w_pp_cutoff=1,interaction = biginter)
sim_prim_mod_altered <- project(params_prim_mod_altered,effort = effort_history_85_95_big)

sim_prim_mod_altered_result <- project(params_prim_mod_altered,effort = effort_history_85_95_big, 
                                       initial_n_pp=sim_prim_mod_altered@n_pp[dim(sim_prim_mod_altered@n_pp)[1],],
                                       initial_n=sim_prim_mod_altered@n[dim(sim_prim_mod_altered@n)[1],,])
plot(sim_prim_mod_altered_result) 

getBiomass(sim_prim_mod_altered_result)[44,]
####################

basicpar <- log10(c(1,10*capacity_imp,(params_data_altered$r_max[(no_traits+1):(no_traits+12)])))


runusingpar <- function(parr=basicpar){
  properpar <- 10^parr
  capacity <- properpar[2]
  rmax_trait <- properpar[1]*min(prim_mod$r_max[1:no_traits])
  params_data_altered <- prim_mod
  params_data_altered$r_max[1:no_traits] <- rmax_trait*prim_mod$r_max[1:no_traits]/min(prim_mod$r_max[1:no_traits])
  params_data_altered$r_max[(no_traits+1):(no_traits+12)] <- properpar[3:14]
  params_prim_mod_altered <- MizerParams(params_data_altered,
                                         kappa=capacity,w_pp_cutoff=1,interaction = biginter)
  sim_prim_mod_altered <- project(params_prim_mod_altered,effort = effort_history_85_95_big)
  sim_prim_mod_altered_result <- project(params_prim_mod_altered,effort = effort_history_85_95_big, 
                                         initial_n_pp=sim_prim_mod_altered@n_pp[dim(sim_prim_mod_altered@n_pp)[1],],
                                         initial_n=sim_prim_mod_altered@n[dim(sim_prim_mod_altered@n)[1],,])
  return(sim_prim_mod_altered_result)
}
plot(runusingpar())

somesim <- runusingpar()
vv <- log((getYield(somesim)+10^(-10))*10^(-6))[,((no_traits+1):(no_traits+12))]
# sum of squares
sum((vv[(dim(vv)[1]+1-dim(landings)[1]):dim(vv)[1],]-log(10^(-10)+landings[,]))^2)

#######

minme <- function(parr=basicpar){
  somesim <- runusingpar(parr)
  vv <- log((getYield(somesim)+10^(-10))*10^(-6))[,((no_traits+1):(no_traits+12))]
  return(sum((vv[(dim(vv)[1]+1-dim(landings)[1]):dim(vv)[1],]-log(10^(-10)+landings[,]))^2))
}
minme()
# USE THIS ###########
#op <- optim(par=basicpar, fn=minme, method = "SANN", control = list(maxit = 1000))
op$par
#[1]  0.7172174 12.6982774 12.6719398 10.8196778 14.7854789 11.1804614
#[7] 10.9234074 10.9920829 10.6595965 10.4257599 13.7506109 13.1286795
#[13] 11.8676254 12.7580534
basicpar
#######################################

basicpar2 <- log10(c(1,50*capacity_imp,(params_data_altered$r_max[(no_traits+1):(no_traits+12)])))
#basicpar2[13] <- 12
simw <- runusingpar(basicpar2)
getBiomass(simw)[44,]

params_imp@species_params$species

################# This will do for a 9 trait system

opfound <- c(0.7172174, 12.6982774, 12.6719398, 10.8196778,
             14.7854789, 11.1804614, 10.9234074, 10.9920829,
             10.6595965, 10.4257599, 13.7506109, 13.1286795,
             11.8676254, 12.7580534)
simww <- runusingpar(opfound)
getBiomass(simww)[44,]
plot(simww)

# note this works when max trait size is 5g, and plankton cutoff is 1g