library(deSolve)
######### returne abundance test
#params <- set_scaling_model()
params <- set_scaling_model()
params@A[] <- NA
params@A[length(params@A)] <- 1
retune_abundance(params)

######### get scaling model

params <- set_scaling_model(max_w_inf = 5000,alpha = 0.6)
params@species_params$r_max <- params@species_params$w_mat
params@species_params$r_max[] <- 10^50
params@A[] <- NA

######### add mullet
# some data from fishbase at 
# http://www.fishbase.org/summary/Mullus-barbatus+barbatus.html
# some parameter info is in this table
# https://www.dropbox.com/s/iqiydcrxqrx0k0w/paramsTable.jpg?dl=0
# length to weight conversion constants from 
# http://www.fishbase.org/popdyn/LWRelationshipList.php?ID=790&GenusName=Mullus&SpeciesName=barbatus+barbatus&fc=332
a_m <- 0.0085
b_m <- 3.11
# asymptotic length from
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=790&GenusName=Mullus&SpeciesName=barbatus+barbatus&fc=332
L_inf_m <- 24.3
# length at maturity from 
# http://www.fishbase.org/summary/Mullus-barbatus+barbatus.html
L_mat <- 11.1
species_params <- data.frame(
    species = "mullet",
    w_min = 0.001, # mizer's default egg weight, used in NS
    # w_inf = 251.94, #is the old value we used. Where is it from ? It differs to below
    w_inf = a_m*L_inf_m^b_m, # from fishbase
    # w_mat = 16.48, #is the old value we used. Where is it from ? It differs to below
    w_mat = a_m*L_mat^b_m, # from fishbase
    h = NA, # will compute this later
    #ks = 4,
    ks = NA, # unknown, so we setup mizer's default of ks=0.2*h below
    beta = 283, # = beta_gurnard from North sea. Silvia says gurnard is similar.
    sigma = 1.8, # = sigma_gurnard from North sea. Silvia says gurnard is similar.
    z0 = 0,
    #alpha = 0.4, # unknown, set same as set_scaling default. Normal mizer default=0.6
    alpha = 0.6, # unknown, set same as set_scaling default. Normal mizer default=0.6
    erepro = 0.1, # unknown
    sel_func = "knife_edge", # not used but required
    knife_edge_size = 100, # we can choose
    gear = "knife_edge_gear",
    k = 0,
    gamma = NA,
    w_min_idx = NA,
    r_max = 10^50,
    k_vb = 0.6,
    aa = a_m,
    bb = b_m
)
# k_vb is from 
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=790&GenusName=Mullus&SpeciesName=barbatus+barbatus&fc=332
fc <- 0.2/species_params$alpha
species_params$h <- 3*species_params$k_vb*(species_params$w_inf^(1/3))/(
    species_params$alpha*params@f0*(1-fc/params@f0))
# The above definition of h differs from old eqn (8.1) from vignette, see notes
# https://www.dropbox.com/s/j9yajrx3d2t0zub/hDefinition.jpg?dl=0
species_params$ks <- 0.2*species_params$h # mizer's default setting
ae <- sqrt(2*pi) * species_params$sigma * species_params$beta^(params@lambda-2) * exp((params@lambda-2)^2 * species_params$sigma^2 / 2)
species_params$gamma <- (species_params$h / (params@kappa * ae)) * (params@f0 / (1 - params@f0))
species_params$w_min_idx <- sum(params@w<=species_params$w_min)
#@ params_out <- add_species(params, species_params, mult = 5.5 * 10 ^ (8))
params_out <- add_species(params, species_params, biomass = 3070953023)
sim <- project(params_out, t_max = 5, effort = 0)
plot(sim)

############# add hake 
# Merluccius merluccius  (European hake)
# http://www.fishbase.org/summary/Merluccius-merluccius.html
#! Currently hake and mullet are both using the same feeding level. What to do about it ?
# length to weight conversion: w=a*L^b, a = 0.0046, b = 3.12
# http://www.fishbase.org/popdyn/LWRelationshipList.php?ID=30&GenusName=Merluccius&SpeciesName=merluccius&fc=184
a <- 0.0046
b <- 3.12
# characteristic weights from
# http://www.fishbase.org/Reproduction/MaturityList.php?ID=30&GenusName=Merluccius&SpeciesName=merluccius&fc=184
# http://www.fishbase.org/graph/graphLengthFM01.php?RequestTimeout=50000&ID=30&genusname=Merluccius&speciesname=merluccius&fc=184&gm_lm=29.832069860776&gm_loo=81.220460002349
L_inf <- 81.2
L_mat <- 29.83
# Some information below is from Richard Law's document (RLD) at
# https://www.dropbox.com/s/g701wgcnhr12qpg/species%20%282%29.pdf?dl=0

species_params <- data.frame(
    species = "hake",
    w_min = 0.001, # mizer default
    w_inf = a*L_inf^b, # from fishbase
    w_mat = a*L_mat^b, # from fishbase
    h = NA, # will compute this later
    # ks = 1/2, # unknown, mizer default =0.2*h
    ks = NA, # defined later as mizer default ks=0.2*h
    beta = exp(2.4), #RLD and Blanchard thesis p 88
    sigma = 1.1, #RLD and Blanchard thesis p 88
    z0 = 0,
    #alpha = 0.4, # unknown, set same as set_scaling default. Normal mizer default=0.6
    alpha = 0.6, # unknown, set same as set_scaling default. Normal mizer default=0.6
    erepro = 0.1, # unknown
    sel_func = "knife_edge", # not used but required
    knife_edge_size = 100, # can choose
    gear = "knife_edge_gear",
    k = 0,
    gamma = NA,
    w_min_idx = NA,
    r_max = 10^50, #why do I need r_max after combining before
    k_vb = 0.1, # from FB website below
    aa = a,
    bb = b
)
#k_vb <- 0.1 # from FB website below
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=30&GenusName=Merluccius&SpeciesName=merluccius&fc=184
fc <- 0.2/species_params$alpha
species_params$h <- 3*species_params$k_vb*(species_params$w_inf^(1/3))/(species_params$alpha*params@f0*
                                                             (1-fc/params@f0))
species_params$ks <- 0.2*species_params$h # mizer's default setting
ae <- sqrt(2*pi) * species_params$sigma * species_params$beta^(params@lambda-2) * exp((params@lambda-2)^2 * species_params$sigma^2 / 2)
species_params$gamma <- (species_params$h / (params@kappa * ae)) * (params@f0 / (1 - params@f0))
species_params$w_min_idx <- sum(params@w<=species_params$w_min)


#@ params_out_2 <- add_species(params_out, species_params, mult = 5.5 * 10 ^ (8))
params_out_2 <- add_species(params_out, species_params, biomass = 427977180)
sim <- project(params_out_2, t_max = 5, effort = 0)
plot(sim)

#### mullet growth curve
mysp <- 12
#gNS <- getEGrowth(params_out_2, sim@n[dim(sim@n)[1], , ], sim@n_pp[dim(sim@n_pp)[1], ])[mysp,]
gNS <- getEGrowth(params_out_2, sim@n[1, , ], sim@n_pp[dim(sim@n_pp)[1], ])[mysp,]

g_fnNS <- approxfun(params_out_2@w, gNS)
myodefunNS <- function(t, state, parameters){
    return(list(g_fnNS(state)))
}
ageNS <- (0:20)
mullet_weight <- ode(y = params_out_2@species_params$w_min[mysp], times = ageNS, func = myodefunNS, parms = 1)[,2]
plot(ageNS,params_out_2@species_params$aa[mysp]*(L_inf_m*(1-exp(-params_out_2@species_params$k_vb[mysp]*ageNS)))^params_out_2@species_params$bb[mysp],type="l")
lines(ageNS,mullet_weight,col="red",lty=2)

params_out_2@species_params$w_min[mysp]
#### hake growth curve
mysp <- 13
#gNS <- getEGrowth(params_out_2, sim@n[dim(sim@n)[1], , ], sim@n_pp[dim(sim@n_pp)[1], ])[mysp,]
gNS <- getEGrowth(params_out_2, sim@n[1, , ], sim@n_pp[1, ])[mysp,]

g_fnNS <- approxfun(params_out_2@w, gNS)
myodefunNS <- function(t, state, parameters){
    return(list(g_fnNS(state)))
}
ageNS <- (0:100)
hake_weight <- ode(y = params_out_2@species_params$w_min[mysp], times = ageNS, func = myodefunNS, parms = 1)[,2]
#plot(ageNS,a*(L_inf*(1-exp(-k_vb*ageNS)))^b,col="blue",type="l")
plot(ageNS,params_out_2@species_params$aa[mysp]*(L_inf*(1-exp(-params_out_2@species_params$k_vb[mysp]*ageNS)))^params_out_2@species_params$bb[mysp],type="l",col="blue")
lines(ageNS,hake_weight,col="purple",lty=2)

################ investigate the effect of changing fishing gears #############
eff <- 0.15
sim <- project(params_out_2, t_max = 5, effort = eff)
plot(sim)
gyA <- getYield(sim)
gyA[dim(gyA)[1],]

params_out_2B <- params_out_2
# change the knife edge on hake to 50g
params_out_2B@selectivity[1,13,] <- 1*params_out_2@w>50
sim <- project(params_out_2B, t_max = 5, effort = eff)
plot(sim)
gyB <- getYield(sim)
gyB[dim(gyB)[1],]

# notice there is a slight increase in the catch of mullet 
# in the latter system, because we are fishing its predator, hake, 
# with a more disruptive gear. 

# #18 #24 #29 Have got code that holds hake and mullet (pushed to inst/mullet_hake.R in adsp branch). 
# Next I want to retune the abundance multipliers to be more reasonable, and do some experiments with fishing gears.

# #18 #24 #29 filled out where some hake and mullet parameters are from
# am having difficulty reproducing the mullet characteristic weights I used before.
# also, using the default ks=0.2*h the mullet seems to have too much energy, 
# and curls up, while the hake does not seem to have enough energy and its 
# biomass curve slopes down. Need to figure out what I can tune. 

# #18 #24 #29 added growth curves. Having issues because I want to make max_w_inf=5000 
# in the initial set_scaling model, but when I do, the code breaks.

# #18 #24 #29 Have fixed wrapper functions, so that set_scaling has an integer number 
# of weight bins, and now it works with general max_w_inf, so now I can run the 
# system properly, but it looks like the hake is still growing too slow. 
# Next job is to plot the true VB growth curves alongside.

# #18 #24 #29 Plotted growth curves vs VB curves (with t0 =0). Agreement is not good. 
# also corrected a typo about specifying W_mat for mullet.

# #18 #24 #29 Modified definition of h to account for fc. Now growth curves look 
# better, except that hake has stunted growth.

# #18 #24 #29 Now I use the initial conditions (rather than state after sim) 
# to define the growth rates used for growth curves.

# #18 #24 #29 Still dont understand why I could not match the Hake growth curve 
# when I used the original settings. Now I run the code again with all alpha=0.6 
# (which is the mizer default), and the growth curves seem to match much better with 
# the VB growth curves now. There is a slight discrepency because w_inf does not 
# lie upon a grid point. I guess the w_inf etc. should be coerced to lie upon 
# a grid point. Also, now with alpha =0.6, I do not understand why now it keeps 
# repeatedly warning "Note: Negative background mortality rates overwritten with zeros"
# I guess this warning message should be suppressed after first print.

# #18 #24 #29 corrected previous typo with setup of alpha.

# #21 Added very basic fishing gears, to show slight increase in 
# the yield of mullet that comes as a knock-on effect of changing 
# to use a more aggressive fishing gear on the hake. 

# #53 added $aa, $bb and $k_vb to species params dataframe, and added column presence 
# matcher in add_species for these quantities

