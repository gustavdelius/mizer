library(deSolve)
######### returne abundance test
#params <- set_scaling_model()
params <- set_scaling_model()
params@A[] <- NA
params@A[length(params@A)] <- 1
multipliers <- retune_abundance(params)

######### get scaling model
rfac <- 2
params <- set_scaling_model(max_w_inf = 5000,knife_edge_size = 10^8,kappa = 2*10^10,rfac=rfac)
#params@species_params$r_max <- params@species_params$w_mat
#params@species_params$r_max[] <- 10^50
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
    beta = 283, # = beta_gurnard from North sea. Silvia says gurnard is similar.
    sigma = 1.8, # = sigma_gurnard from North sea. Silvia says gurnard is similar.
    z0 = 0,
    alpha = 0.4, # unknown, mizer default=0.6
    erepro = 0.1, # unknown
    sel_func = "knife_edge", # not used but required
    knife_edge_size = 100, # we can choose
    gear = "knife_edge_gear",
    k = 0,
    r_max = 10^50,
    k_vb = 0.6,
    a = a_m,
    b = b_m
)
# k_vb is from 
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=790&GenusName=Mullus&SpeciesName=barbatus+barbatus&fc=332
params_out <- add_species(params, species_params, biomass = 2*10^9, min_w_observed = species_params$w_mat,rfac = rfac)
#params_out@species_params$erepro[12] <- (rfac / (rfac - 1)) * params_out@species_params$erepro[12]
#params_out@species_params$r_max[12] <-
#    (rfac - 1) * getRDI(params_out, params_out@initial_n, params_out@initial_n_pp)[12,1]


sim <- project(params_out, t_max = 250, effort = 0)
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
    beta = exp(2.4), #RLD and Blanchard thesis p 88
    sigma = 1.1, #RLD and Blanchard thesis p 88
    z0 = 0,
    alpha = 0.6, # unknown, using mizer default=0.6
    erepro = 0.1, # unknown
    sel_func = "knife_edge", # not used but required
    knife_edge_size = 100, # can choose
    gear = "knife_edge_gear",
    k = 0,
    r_max = 10^50, #why do I need r_max after combining before
    k_vb = 0.1, # from FB website below
    a = a,
    b = b
)
#k_vb <- 0.1 # from FB website below
# http://www.fishbase.org/popdyn/PopGrowthList.php?ID=30&GenusName=Merluccius&SpeciesName=merluccius&fc=184
params_out_2 <- add_species(params_out, species_params,  biomass = 2*10^9, min_w_observed = species_params$w_mat,rfac=rfac)
#params_out_2@species_params$erepro[13] <- (rfac / (rfac - 1)) * params_out_2@species_params$erepro[13]
#params_out_2@species_params$r_max[13] <-
#    (rfac - 1) * getRDI(params_out_2, params_out_2@initial_n, params_out_2@initial_n_pp)[13,1]

sim <- project(params_out_2, t_max = 55, effort = 0)
plot(sim)
params_out_2@species_params$erepro
#### mullet growth curve
mysp <- 12
gNS <- getEGrowth(params_out_2, sim@n[1, , ], sim@n_pp[dim(sim@n_pp)[1], ])[mysp,]
g_fnNS <- approxfun(params_out_2@w, gNS)
myodefunNS <- function(t, state, parameters){
    return(list(g_fnNS(state)))
}
ageNS <- (0:20)
mullet_weight <- ode(y = params_out_2@species_params$w_min[mysp], times = ageNS, func = myodefunNS, parms = 1)[,2]
plot(ageNS,params_out_2@species_params$a[mysp]*(L_inf_m*(1-exp(-params_out_2@species_params$k_vb[mysp]*ageNS)))^params_out_2@species_params$b[mysp],type="l")
lines(ageNS,mullet_weight,col="red",lty=2)

#### hake growth curve
mysp <- 13
gNS <- getEGrowth(params_out_2, sim@n[1, , ], sim@n_pp[1, ])[mysp,]
g_fnNS <- approxfun(params_out_2@w, gNS)
myodefunNS <- function(t, state, parameters){
    return(list(g_fnNS(state)))
}
ageNS <- (0:100)
hake_weight <- ode(y = params_out_2@species_params$w_min[mysp], times = ageNS, func = myodefunNS, parms = 1)[,2]
plot(ageNS,params_out_2@species_params$a[mysp]*(L_inf*(1-exp(-params_out_2@species_params$k_vb[mysp]*ageNS)))^params_out_2@species_params$b[mysp],type="l",col="blue")
lines(ageNS,hake_weight,col="purple",lty=2)

######################## all growth curves

mysp <- 1
gNS <- getEGrowth(params_out_2, sim@n[1, , ], sim@n_pp[1, ])[mysp,]
g_fnNS <- approxfun(params_out_2@w, gNS)
myodefunNS <- function(t, state, parameters){
    return(list(g_fnNS(state)))
}
ageNS <- (0:100)
weight <- ode(y = params_out_2@species_params$w_min[mysp], times = ageNS, func = myodefunNS, parms = 1)[,2]
plot(ageNS,weight,type="l",ylim=c(0,5000))
for (mysp in 2:11){
    gNS <- getEGrowth(params_out_2, sim@n[1, , ], sim@n_pp[1, ])[mysp,]
    g_fnNS <- approxfun(params_out_2@w, gNS)
    myodefunNS <- function(t, state, parameters){
        return(list(g_fnNS(state)))
    }
    weight <- ode(y = params_out_2@species_params$w_min[mysp], times = ageNS, func = myodefunNS, parms = 1)[,2]
    lines(ageNS,weight)
}



################ investigate the effect of changing fishing gears #############
t_max <- 20
eff <- 0.15
dt <- 0.01

############ control 
params_out_2_control <- params_out_2

#- control net: L50= 16.16 cm TL; sd=0.462
#- experimental net: L50= 20.50; sd= 0.331

L50 <- 16.16
sig <- 0.462
# for mullet
mysp <- 12
len <- (params_out_2_control@w/params_out_2_control@species_params$a[mysp])^(1/params_out_2_control@species_params$b[mysp])
params_out_2_control@selectivity[1,mysp,] <- 1/(1+exp(-(len-L50)/sig))
# for hake
mysp <- 13
len <- (params_out_2_control@w/params_out_2_control@species_params$a[mysp])^(1/params_out_2_control@species_params$b[mysp])
params_out_2_control@selectivity[1,mysp,] <- 1/(1+exp(-(len-L50)/sig))
sim <- project(params_out_2_control, t_max = t_max, effort = eff, dt = dt)
plot(sim)
gy_control <- getYield(sim)
############################### T90
params_out_2_t90 <- params_out_2
L50 <- 20.50 
sig <- 0.331
# for mullet
mysp <- 12
len <- (params_out_2_t90@w/params_out_2_t90@species_params$a[mysp])^(1/params_out_2_t90@species_params$b[mysp])
params_out_2_t90@selectivity[1,mysp,] <- 1/(1+exp(-(len-L50)/sig))
# for hake
mysp <- 13
len <- (params_out_2_t90@w/params_out_2_t90@species_params$a[mysp])^(1/params_out_2_t90@species_params$b[mysp])
params_out_2_t90@selectivity[1,mysp,] <- 1/(1+exp(-(len-L50)/sig))
sim <- project(params_out_2_t90, t_max = t_max, effort = eff, dt = dt)
plot(sim)
gy_t90 <- getYield(sim)

dim(gy_control[])

# comparision of gears on mullet
mysp <- 12
plot(gy_control[,mysp],type="l",ylim = c(min(c(gy_control[,mysp],gy_t90[,mysp])),max(c(gy_control[,mysp],gy_t90[,mysp]))))
lines(gy_t90[,mysp],col="red")
# comparision of gears on hake (I guess Francesc's data only applies for hake)
mysp <- 13
plot(gy_control[,mysp],type="l",ylim = c(min(c(gy_control[,mysp],gy_t90[,mysp])),max(c(gy_control[,mysp],gy_t90[,mysp]))))
lines(gy_t90[,mysp],col="red")


############################### 

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

# #53 added $a, $b and $k_vb to species params dataframe, and added column presence 
# matcher in add_species for these quantities

# #53 Have moved the part of the mullet_hake.R code that fills out 
# missing info about h,ks and gamma to the wrapper functions, so the user does 
# not have to confront it. Next I have to add warnings like "Note: 	No gamma column in species data frame so using f0, h, beta, sigma, 
#lambda and kappa to calculate it.", and clean up mullet_hake and wrapper_functions.R

# #53 Added warnings for when h,ks and gamma are missing from new species, and must be 
# calculated. Note that this computation of h uses the (1-fc/params@f0) term, 
# which I don't understand. Also, the way of getting h that I 
# penned the other day involves the assumption that k_s=0.2*h 
# which may not be the case. Can I replace is.null(species_params$a) 
#with !("a" %in% colnames(species_params)) ? Will this still suffer from 
#the same a/alpha confusion, or will this allow me to call them a and b, rather 
#then aa and bb ?

# #53 Got rid of aa and bb, by using more sensible ways of checking whether column names are 
# present.

# #53 Cleaned wrapper functions, and added min_w_observed, for biomass setting. 
# Next I need to clean mullet_hake.R, and add a clean add_species example.

# #53 cleaned mullet_hake.R, and added a clean add_species example.

# #53 Changed alpha back to default (0.4) for background, for the background species. 
# Also, now calculating h under the assumption that ks = 0.2*h. Updated hellp examples,
# and printed erepro. Have also changed the alpha value of mullet back to 0.4. The 
# erepro values for the added species seem too low, but they dont seem to change 
# much when we vary alpha.

# #21 Setup fishing gears according to Francesc's recomendations. I think he 
# only gave me the std and L50 for hake. I shall ask him for them for mullet, 
# as well as which a and b he used. I don't know which a and b to use for 
# the background species, so I just did not fish them in these experiments.

# #21 corrected a bug in the way w_min_observed was used in add_species

# #21 Readjusted the biomass of hake and mullet in the input, so the 
# initial SSB of both hake and cod are around 2000 metric tons. Lowered the 
# value of kappa from 7*10^10 to 2*10^10 so that the added species 
# have similar abundance curves to the background species. We just discovered that 
# co-existance is actually not happening.

# #18 #24 #29 #53 Since we found out the hake and mullet are not coexisting, 
# I have tried to set rmax=10. There was a problem with the way RDI was used to 
# set rmax before (rmax was made as a matrix, so rbind broke, so I had to change it to 
# getRDI()[,1]), now this error is gone, but maybe there is now an issue with how epsilon and 
# rmax are being setup for the new species. Perhaps I should delete these new 
# fragments from mullet_hake.R and put something extra in add_species to rework these 
# quantities at the end, like we did at the end of set_scaling.


# #18 #24 #29 #53 I replaced the code at the end of add species to reset erepro 
# and rmax for all species according to rfac in the same way as was done at the end of 
# add species. Unfortunately this does not get rid of the wierd dynamics.

# #18 #24 #29 #53 Changed rfac to 2, now the system seems stable. Made growth curve plots 
# for the background species
