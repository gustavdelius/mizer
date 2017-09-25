# makes output1
load("for_richard.Rdata")
params_data$species

# mikes table, species ordering
mikes_species_list <- c("Sprat",
  "Sandeel",
  "N.pout",
  "Dab",
  "Herring",
  "Gurnard",
  "Sole",
  "Whiting",
  "Plaice",
  "Haddock",
  "Saithe",
  "Cod"
  )

getmsindex <- function(s){
  matcher <- match(mikes_species_list,s)==1
  matcher[is.na(matcher)] <- FALSE
  return((1:length(matcher))[matcher])
}
getmsindex("Cod")

#################

library(mizer)
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
#load("Fmat.RData")
#Fmat <- Fmat[,11:(dim(Fmat)[2]-1)]
load("Landings.RData")
#params_data$catchability <- as.numeric(colMeans((f_history)[19:29,]))
L <- t(landings)
L[is.na(L)] <- 0
log_landings <- log10(10^(-10)+L[18:(dim(L)[1]-1),])
#####
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
#params_data$catchability <- as.numeric(f_history["1990",])
params_data$catchability <- as.numeric(colMeans((f_history)[19:29,]))
params <- MizerParams(params_data, inter, kappa = 9.27e10)
#relative_effort <- sweep(f_history,2,f_history["1990",],"/")
relative_effort <- sweep(f_history,2,colMeans((f_history)[19:29,]),"/")

load("Fmat.RData")
Fmat <- Fmat[,11:(dim(Fmat)[2]-1)]
Fmat_proper <- Fmat
for (i in (1:12)){
  Fmat_proper[i,] <- 
    Fmat[getmsindex(params_data$species[i]),]
}

relative_effort <- sweep(Fmat_proper,2,colMeans((f_history)[19:29,]),"/")

relative_effort[as.character(1988:1992),]
initial_effort <- matrix(relative_effort[1,],byrow=TRUE, nrow=100,
                         ncol=ncol(relative_effort), dimnames = list(1867:1966))
#relative_effort <- rbind(initial_effort,relative_effort)

load("for_richard.RData")
# find entry with Haddock in it, and write in the Rmax manually, while checking against for_richard.Rdata

### more accurate part specification

getnsindex <- function(s){
  matcher <- match(params_data$species,s)==1
  matcher[is.na(matcher)] <- FALSE
  return((1:length(matcher))[matcher])
}
getnsindex("Haddock")
getnsindex("Cod")

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

################## output1

#columns 1:12 are log(Rmax) (notice the new ordering, see paper for ordering)
#column 13 is log(carrying capacity)
#columns 14:25 are the fishing mortality rates, relative to baseF I believe,
#that we used to run the model into equilibrium (again in the order of the paper)
#column 26 is the fishing mortality, again relative to baseF of Norway
#pout in 2005 (it is 0 in Fmat but there are landings that year,
#albeit no very many)
lenout <- dim(output1)[1]
ln_r_max_m <- colMeans(output1[10000:lenout,1:12])
ln_r_max_m_proper <- r_max_m
for (i in (1:12)){
  ln_r_max_m_proper[i] <- ln_r_max_m[getmsindex(params_data$species[i])]
}
ln_cap_proper <- mean(output1[10000:lenout,13])

starting_effort_m <- colMeans(output1[10000:lenout,14:25])
starting_effort_m_proper <- starting_effort_m
for (i in (1:12)){
  starting_effort_m_proper[i] <- starting_effort_m[getmsindex(params_data$species[i])]
}
names(starting_effort_m_proper) <- params_data$species


rmax <- exp(ln_r_max_m_proper)

capacity <- exp(ln_cap_proper)

dd <- params_data
dd$r_max <- rmax[1:(length(dd$r_max))]

# ! use basic or new, are they much different ? not much
params <- MizerParams(dd, interaction = inter, kappa=capacity)


####################################

#simini <- project(params, effort = relative_effort, dt = 0.1, t_save =1)
simini <- project(params, effort = starting_effort_m_proper/colMeans((f_history)[19:29,]), dt = 0.25, t_save =1, t_max=300)



sim <- project(params, effort = t(relative_effort), dt = 0.25, t_save =1, initial_n=simini@n[dim(simini@n)[1],,],initial_n_pp=simini@n_pp[dim(simini@n_pp)[1],])

vv <- log10((10^(-10))+getYield(sim))


j <- 12
params_data$species[j]
plot(log_landings[,j])
lines(vv[(dim(vv)[1]+1-dim(log_landings)[1]):dim(vv)[1],j]-6)


#######################
# now we need to repermute rmax, and use Fmat instead
