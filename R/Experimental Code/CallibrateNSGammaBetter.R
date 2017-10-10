library(mizer)
library(plyr)
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
load("Landings.RData")
landings <- t(landings)
landings[is.na(landings)] <- 0
landings <- landings[18:(dim(landings)[1]-1),]

landings_data90 <- landings[24,]

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
params_data$catchability <- as.numeric(f_history["1990",])
#params <- MizerParams(params_data, inter, kappa = 9.27e10)
dd <- params_data
dd$r_max <- rep(10^49,12)
capacity <- exp(25.210)

params <- MizerParams(dd, interaction = inter, kappa=capacity)
epsi <- 0.03
kappaR2 <- capacity
lambda2 <- 2+0.8-(2/3)
source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
source("./R/Experimental Code/projectmodPREYSWITCH2.R")
simini <- project(params, effort = 1, dt = 0.1, t_save =1)
first_init_n <- simini@n[dim(simini@n)[1],,]
first_init_n_pp <- simini@n_pp[dim(simini@n_pp)[1],]
sim <- project(params, effort = 1, dt = 0.1, t_save =1, initial_n=first_init_n,initial_n_pp=first_init_n_pp)


mypar <- c(log(88830864522),log(c(3.591300e-11, 1.692620e-11, 9.573345e-11, 1.264439e-11, 5.229280e-11, 8.319884e-11,
                                  3.506987e-11, 3.319288e-11, 3.204404e-11, 4.544627e-11, 1.798108e-10, 1.385290e-10)))

init_n <- first_init_n
init_n_pp <- first_init_n_pp


runnit <- function(par=mypar){
  dd <- params_data
  dd$r_max <- rep(10^49,12)
  dd$gamma <- exp(par[2:13])
  params <- MizerParams(dd, interaction = inter, kappa=exp(par[1]))
  #epsi <- 0.15
  kappaR2 <- exp(par[1])
  #lambda2 <- 2+0.8-(2/3)
  source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
  source("./R/Experimental Code/projectmodPREYSWITCH2.R")
  sim <- project(params, effort = 1, dt = 0.1, t_save =1, initial_n=init_n,initial_n_pp=init_n_pp,t_max=100)
  init_n <<- sim@n[dim(sim@n)[1],,]
  init_n_pp <<- sim@n_pp[dim(sim@n_pp)[1],]
  return(sim)
}

###########################

simm <- runnit()
plot(runnit())

gy <- getYield(simm)
gy[dim(gy)[1],]*10^(-6)
landings_data90

minmesimpler <- function(par=mypar){
  #if ((all(-26<par[2:13]))&(all(par[2:13]< (-22)) )){
    if ((all(-26<par[2:13]))&(all(par[2:13]< (-20)) )){   
    simm <- runnit(par)
    gy <- getYield(simm)
    
    return(
      sum((log(gy[dim(gy)[1],]*10^(-6)+10^(-10))-log(landings_data90+10^(-10)))^2)
    )
  }
  else {
    return(10^49)
  }
}

opout <- c(25.32141, -23.48951, -25.39341, -23.11366, -25.09401, -23.44043, -23.19967, -24.57311,
            -24.47119, -24.35656, -24.45090, -23.16593, -22.53987)
opout <- c(24.97661, -23.70650, -25.21430, -22.86702, -24.72171, -23.70656, -23.24409, -23.79086, -22.43615, -22.09472,
           -23.50825, -22.43677, -22.78187)
# this gives a good fit when chi=0.15, except it over estimates cod

opout <- c(25.12258, -23.24448, -25.41143, -22.87238, -24.77804, -23.42233, -22.70323, -23.65563, -22.45981, -21.38096, -23.28547,
           -22.43764, -22.81616)

op <- optim(par=opout, fn=minmesimpler, method = "SANN", control = list(maxit = 5000))
opout <- op$par
opout
mypar
simm <- runnit(opout)
plot(simm)
gy <- getYield(simm)
log(gy[dim(gy)[1],]*10^(-6)+10^(-10))
log(landings_data90+10^(-10))

