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
kappaR2 <- capacity/12
lambda2 <- 2+0.8-(2/3)
#source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
#source("./R/Experimental Code/projectmodPREYSWITCH2.R")
source("./R/Experimental Code/project_methods_small_change.R")

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
  kappaR2 <- exp(par[1])/12
  #lambda2 <- 2+0.8-(2/3)
  #source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
  #source("./R/Experimental Code/projectmodPREYSWITCH2.R")
  source("./R/Experimental Code/project_methods_small_change.R")
  
  sim <- project(params, effort = 1, dt = 0.1, t_save =1, initial_n=init_n,initial_n_pp=init_n_pp,t_max=75)
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
  if ((all(-26<par[2:13]))&(all(par[2:13]< (-21)) )){   
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


opout <- c(23.72777, -22.48480, -24.05474, -21.34615, -23.67116, -22.75233, -22.11081, -23.05727, -22.30578, -22.51287,
           -22.27092, -21.61280, -21.41763)

### # # # op <- optim(par=opout, fn=minmesimpler, method = "Nelder-Mead", control = list(maxit = 1000))
opout <- op$par
opout
mypar
simm <- runnit(opout)
plot(simm)
gy <- getYield(simm)
log(gy[dim(gy)[1],]*10^(-6)+10^(-10))
log(landings_data90+10^(-10))

plot(log(gy[dim(gy)[1],]*10^(-6)+10^(-10)))
lines(log(landings_data90+10^(-10)))

#####################################################
#####################################################
opout <- c(23.72777, -22.48480, -24.05474, -21.34615, -23.67116, -22.75233, -22.11081, -23.05727, -22.30578, -22.51287,
           -22.27092, -21.61280, -21.41763)
simmWorks <- runnit(opout)
minmesimpler(opout)
plotSpectra(simmWorks,ylim=c(10^(-4),10^12))
simmWorks@n[dim(simmWorks@n)[1],,length(params@w[params@w<36])]
params@species_params


simmWorks@n[dim(simmWorks@n)[1],,length(params@w[params@w<21])]

nn <- simmWorks@n[dim(simmWorks@n)[1],,]
vz <- nn[2,length(params@w[params@w<22]):dim(simmWorks@n)[3]]
nn[2,length(params@w[params@w<22]):dim(simmWorks@n)[3]] <- rep(0,length(vz))


parTY <- c(23.72777, -22.48480, -24.05474, -21.34615, -23.67116, -22.75233, -22.11081, -23.05727, -22.30578, -22.51287,
           -22.27092, -21.61280, -21.41763)


dd <- params_data
dd$r_max <- rep(10^49,12)
dd$gamma <- exp(parTY[2:13])
params <- MizerParams(dd, interaction = inter, kappa=exp(parTY[1]))
#epsi <- 0.15
kappaR2 <- exp(parTY[1])/12
#lambda2 <- 2+0.8-(2/3)
#source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
#source("./R/Experimental Code/projectmodPREYSWITCH2.R")
source("./R/Experimental Code/project_methods_small_change.R")

sim <- project(params, effort = 1, dt = 0.1, t_save =1, initial_n=nn,initial_n_pp=init_n_pp,t_max=250)
plot(sim)
gy <- getYield(sim)
log(gy[dim(gy)[1],]*10^(-6)+10^(-10))
log(landings_data90+10^(-10))

sum((log(gy[dim(gy)[1],]*10^(-6)+10^(-10))-log(landings_data90+10^(-10)))^2)

minmesimpler(opout)

################################################
################################################

# h
# vb best fit
# halting algo
# data from mike ?

mypar <- c(c(25.21000, -24.04992, -24.80216, -23.06945, -25.09381, -23.67416, -23.20979,
           -24.07368, -24.12869, -24.16391, -23.81449, -22.43912, -22.69994),
           c(18.2027, 27.5160, 32.83924, 35.03807, 30.67834, 28.53952, 22.55847,
             19.37727, 14.62366, 34.23910, 61.58170, 37.39168))
           
           
           
runnitWmore <- function(par=mypar){
  dd <- params_data
  dd$r_max <- rep(10^49,12)
  dd$gamma <- exp(par[2:13])
  #dd$h <- par[14:25]
  params <- MizerParams(dd, interaction = inter, kappa=exp(par[1]))
  #epsi <- 0.15
  kappaR2 <- exp(par[1])/12
  #lambda2 <- 2+0.8-(2/3)
  #source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
  #source("./R/Experimental Code/projectmodPREYSWITCH2.R")
  source("./R/Experimental Code/project_methods_small_change.R")
  
  sim <- project(params, effort = 1, dt = 0.1, t_save =1, initial_n=init_n,initial_n_pp=init_n_pp,t_max=75)
  init_n <<- sim@n[dim(sim@n)[1],,]
  init_n_pp <<- sim@n_pp[dim(sim@n_pp)[1],]
  return(sim)
}

library(grid)
library(methods)
library(reshape2)
library(mvtnorm)
library(deSolve)

L_egg <- (params@w[1]/params_data$a)^(1/params_data$b)
W_egg <- params@w[1]
L_egg <- (W_egg/params_data$a)^(1/params_data$b)
gett0 <- function(Linf,k){
  return(log(1-L_egg/Linf)/k)
}
params_data_classic <- read.csv("./vignettes/NS_species_params.csv")

L_inf_from_params <- (params_data_classic$w_inf/params_data$a)^(1/params_data$b)
k_from_params <- params_data_classic$k_vb




minmesimplerWmore <- function(par=mypar){
  if ((all(-26<par[2:13]))&(all(par[2:13]< (-21)) )){   
    sim <- runnitWmore(par)
    gy <- getYield(sim)
    # # #
    ggB <- getEGrowth(params, sim@n[dim(sim@n)[1],,], sim@n_pp[dim(sim@n_pp)[1],])
    sols <- as.list(1:12)
    mytimes <- seq(0.5,50,0.1)
    weightsB <- matrix(0,nrow=length(mytimes),ncol=dim(sim@n)[2])
    lengthsB <- weightsB
    L_inf_best_fit <- rep(.7,12)
    k_best_fit <- rep(70,12)
    for (i in (1:12)){
      # scary try, because least squares will go wrong sometimes
      try({
        gini <- approxfun(params@w, ggB[i,])
        myodefun <- function(t, state, parameters){
          return(list(gini(state)))
        }
        weightsB[,i] <- ode(y = params@w[1], times = mytimes, func = myodefun, parms = 1)[,2]
        lengthsB[,i] <- (weightsB[,i]/params_data$a[i])^(1/params_data$b[i])
        fitTypB<-nls(Y~vbTyp(X,Linf,k),data=data.frame(X=mytimes, Y=lengthsB[,i]),start=list(Linf=L_inf_from_params[i],k=k_from_params[i]))
        L_inf_best_fit[i] <- coef(fitTypB)[["Linf"]]
        k_best_fit[i] <- coef(fitTypB)[["k"]]
      }, silent = T)
    }
    L_inf_ss <- sum((L_inf_best_fit-L_inf_from_params)^2)
    k_ss <- sum((k_best_fit-k_from_params)^2)
    
    # # #
    return(
      L_inf_ss+k_ss+sum((log(gy[dim(gy)[1],]*10^(-6)+10^(-10))-log(landings_data90+10^(-10)))^2)
    )
  }
  else {
    return(10^49)
  }
}

minmesimplerWmore() 

runnitWmore(mypar)

sm <- runnit(mypar)


#opWmore <- optim(par=mypar, fn=minmesimplerWmore, method = "SANN", control = list(maxit = 100))

#opWmore <- optim(par=mypar, fn=minmesimplerWmore, method = "Nelder-Mead", control = list(maxit = 3))

opWmore$par

oppar <- c(28.67780, -22.76756, -23.31577, -24.52235, -23.30139, -22.55133, -23.33563,
-23.64171, -23.16091, -21.94299, -22.22917, -24.26757, -23.04357,  17.20392,
28.21609,  34.51110,  37.56373,  27.83461,  27.24384,  19.95958,  19.83727,
13.86937,  30.70238,  60.83273,  36.04262)

opWmore$par
#
oppar <- c(28.67780, -22.76756, -23.31577, -24.52235, -23.30139, -22.55132, -23.33563, -23.64171,
-23.16090, -21.94299, -22.22919, -24.26759, -23.04357,  17.20392,  28.21610,  34.51111,
43.64696,  27.83462,  27.24385,  19.95958,  19.83727,  13.86937,  30.70239,  60.83273, 36.04262)

opWmore <- optim(par=oppar, fn=minmesimplerWmore, control = list(maxit = 2000))


#opWmore <- optim(par=oppar, fn=minmesimplerWmore, method = "SANN", control = list(maxit = 100))



sm <- runnitWmore(opWmore$par)
plot(sm)
minmesimplerWmore(opWmore$par)
