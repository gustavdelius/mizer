library(devtools)
#install_github("gustavdelius/mizer", ref = "master")
#install_github("drfinlayscott/mizer", ref = "master")


library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library(mizer)
#setwd("C:/Users/richa/Dropbox/Job/YorkJob/JulyMizerDB/mizer-fftclean2")


source("R/project.R")
source("R/MizerSim-class.R")


params_data <- read.csv("./vignettes/NS_species_params.csv")

inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")

####
params <- MizerParams(params_data, interaction = inter, no_w = 100)

ptm <- proc.time() 
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
Tsim <- (proc.time() - ptm)[[3]]

plot(sim)

param1 <- params

test<-project(param1,effort=0.5,t_max = 3, dt = 0.25, t_save = 1) ## run 3 years

dim(test@growth)

plot(test@growth[,1,1])
plot(test@growth[1,1,])

## 

#####

param1@w

gmid <- approxfun(param1@w, test@growth[dim(test@growth)[1],1,]) 

T <- 5

times <- seq(0, 100, by = 0.01)

#####


myodefun <- function(t, state, parameters){
 return(list(gmid(state)))
}

times <- seq(0, 10, by = 0.01)

library(deSolve)
out <- ode(y = param1@w[1], times = times, func = myodefun, parms = 1)
head(out)

plot(out[,1],out[,2])

######

T <- 100
dt <- 0.01
Nt <- floor(T/dt)
times <- (1:Nt)*dt
wlist <- times
wlist[1] <- param1@w[1]
for (t in (1:(Nt-1))){
  wlist[t+1] <- wlist[t] + dt*gmid(wlist[t])
}
plot(times,wlist)

###############
plot(param1@w, test@growth[dim(test@growth)[1],1,], log="xy") 
seqq <- seq(10^(-3), 10^3, 0.1)
points(seqq,sapply(seqq, gmid),log="xy",col="Red")


approxfun(param1@w, test@growth[dim(test@growth)[1],1,]) 


BuiltInGrowth <- getEGrowth(param1, test@n[dim(test@n)[1],,], test@n_pp[dim(test@n_pp)[1],])

gmidbuilt <- approxfun(param1@w, BuiltInGrowth[1,])

plot(param1@w, test@growth[dim(test@growth)[1],1,], log="xy") 
seqq <- seq(10^(-3), 10^3, 0.1)
points(seqq,sapply(seqq, gmid),log="xy",col="Red")

points(param1@w, BuiltInGrowth[1,], log="xy", col="Blue")


###############

approxfun(x, y = NULL,       method = "linear",
          yleft, yright, rule = 1, f = 0, ties = mean)

hh <- approxfun(test@growth[,1,1], y = NULL,       method = "linear",
          1, 3, rule = 1, f = 0, ties = mean)
hh(.2)

###

x <- 1:10
y <- rnorm(10)
par(mfrow = c(2,1))
plot(x, y, main = "approx(.) and approxfun(.)")
points(approx(x, y), col = 2, pch = "*")
points(approx(x, y, method = "constant"), col = 4, pch = "*")

f <- approxfun(x, y)
