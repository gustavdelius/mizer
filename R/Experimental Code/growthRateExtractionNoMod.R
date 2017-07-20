library(devtools)
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library(mizer)
library(deSolve)
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
param1 <- MizerParams(params_data, interaction = inter, no_w = 100)
test<-project(param1,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)
plot(test)

# size history outputs weight vs time curve, the inputs are mizerparams, mizersim, 
# and the times one wishes to obtain the weights at.

sizeHistory <- function(param1, test, times=seq(0,1,0.1) ){
ns <- dim(test@n)[2]
gg <- getEGrowth(param1, test@n[dim(test@n)[1],,], test@n_pp[dim(test@n_pp)[1],])
sols <- as.list(1:ns)
MM <- matrix(0,nrow=length(times),ncol=dim(test@n)[2])
for (i in (1:ns)){
  gini <- approxfun(param1@w, gg[i,])
  
  myodefun <- function(t, state, parameters){
    return(list(gini(state)))
  }
  out <- ode(y = param1@w[1], times = times, func = myodefun, parms = 1)
  MM[,i] <- out[,2]
}
return(MM)
}

# output is MM where MM[t,i] is the size of a species i fish at time times[t], where the
# fish is egg sized at times[1], and the growth rates used to computed using the final state 
# of the mizersim object test, that was generated using mizerparams object param1

#################
mytimes <- seq(0,100,0.1)
MM <- sizeHistory(param1,test,times=mytimes)
plot(mytimes,MM[,11])

hh<-param1@species_params
#species
hh[,1]
# empirical max weight
hh[,2]
# my prediction of weight after 10 yrs
MM[dim(MM)[1],]

plot(hh[,2])
points(MM[dim(MM)[1],], col="Red")

hh[,1][[11]]

hh[11,2]

plot(mytimes,MM[,11])


MM[dim(MM)[1],]-hh[,2]

for (i in 1:12){
  plot(mytimes,MM[,i])
}


####

GG<-  getEGrowth(param1, test@n[dim(test@n)[1],,], test@n_pp[dim(test@n_pp)[1],])
dim(GG)
GG[1,]
plot(param1@w,GG[1,],log="xy")

GG[1,sum(GG[1,]>0)]
param1@w[sum(GG[1,]>0)]
param1@w[sum(GG[1,]>0)+1]

#######

MM <- sizeHistory(param1,test,times=mytimes)

dim(MM)

plot(mytimes,MM[,1])

length(param1@w)

length(MM[1:1000,1])

dim(param1@species_params)

param1@species_params[1,]

param1

####

weightTo

a <- 0.058
b <- 3.08

w <- a*L^b

w/a = L^b

exp(log(MM[dim(MM)[1],1]/a)/b)

a*8^b

param1@species_params


plot(mytimes,MM[,4])

a <- 0.006
b <- 3.05
exp(log(MM[dim(MM)[1],4]/a)/b)


MM[,4]

plot(param1@w,test@n[4,4,],log="xy")
plot(test)

