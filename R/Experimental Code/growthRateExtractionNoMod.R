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


plot(mytimes,MM[,4])

weight <- MM[,4]
lengths <- sapply(MM[,4], function(x) exp(log(x/a)/b))

#winf, gamma, h

plot(mytimes,lengths)


################

param1@species_params[4,]

#w inf
param1@species_params[4,2]

# h 
param1@species_params[4,14]

# gamma
param1@species_params[4,15]

lengthlist <- function(winf, h, gamma, param1, mytimes = seq(0,100,0.1)){
  param1b <- param1
  param1b@species_params[4,2] <- winf
  param1b@species_params[4,14] <- h
  param1b@species_params[4,15] <- gamma
  test<-project(param1b,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)
  MM <- sizeHistory(param1b,test,times=mytimes)
  lengths <- sapply(MM[,4], function(x) exp(log(x/a)/b))
  return(lengths)
  
}

lengthlist(winf,h,gamma,param1)==lengthlist(winf+500,h,gamma,param1)

plot(mytimes,lengthlist(winf+500,h,gamma,param1))

plot(mytimes,lengthlist(winf,h,gamma,param1)-lengthlist(winf,h,gamma+10^3,param1))


#w inf
winf <- param1@species_params[4,2]

h <- param1@species_params[4,14]

gamma <- param1@species_params[4,15]


plot(mytimes,lengthlist(winf,h,gamma,param1))

points(mytimes,lengthlist(winf+5,h,gamma,param1),col="Red")

tail(lengthlist(winf,h,gamma,param1))

tail(lengthlist(winf+50,h,gamma,param1))

#FSA
package


########################

param1b <- param1
param1b@species_params[4,2] <- winf+500
param1b@species_params[4,14] <- h
param1b@species_params[4,15] <- gamma
test<-project(param1b,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)
MM <- sizeHistory(param1b,test,times=mytimes)
lengths <- sapply(MM[,4], function(x) exp(log(x/a)/b))
lengths

param1==param1b


MM[,4]


sizeHistory(param1b,test,times=mytimes)==sizeHistory(param1,test,times=mytimes)
param1b@species_params==param1@species_params


test1<-project(param1,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)
test1b<-project(param1b,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)
test1@n==test1b@n

head(getEGrowth(param1, test1@n[dim(test1@n)[1],,], test1@n_pp[dim(test1@n_pp)[1],])[4,])

head(getEGrowth(param1b, test1b@n[dim(test1b@n)[1],,], test1b@n_pp[dim(test1b@n_pp)[1],])[4,])


############

# fix length list for time variable 

dats <- data.frame(X=mytimes, Y= lengths)

vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*X)))}

# obs_Linf = 35 cm

obs_Linf <- 29.3
obs_k <- 0.6
fitTyp<-nls(Y~vbTyp(X,Linf,k),data=dats,start=list(Linf=obs_Linf,k=obs_k))

load("Herring.Rdata")
out

library(mvtnorm)

SIGMA<-cov(out[,1:2])

MU<-colMeans(out[,1:2])



## make x the mizer prediction for the VBGF parameters

dmvnorm(c(.58, 301),MU,SIGMA,log=T)

dmvnorm(MU,MU,SIGMA,log=T)


head(out)

#col = k
# col2 = linf
# col 3= t0 ignore
#col 4 = sigma (squred), ignore

####
# set up fishing rates and make size history based on dynamic 
# growth rates

# get likelyhood function by dmvnorm (by calculating sigma and mu for real data points from
# p(theta|k), as mike sent)

# use vbTyp to get VB parameters that best fit data when mizer is run 
# on theta, and evaluate dmvnorm for these parameters, this is our likelihood

# write mcmc for this stuff


