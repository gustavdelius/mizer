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



sizeHistory <- function(param1, times=seq(0,1,0.1) ){
  test<-project(param1,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)
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


winf <- param1@species_params[4,2]
h <- param1@species_params[4,14]
gamma <- param1@species_params[4,15]

param1b <- param1
param1b@species_params[4,2] <- winf
param1b@species_params[4,14] <- h
param1b@species_params[4,15] <- gamma+10^(20)

param1b <- param1
XX <- param1b@species_params
XX[4,2] <- winf
XX[4,14] <- h
XX[4,15] <- gamma+10^(20)
param1b@species_params <- XX


mytimes <- seq(0,100,0.1)
tail(sizeHistory(param1,mytimes)[,4])


tail(sizeHistory(param1b,mytimes)[,4])


test<-project(param1,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)
plot(test)

plot(project(param1b,effort=0.5,t_max = 3, dt = 0.25, t_save = 1))


##############

params_data[4,7]

params_data2 <- params_data
params_data2[4,7] <- winf +500


##########

param1b <- param1
XX <- param1b@species_params
XX[4,2] <- winf*1
XX[4,14] <- h
XX[4,15] <- gamma*1
param1b@species_params <- XX

test1 <- project(param1,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)
test1b <- project(param1b,effort=0.5,t_max = 3, dt = 0.25, t_save = 1)

dim(test1@n)

plot(param1@w, test1@n[4,4,], log="xy")
points(param1@w, test1b@n[4,4,], col="Red")

sum(test1@n[4,4,]==test1b@n[4,4,])
#############

# only altering winf seems to make a difference

###########

plot(mytimes, sizeHistory(param1,mytimes)[,4])
points(mytimes, sizeHistory(param1b,mytimes)[,4], col="Red")

sizeHistory(param1b,mytimes)[,4]==sizeHistory(param1,mytimes)

tail(sizeHistory(param1,mytimes)[,4])

tail(sizeHistory(param1b,mytimes)[,4])

param1@search_vol

param1@activity

##############

params_data3 <- params_data
# set herring's max size to be larger
params_data3[4,2] <- 500
param3 <- MizerParams(params_data3, interaction = inter, no_w = 100)

plot(mytimes, sizeHistory(param1,mytimes)[,4])
points(mytimes, sizeHistory(param3,mytimes)[,4], col="Red")


#########################
param1@species_params[,3]

colnames(param1@species_params)

##################

params_data4 <- params_data
colnames(params_data4)

#rbind(params_data4,param1@species_params[,14])

head(param1@species_params)

param1@species_params[,14]

params_data4[["h"]] <- param1@species_params[,14]
params_data4[["gamma"]] <- param1@species_params[,15]*10^(10)

params_data4


# # # # 

param4 <- MizerParams(params_data4, interaction = inter, no_w = 100)

plot(mytimes, sizeHistory(param1,mytimes)[,4])
points(mytimes, sizeHistory(param4,mytimes)[,4], col="Red")

###############



