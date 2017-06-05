library(mizer)
library(FME)
library(ggplot2)
library(reshape2)
library(deSolve)
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library("plot3D")
library(rgl)
library("plot3Drgl")
library(optimx)
set.seed(123)
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")

model <- function(paras)
{
  capacity <- 10^(paras[1])
  rmax <- 10^(paras[2:(1+length(paras))])
  dd <- params_data
  dd$r_max <- rmax[1:(length(dd$r_max))]
  params <- MizerParams(dd, interaction = inter, kappa=capacity)
  sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save =1)
  return(log10(getYield(sim)))
}
#model(mypar)
#params <- MizerParams(dd, interaction = inter, kappa=10^(mypar[1]))
#sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
#plot(sim)


mypar <- c(11,log10(params_data$r_max))
mysd <- 1
XX <- model(mypar)
Y <-XX+rnorm(prod(dim(XX)),mean = 0,sd=mysd)

#somepar <- c(10,log10(params_data$r_max))
#sum((model(somepar)-Y)^2/(mysd^2))

# compute model cost (=-2loglikelihood)
fish_model_cost <- function(par=mypar,YY=Y, sd=mysd){
  return(sum((model(par)-YY)^2/(sd^2)))
}
fish_model_cost()

fish_model_res <- function(par = mypar[1:2]){
  usedpar <- mypar
  usedpar[1:2] <- par
  return(fish_model_cost(usedpar))
}
mypar


fish_model_res(par=c(11,12.868))


####################

dp <- 0.1
log_carrying_capacity <- seq(10,12,by=dp)
log_r_max <- seq(10,13,by=dp)
time_pts <- as.numeric(rownames(Y))

# make array to store time series data for many points
RD <- array(0,c(length(log_carrying_capacity),length(log_r_max),length(time_pts),dim(Y)[2]))
for (i in (1:length(log_carrying_capacity))){
  for (j in (1:length(log_r_max))){
    usedpar <- mypar
    usedpar[1:2] <- c(log_carrying_capacity[i],log_r_max[j])
    RD[i,j, ,] <- model(usedpar)
  }
}

time_series_to_cost <- function(timedata=Y,YY=Y, sd=mysd){
  return(sum((timedata-YY)^2/(sd^2)))
}
RDC <- array(0, c(length(log_carrying_capacity),length(log_r_max)))
for (i in (1:length(log_carrying_capacity))){
  for (j in (1:length(log_r_max))){
    RDC[i,j] <- time_series_to_cost(timedata=RD[i,j, ,],YY=Y, sd=mysd)
  }
}

########### save data ########
EvolutionData2 <- list(log_carrying_capacity, log_r_max, time_pts, RD, RDC)
names(EvolutionData2) <- c("log_carrying_capacity","log_r_max",
                           "time_pts", "RD", "RDC")
save(EvolutionData2, file="EvolutionData2.RData")
#load("EvolutionData2.RData")

########## plots ###########

persp3D(x=log_carrying_capacity, y=log_r_max, z=exp(-RDC/2),xlab = "log_carrying_capacity",
        ylab = "log_r_max", zlab = "likelihood")

plotrgl()

contour(x=log_carrying_capacity, y=log_r_max, z=exp(-RDC/2),xlab = "log_carrying_capacity",
        ylab = "log_r_max")

############# run MCMC #############

MCMC <- modMCMC(f = fish_model_res, p = c(12,7),
                niter = 100, jump = 0.1, updatecov = 10, lower=c(0.1,0.1),
                upper = c(20,20))
plot(MCMC)
summary(MCMC)
hist(MCMC)

########### reading part (uncomment and run to recover previously generated data) #####

# load("EvolutionData2.RData")
#log_carrying_capacity <- EvolutionData2$log_carrying_capacity
#log_r_max <- EvolutionData2$log_r_max
#RDC <- EvolutionData2$RDC

#persp3D(x=log_carrying_capacity, y=log_r_max, z=exp(-RDC/2),xlab = "log_carrying_capacity",
#        ylab = "log_r_max", zlab = "likelihood")

