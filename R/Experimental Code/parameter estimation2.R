library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
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
dd <- params_data[1,]
dd$r_max <- 10^13
params <- MizerParams(dd, kappa=10^8)


# the output of model(c(log10_capacity, log10_rmax)) is a length 292 time series of 
# the log biomass of sprat over 30 years
model <- function(paras)
{
  capacity <- 10^(paras[1])
  rmax <- 10^(paras[2])
  dd <- params_data[1,]
  dd$r_max <- rmax
  params <- MizerParams(dd, kappa=capacity)
  sim <- project(params, effort = 0, t_max = 30, dt = 0.1, t_save =0.1, no_w=1000)
  return(log10(getBiomass(sim)))
}

# time series plotter
show_series <- function(h, isline=0){
  if (isline==0){
    return(plot(x=as.numeric(rownames(h)),y=h, xlab="time", ylab="log biomass"))
  }else
    return(lines(x=as.numeric(rownames(h)),y=h, xlab="time", ylab="log biomass", col="Red"))
}


show_series(model(c(13,8)))

# next we generate the empirical data
mypar <- c(13,8)
mysd <- 0.05
Y <- model(mypar) + rnorm(length(model(mypar)),mean = 0,sd=mysd)

# plot empirical data
show_series(Y)

# compute model cost (=-2loglikelihood)
fish_model_cost <- function(par=mypar,YY=Y, sd=mysd){
  return(sum((model(par)-YY)^2/(sd^2)))
}

fish_model_cost(par = c(8,8))
fish_model_cost(par = mypar)

####### plot model cost

log_carrying_capacity <- (1:15)
cost_vals <- sapply(log_carrying_capacity, function(cap) fish_model_cost(par=c(cap,8)))
plot(x=log_carrying_capacity, y=cost_vals)

log_r_max <- (1:15)
costt_vals <- matrix(0, nrow=length(log_carrying_capacity), ncol = length(log_r_max))
for (i in (1:length(log_carrying_capacity))){
  for (j in (1:length(log_r_max))){
    costt_vals[i,j] <- fish_model_cost(par=c(log_carrying_capacity[i],log_r_max[j]))
  }
}

# plots of how the model cost varies over parameter space

persp3D(x=log_carrying_capacity, y=log_r_max, z=costt_vals,xlab = "log_carrying_capacity",
        ylab = "log_r_max", zlab = "log_model_cost")

contour(x=log_carrying_capacity, y=log_r_max, z=costt_vals, xlab = "log_carrying_capacity",
        ylab = "log_r_max")


# A plot using the log of the model cost reveals that the model cost has a global minimum at 
# params = c(log_carrying_capacity, log_r_max) = mypar = c(13,8) 

persp3D(x=log_carrying_capacity, y=log_r_max, z=log(costt_vals),xlab = "log_carrying_capacity",
        ylab = "log_r_max", zlab = "log_model_cost")

contour(x=log_carrying_capacity, y=log_r_max, z=log(costt_vals),xlab = "log_carrying_capacity",
        ylab = "log_r_max")



