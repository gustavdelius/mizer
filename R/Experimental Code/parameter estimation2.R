
#########

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
sim <- project(params, effort = 0, t_max = 10, dt = 0.1, t_save = .1)
#plot(sim)
DataY <- getBiomass(sim)
DataY

###

model <- function(paras)
{
  capacity <- 10^(paras[1])
  rmax <- 10^(paras[2])
  dd <- params_data[1,]
  dd$r_max <- rmax
  params <- MizerParams(dd, kappa=capacity)
  sim <- project(params, effort = 0, t_max = 30, dt = 0.1, t_save =0.1, no_w=1000)
  return(getBiomass(sim))
}

###log10
#h <- model(c(log(10^13), log(10^8)))

#h9 <- model(c(log(10^14), log(10^8)))



show_series <- function(h, isline=0){
  if (isline==0){
    return(plot(x=as.numeric(rownames(h)),y=h, xlab="time", ylab="biomass"))
  }else
    return(lines(x=as.numeric(rownames(h)),y=h, xlab="time", ylab="biomass", col="Red"))
}
h <- model(c(13,8))
h9 <- model(c(8,8))

show_series(h)
show_series(h9,isline=1)

# (carrying_capacity, r_max)
default_params <- c(13,8)
mysd <- 1

### generate training/empirical data using the time series from default params, plus 
#gaussian noise with standard deviation mysd 

length(model(c(13,8)))

u <- model(c(13,8))
DataY <- u+rnorm(length(u),mean = 0, sd=mysd)
show_series(DataY)

fish_model_cost <- function(par=default_params,Y=DataY, sd=mysd){
   return(sum((model(par)-Y)^2/(sd^2)))
}

#log_r_max =8
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
  
persp3D(x=log_carrying_capacity, y=log_r_max, z=log(costt_vals))


contour(x=log_carrying_capacity, y=log_r_max, z=log(costt_vals), xlab = "log_carrying_capacity",
ylab = "log_r_max")

costt_vals[8,13]
costt_vals[13,8]

lines(h9)
DataY
# output log(time series)
# normal errors
# interpolation
# modmcmc
