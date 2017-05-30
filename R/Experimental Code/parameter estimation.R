
#########

library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library(mizer)

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
  capacity <- exp(paras[1])
  rmax <- exp(paras[2])
  dd <- params_data[1,]
  dd$r_max <- rmax
  params <- MizerParams(dd, kappa=capacity)
  sim <- project(params, effort = 0, t_max = 30, dt = 0.1, t_save =0.1)
  return(getBiomass(sim))
}

###
h <- model(c(log(10^13), log(10^8)))

h9 <- model(c(log(10^14), log(10^8)))


plot(h)
lines(h9)
DataY
# output log(time series)
# normal errors
# interpolation
# modmcmc
