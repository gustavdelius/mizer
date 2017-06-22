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

load("Fmat.RData")
Fmat <- Fmat[,11:(dim(Fmat)[2]-1)]
load("Landings.RData")
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")

mypar <- c(11,log10(params_data$r_max))
mysd <- 1

# initialize by running from 1957 to 2010, to generate close to realistic initial conditions 
# myn to be used for later runs

parass <- mypar
capacity <- 10^(parass[1])
rmax <- 10^(parass[2:(1+length(parass))])
dd <- params_data
dd$r_max <- rmax[1:(length(dd$r_max))]
params <- MizerParams(dd, interaction = inter, kappa=capacity)
sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =1)
#plot(sim)
myn <- sim@n[dim(sim@n)[1],,]
sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=myn)
#plot(sim)

##### prepare landings data

L <- t(landings)
L[is.na(L)] <- 0
log_landings <- log10(10^(-10)+L[18:(dim(L)[1]-1),])
Y <- log_landings


################
model <- function(paras)
{
  capacity <- 10^(paras[1])
  rmax <- 10^(paras[2:(1+length(paras))])
  dd <- params_data
  dd$r_max <- rmax[1:(length(dd$r_max))]
  params <- MizerParams(dd, interaction = inter, kappa=capacity)
  sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=myn)
  return(log10((10^(-10))+getYield(sim)))
}


# make a test time series v using PDE model

v <- model(mypar)
names(dimnames(v))[names(dimnames(v))=="sp"] <- "Species"
ym <- melt(v)
ggplot(ym) + geom_line(aes(x=time,y=value, colour=
                             Species, linetype=Species)) + scale_y_continuous(
                               name="Log Yield") + scale_x_continuous(name="Time")
v-log_landings

# make function to evaluate model cost for a particular parameter choice

fish_model_cost <- function(par=mypar,YY=log_landings, sd=mysd){
  return(sum((model(par)-YY)^2/(sd^2)))
}
fish_model_cost()

#MCMC <- modMCMC(f = fish_model_cost, p = rep(9, length(mypar)),
#                niter = 10, jump = 0.1, updatecov = 500,
#                lower = rep(0, length(mypar)), upper = rep(50, length(mypar)))

MCMC <- modMCMC(f = fish_model_cost, p = rep(9, length(mypar)),
                niter = 100000, jump = 0.1, updatecov = 100,
                lower = rep(0, length(mypar)), upper = rep(50, length(mypar)))


MCMCrun1d <- MCMC
save(MCMCrun1d, file="MCMCrun1d.RData")

MCMCrun1$pars
#hist(MCMCrun1)

plot(MCMCrun1$pars[,1])
