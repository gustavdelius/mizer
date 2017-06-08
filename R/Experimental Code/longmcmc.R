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

mypar <- c(11,log10(params_data$r_max))
mysd <- 1
XX <- model(mypar)
Y <- XX+rnorm(prod(dim(XX)),mean = 0,sd=mysd)

show_time_series <- function(dataa =Y){
  y <- dataa
  names(dimnames(y))[names(dimnames(y))=="sp"] <- "Species"
  ym <- melt(y)
  #return(ggplot(ym) + geom_line(aes(x=time,y=value, colour=Species, linetype=Species)) + scale_y_continuous(trans="log10", name="Log Yield") + scale_x_continuous(name="Time"))
  return(ggplot(ym) + geom_line(aes(x=time,y=value, colour=Species, linetype=Species)) + scale_y_continuous( name="Log Yield") + scale_x_continuous(name="Time"))
  #return(ggplot(ym) + geom_point(aes(x=time,y=value, colour=Species, pointtype=Species)) + scale_y_continuous( name="Log Yield") + scale_x_continuous(name="Time"))
}
#show_time_series(Y)

fish_model_cost <- function(par=mypar,YY=Y, sd=mysd){
  return(sum((model(par)-YY)^2/(sd^2)))
}
fish_model_cost()

#myOP<-optim(fn = fish_model_cost, par=mypar,YY=Y, sd=mysd, method="BFGS")

MCMC <- modMCMC(f = fish_model_cost, p = rep(9, length(mypar)),
                niter = 10000, jump = 0.1, updatecov = 500)

MCMC$pars[3,]

MCMCplay <- MCMC
save(MCMCplay, file="MCMCplay.RData")

#load("MCMCplay.RData")
#plot(MCMCplay$pars[,1])
#summary(MCMCplay)
