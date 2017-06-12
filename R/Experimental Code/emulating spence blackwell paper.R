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

y <- model(mypar)

t(Fmat)


dim(landings)


landings[1,]


#####
dataa <- landings
y <- dataa
names(dimnames(y))[names(dimnames(y))=="sp"] <- "Species"
ym <- melt(y)
#return(ggplot(ym) + geom_line(aes(x=time,y=value, colour=Species, linetype=Species))
#+ scale_y_continuous(trans="log10", name="Log Yield") +
#scale_x_continuous(name="Time"))
return(ggplot(ym) + geom_line(aes(x=time,y=value, colour=
                                    Species, linetype=Species)) + scale_y_continuous(
  name="Log Yield") + scale_x_continuous(name="Time"))
ggplot(ym)

ggplot(ym) + geom_line(aes(x=time,y=value, colour=
                             Species, linetype=Species)) + scale_y_continuous(
                               name="Log Yield") + scale_x_continuous(name="Time")

YL <- landings
dim(YL)
dim(y)
t(Fmat)

#### cut off matrix to start from 1967
