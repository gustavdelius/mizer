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
Fmat2 <- Fmat
Fmat2[1,] <- Fmat[1,]
Fmat2[2,] <- Fmat[2,]
Fmat2[3,] <- Fmat[3,]
Fmat2[4,] <- Fmat[5,]
Fmat2[5,] <- Fmat[4,]
Fmat2[6,] <- Fmat[8,]
Fmat2[7,] <- Fmat[7,]
Fmat2[8,] <- Fmat[6,]
Fmat2[9,] <- Fmat[9,]
Fmat2[10,] <- Fmat[10,]
Fmat2[11,] <- Fmat[12,]
Fmat2[12,] <- Fmat[11,]
Fmat <- Fmat2

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
sim_ini <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=myn)
sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=sim_ini@n[dim(sim_ini@n)[1],,],initial_n_pp=sim_ini@n_pp[dim(sim_ini@n_pp)[1],])



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
  sim_ini <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=myn)
  sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, initial_n=sim_ini@n[dim(sim_ini@n)[1],,],initial_n_pp=sim_ini@n_pp[dim(sim_ini@n_pp)[1],])
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

optt <- optim(par=rep(9, length(mypar)),
             fish_model_cost,method="L-BFGS-B",lower=rep(0,13),upper = rep(50,13),control = list(maxit = 3))

optt <- optim(par=rep(9, length(mypar)),
              fish_model_cost,method="SANN",lower=rep(0,13),upper = rep(50,13),control = list(maxit = 3))

#MCMC <- modMCMC(f = fish_model_cost, p = rep(9, length(mypar)),
#                niter = 10, jump = 0.1, updatecov = 500,
#                lower = rep(0, length(mypar)), upper = rep(50, length(mypar)))

MCMC <- modMCMC(f = fish_model_cost, p = rep(9, length(mypar)),
                niter = 3000, jump = 0.1, updatecov = 100,
                lower = rep(0, length(mypar)), upper = rep(50, length(mypar)))


MCMCrun11 <- MCMC
save(MCMCrun11, file="MCMCrun11.RData")

MCMCrun1$pars
#hist(MCMCrun1)

plot(MCMCrun11$pars[,1])

mean_vec <- (1:dim(MCMCrun11$pars)[2])
sd_vec <- (1:dim(MCMCrun11$pars)[2])
for (i in (1:dim(MCMCrun11$pars)[2])){
mean_vec[i] <- mean(MCMCrun11$pars[,i])  
sd_vec[i] <- sd(MCMCrun11$pars[,i])  

}


