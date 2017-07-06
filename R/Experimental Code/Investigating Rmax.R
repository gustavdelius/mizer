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

###

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


###

paras <- mypar
paras[2]
capacity <- 10^(paras[1])
rmax <- 10^(paras[2:(1+length(paras))])
dd <- params_data
dd$r_max <- rmax[1:(length(dd$r_max))]
params <- MizerParams(dd, interaction = inter, kappa=capacity)
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save =1)

loglandings <- log10(getYield(sim))

loglandings[10,1]

plot(sim)

dim(sim@n)

plot(params@w,sim@n[11,1,],log="xy")

#####
R_pts <- (1:20)
loglandingsout <- R_pts
for (g in R_pts){
    paras <- mypar
    paras[2] <- g
    capacity <- 10^(paras[1])
    rmax <- 10^(paras[2:(1+length(paras))])
    dd <- params_data
    dd$r_max <- rmax[1:(length(dd$r_max))]
    params <- MizerParams(dd, interaction = inter, kappa=capacity)
    sim <- project(params, effort = 1, t_max = 30, dt = 0.1, t_save =1)
    myn <- sim@n[dim(sim@n)[1],,]
    sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save =1,initial_n=myn)
    loglandingsout[g] <- log10(getYield(sim))[10,1]
}

################

R_pts_c <- R_pts[1:10]
loglandingsout_c <- loglandingsout[1:10]
mylm <- lm(loglandingsout_c~R_pts_c)
summary(mylm)
#####

plot(R_pts_c,loglandingsout_c)
lines(R_pts_c,coefficients(mylm)[1]+coefficients(mylm)[2]*R_pts_c)

coefficients(mylm)

#########

R_pts <- (1:20)
loglandingsout <- R_pts
ndata <- matrix(0, nrow=length(R_pts),ncol=100)
for (g in R_pts){
    paras <- mypar
    paras[2] <- g
    capacity <- 10^(paras[1])
    rmax <- 10^(paras[2:(1+length(paras))])
    dd <- params_data
    dd$r_max <- rmax[1:(length(dd$r_max))]
    params <- MizerParams(dd, interaction = inter, kappa=capacity)
    sim <- project(params, effort = 1, t_max = 30, dt = 0.1, t_save =1)
    myn <- sim@n[dim(sim@n)[1],,]
    sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save =1,initial_n=myn)
    loglandingsout[g] <- log10(getYield(sim))[10,1]
    myn2 <- sim@n[dim(sim@n)[1],,]
    ndata[g,] <- myn2[1,]
}


plot(log(params@w),log(ndata[15,]),ylim=c(1,40))
points(log(params@w),log(ndata[5,]),col="RED")
points(log(params@w),log(ndata[10,]),col="BLUE")

plot(R_pts,log(ndata[,1]))

