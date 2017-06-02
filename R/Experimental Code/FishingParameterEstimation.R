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

####
params <- MizerParams(params_data, interaction = inter, no_w = 1000)

ptm <- proc.time() 
sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)
Tsim <- (proc.time() - ptm)[[3]]

plot(sim)

# test different starting points

n<- sim@n[11, , ]
dim(n)

plot(n[1,],log="xy")


sim2 <- project(params, effort = 0.5, t_max = 100, dt = 0.1, t_save = 1, initial_n = n
)

plot(sim2)


head(n[1,])

n[1,] <- runif(length(n[1,]),0,1)

sim3 <- project(params, effort = 0.5, t_max = 100, dt = 0.1, t_save = 1, initial_n = n
)

plot(sim3)
