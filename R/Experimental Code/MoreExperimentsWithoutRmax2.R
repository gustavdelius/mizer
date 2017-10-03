library(mizer)
library(plyr)

epsi <- 0.1
# stable state is reached when we instead use
#epsi<-0.5
kappaR2 <- 10^(11)
lambda2 <- 2+0.9-(2/3)
source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
source("./R/Experimental Code/projectmodPREYSWITCH2.R")
params_data <- read.csv("./vignettes/NS_species_params.csv")
reduced_data <- params_data
reduced_data$r_max <- rep(10^(90),length(reduced_data$r_max))
params <- MizerParams(reduced_data)
sim <- project(params, effort = 1, t_max = 200, dt = 0.1, t_save = 1)
ini_n <- sim@n[dim(sim@n)[1],,]
ini_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]

sim2 <- project(params, effort = 1, t_max = 50, dt = 0.01, t_save = 1,initial_n=ini_n,initial_n_pp=ini_n_pp)
plot(sim2)

gb <- getBiomass(sim2)
dim(gb)

biomassform <- sweep(sim2@n,3,sim@params@w,"*")
plot(sim@params@w,biomassform[51,11,],log="xy")


library("rgl")
library("plot3D")
library("plot3Drgl")
persp3D(x=1:51,y=log(sim@params@w),z=log(biomassform[,11,]),ticktype = "detailed",theta = 140, phi = 30)

persp3d(x=1:51,y=log(sim@params@w),z=log(biomassform[,11,]),ticktype = "detailed",theta = 140, phi = 30,col = "lightblue")


persp(log(biomassform[,11,]))


image2D(log(biomassform[,11,]))


library(animation)
library(anim.plots)

anim.curve(x^t, times=10:50/10, n=20)


saveGIF({
  for(i in 1:100){
    curve(sin(x), from = -5 + (i * 0.05), to = 5 + (i * 0.05), col = "red", ylab = "")
    curve(cos(x), from = -5 + (i * 0.05), to = 5 + (i * 0.05), add = TRUE, col = "blue", ylab = "")
    legend("topright", legend = c("sin(x)", "cos(x)"), fill = c("red", "blue"), bty = "n")
  }
}, interval = 0.1, ani.width = 550, ani.height = 350)
