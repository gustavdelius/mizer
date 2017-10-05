library(mizer)
library(plyr)

epsi <- 0
kappaR2 <- 10^(11)
lambda2 <- 2+0.9-(2/3)
source("./R/Experimental Code/project_methodsmodPREYSWITCH2.R")
source("./R/Experimental Code/projectmodPREYSWITCH2.R")
params_data <- read.csv("./vignettes/NS_species_params.csv")
reduced_data <- params_data[11,]
reduced_data$r_max <- rep(10^(90),length(reduced_data$r_max))
#reduced_data$sigma == 1.3
reduced_data$sigma <- rep(1.3,length(reduced_data$r_max))

params <- MizerParams(reduced_data)
sim <- project(params, effort = 1, t_max = 200, dt = 0.1, t_save = 1)
ini_n <- sim@n[dim(sim@n)[1],,]
ini_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]

sim2 <- project(params, effort = 1, t_max = 50, dt = 0.01, t_save = 1,initial_n=ini_n,initial_n_pp=ini_n_pp)
plot(sim2)


for (i in (2:40)){
#sim3 <- project(params, effort = 1, t_max = i, dt = 0.01, t_save = 1,initial_n=ini_n,initial_n_pp=ini_n_pp)
#plot(sim3)
#jpeg('rplot.jpg')
#i
#gb <- getBiomass(sim3)
#gb[dim(gb)[1],]
  jpeg(paste0("./pics/h",i,".jpg"))
  plot(params@w,sim2@n[i,,]*params@w,log="xy",ylim=c(10^(4),10^(10)),type="l")
  

#plot(sim3)


#plot(x,y)
dev.off()
}

plot(params@w,sim3@n[dim(sim3@n)[1],,]*params@w,log="xy",ylim=c(10^(4),10^(10)),type="l")

v <- sim3@n[dim(sim3@n)[1],,]*params@w
params@w[which.max(v)]
#cod are starving before they reach max size ?
# biomass sloping up
# max is held fixed, that is the max size that it can feed off plankton


params@w[which.max(v)]/reduced_data$beta

paste0(1:12, c("st", "nd"))


gb <- getBiomass(sim2)
dim(gb)

biomassform <- sweep(sim2@n,3,sim@params@w,"*")
plot(sim@params@w,biomassform[51,1,],log="xy")


library("rgl")
library("plot3D")
library("plot3Drgl")
persp3D(x=1:251,y=log(sim@params@w),z=log(biomassform[,1,]),ticktype = "detailed",theta = 140, phi = 30)

persp3d(x=1:51,y=log(sim@params@w)[1:50],z=log(biomassform[1:51,1,1:50]+10^(-110)),ticktype = "detailed",theta = 140, phi = 30,col = "lightblue")


persp(log(biomassform[,1,1:50]))


image2D(log(biomassform[,1,1:70]+10^(-110)))

log(biomassform[,1,]+0.01)
    