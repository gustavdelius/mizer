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

#> reduced_data$w_mat
#[1] 1606

#reduced_data$w_mat <- 800
#reduced_data$w_mat <- 3800
#reduced_data$sigma == 1.3
reduced_data$sigma <- rep(1.3,length(reduced_data$r_max))

params <- MizerParams(reduced_data)
sim <- project(params, effort = 1, t_max = 200, dt = 0.1, t_save = 1)
ini_n <- sim@n[dim(sim@n)[1],,]
ini_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]
max_grower <- 881.6151
ini_n2 <- rep(0,length(ini_n))
ini_n2[params@w<=max_grower] <- ini_n[params@w<=max_grower]
sim2 <- project(params, effort = 1, t_max = 50, dt = 0.01, t_save = 1,initial_n=ini_n2,initial_n_pp=ini_n_pp)
plot(sim2)


for (i in (2:40)){
#sim3 <- project(params, effort = 1, t_max = i, dt = 0.01, t_save = 1,initial_n=ini_n,initial_n_pp=ini_n_pp)
#plot(sim3)
#jpeg('rplot.jpg')
#i
#gb <- getBiomass(sim3)
#gb[dim(gb)[1],]
  jpeg(paste0("./pics/h",i,".jpg"))
  plot(params@w,sim2@n[i,,]*params@w,log="xy",ylim=c(10^(-4),10^(14)),type="l")
  points(params@w,params@w*sim2@n_pp[i,params@w_full>=params@w[1]])
  abline(v=reduced_data$w_mat)
  lines(params@w,params@psi,col="Red")
  lines(params@w,getEGrowth(params,t(as.matrix(sim2@n[i,,])),sim2@n_pp[i,]),col="Blue")
  
  

#plot(sim3)


#plot(x,y)
dev.off()
}

plot(params@w,sim3@n[dim(sim3@n)[1],,]*params@w,log="xy",ylim=c(10^(4),10^(10)),type="l")

plot(params@w,sim2@n[i,,]*params@w,log="xy",ylim=c(10^(4),10^(14)),type="l")
points(params@w,params@w*sim2@n_pp[i,params@w_full>=params@w[1]])
abline(v=reduced_data$w_mat)

##################################

i <- 22
plot(params@w,sim2@n[i,,]*params@w,log="xy",ylim=c(10^(-8),10^(14)),type="l")
points(params@w,params@w*sim2@n_pp[i,params@w_full>=params@w[1]])
abline(v=reduced_data$w_mat)
lines(params@w,params@psi,col="Red")
lines(params@w,getEGrowth(params,t(as.matrix(sim2@n[i,,])),sim2@n_pp[i,]),col="Blue")
#getEGrowth(params,as.matrix(sim2@n[i,,,drop=FALSE]),sim2@n_pp[i,])
gg <- getEGrowth(params,t(as.matrix(sim2@n[i,,])),sim2@n_pp[i,])
params@w[length(gg[1,gg[1,]>0])+1]
#617.8118

# why is there biomass above the max growth rate ?
# RUN SAME CODE IN ORIG MIZER
# INVESTIGATE INITIAL CNDN DEPENDENCE


#####################################

v <- sim2@n[dim(sim2@n)[1],,]*params@w
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
    
###############################

plot(params@w,params@psi,col="Red",log="x")
abline(v=params@species_params$w_mat)

##########################
gg <- getEGrowth(params,t(as.matrix(sim2@n[i,,])),sim2@n_pp[i,])

###############

i <- 22
n <- t(as.matrix(sim2@n[i,,]))
n_pp <- sim2@n_pp[i,]
# Calculate amount E_{a,i}(w) of available food
phi_prey <- getPhiPrey(sim@params, n=n, n_pp=n_pp)
# Calculate amount f_i(w) of food consumed
feeding_level <- getFeedingLevel(sim@params, n=n, n_pp=n_pp, phi_prey=phi_prey)
# Calculate the resources available for reproduction and growth
e <- getEReproAndGrowth(sim@params, n=n, n_pp=n_pp, feeding_level=feeding_level)
# Calculate the resources for reproduction
e_spawning <- getESpawning(sim@params, n=n, n_pp=n_pp, e=e)
# Calculate the growth rate g_i(w)
e_growth <- getEGrowth(sim@params, n=n, n_pp=n_pp, e_spawning=e_spawning, e=e)
plot(params@w,sim2@n[i,,]*params@w,log="xy",ylim=c(10^(-8),10^(14)),type="l")
points(params@w,params@w*sim2@n_pp[i,params@w_full>=params@w[1]])
abline(v=reduced_data$w_mat)
lines(params@w,e_growth,col="Blue")
lines(params@w,e_spawning,col="Green")

max_grower <- params@w[length(gg[1,gg[1,]>0])+1]
#881.6151
max_grower <- 881.6151
ini_n2 <- rep(0,length(ini_n))
ini_n2[params@w<=max_grower] <- ini_n[params@w<=max_grower]