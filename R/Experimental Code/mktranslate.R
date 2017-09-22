load("for_richard.RData")

windup <- colMeans(output1[1000:dim(output1)[1],14:25])

colMeans(output1[1000:dim(output1)[1],1:12])

names(windup) <- c("Sprat",
                   "Sandeel",
                   "N.pout",
                   "Dab",
                   "Herring",
                   "Gurnard",
                   "Sole",
                   "Whiting",
                   "Plaice",
                   "Haddock",
                   "Saithe",
                   "Cod")


baseF <- c(1.0399246, 0.8685833, 1.4873488, 0.4000000,
           0.7665488, 0.4000000, 0.8171929, 1.1308975, 0.6881135, 1.1435399, 
           0.9045996, 0.9882556)
names(baseF) <- names(windup)

library(mizer)
load("Fmat.RData")
Fmat <- Fmat[,11:(dim(Fmat)[2]-1)]
load("Landings.RData")
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
L <- t(landings)
L[is.na(L)] <- 0
log_landings <- log10(10^(-10)+L[18:(dim(L)[1]-1),])
Y <- log_landings
# values copied from table 3 of mike's paper
capacity <- exp(25.210)
rmax <- exp(c(26.7, 26, 30.67,26.56,23.1,26.03, 22.95, 25.38,30.56, 28.375, 22.77,26.92))
dd <- params_data
dd$r_max <- rmax[1:(length(dd$r_max))]
params <- MizerParams(dd, interaction = inter, kappa=capacity)
simini <- project(params, effort = windup*baseF, dt = 0.1, t_save =1, t_max=300)
sim <- project(params, effort = windup*baseF, dt = 0.1, t_save =1, initial_n=simini@n[dim(simini@n)[1],,],initial_n_pp=simini@n_pp[dim(simini@n_pp)[1],])
simini <- project(params, effort = sweep(t(Fmat),2,baseF,"/"), dt = 0.1, t_save =1, t_max=300)
sim <- project(params, effort = sweep(t(Fmat),2,baseF,"/"), dt = 0.1, t_save =1, initial_n=simini@n[dim(simini@n)[1],,],initial_n_pp=simini@n_pp[dim(simini@n_pp)[1],])

vv <- log((10^(-10))+getYield(sim))
j <- 6
plot(vv[,j])
plot(log(10^log_landings[,j]))

params@species_params$species[j]
params

#######################

f_history <- read.csv("./vignettes/NS_f_history.csv", row.names=1)
f_history <- as(f_history, "matrix")
colMeans((f_history)[19:29,])
baseF <-
  [1] 1.0399246 0.8685833 1.4873488 0.4000000 0.7665488 
0.4000000 0.8171929 1.1308975 0.6881135 1.1435399
[11] 0.9045996 0.9882556

# sweep
f_history/colMeans((f_history)[19:29,])

#relative_effort <- sweep(f_history,2,f_history["1990",],"/")
#params_data$catchability <- as.numeric(f_history["1990",])

m_relative_effort <- sweep(f_history,2,colMeans((f_history)[19:29,]),"/")
# compare this with Fmat

params_data$catchability <- as.numeric(colMeans((f_history)[19:29,]))
# add to short parest
# compare with Fmat and landings data
# go through NS model (code sigmoids etc.) and try again