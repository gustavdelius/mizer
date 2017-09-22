library(mizer)
m_relative_effort <- sweep(f_history,2,colMeans((f_history)[19:29,]),"/")
f_history <- read.csv("./vignettes/NS_f_history.csv", row.names=1)
f_history <- as(f_history, "matrix")
# note that
#colMeans((f_history)[19:29,])
# is a permutation of 
#baseF <-
#  [1] 1.0399246 0.8685833 1.4873488 0.4000000 0.7665488 
#0.4000000 0.8171929 1.1308975 0.6881135 1.1435399
#[11] 0.9045996 0.9882556
# what is the connection with Fmat ?

load("Fmat.RData")
Fmat <- Fmat[,11:(dim(Fmat)[2]-1)]
load("Landings.RData")
params_data <- read.csv("./vignettes/NS_species_params.csv")
# should the line below be commented out ?
params_data$catchability <- as.numeric(colMeans((f_history)[19:29,]))

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
m_relative_effort <- sweep(f_history,2,colMeans((f_history)[19:29,]),"/")

#simini <- project(params, effort = t(Fmat), dt = 0.1, t_save =1)
#sim <- project(params, effort = t(Fmat), dt = 0.1, t_save =1, init_n=simini@n[dim(simini@n)[1],,],init_n_pp=simini@n_pp[dim(simini@n_pp)[1],])
simini <- project(params, effort = m_relative_effort, dt = 0.1, t_save =1)
sim <- project(params, effort = m_relative_effort, dt = 0.1, t_save =1, initial_n=simini@n[dim(simini@n)[1],,],initial_n_pp=simini@n_pp[dim(simini@n_pp)[1],])

vv <- log10((10^(-10))+getYield(sim))
j <- 6
plot(vv[,j])
plot(log_landings[,j])

m_relative_effort
Fmat
head(Fmat)
head(t(m_relative_effort))
# this is different to Fmat, but not by that much
