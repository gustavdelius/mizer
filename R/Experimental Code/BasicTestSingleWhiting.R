source('./R/MizerParams-class.r')
source('./R/MizerSim-class.R')
source('./R/project_methods.R')
source('./R/selectivity_funcs.R')
source('./R/summary_methods.R')
source('./R/wrapper_functions.R')
source('./R/project_methods.R')
source('./R/plots.R')
source('./R/project.R')
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library("plot3D")
library("plot3Drgl")

params_data <- read.csv("./vignettes/NS_species_params.csv")
params_data$sel_func <- "sigmoid_length"
params_data$l25 <- c(7.6, 9.8, 8.7, 10.1, 11.5, 26.0, 16.4, 19.8, 11.5,
                                          19.1, 13.2, 35.3)
params_data$l50 <- c(8.1, 11.8, 12.2, 20.8, 17.0, 29.0, 25.8, 29.0, 17.0,
                     24.3, 22.9, 43.6)

params_data$a <- c(0.007, 0.001, 0.009, 0.002, 0.010, 0.006, 0.008, 0.004,
                   0.007, 0.005, 0.005, 0.007)
params_data$b <- c(3.014, 3.320, 2.941, 3.429, 2.986, 3.080, 3.019, 3.198,
                   3.101, 3.160, 3.173, 3.075)
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
#### sing
nme <- rownames(inter)[6]
inter <- as(inter[6,6], "matrix")
rownames(inter) <- nme
colnames(inter) <-nme
####
####
years_run <- 40
disirate <- 1
multip <- 1
intri <- 0

params <- MizerParams(params_data[6,], interaction = inter, no_w = 200, min_landing_weight = 250, disintegration_rate = disirate,
                      fraction_discarded = 1, predation_multiplier = multip, intrinsic_annihilation = intri)
sim <- project(params, effort = 1, t_max = years_run, dt = 0.1, t_save = 1)
nDisc <- sim@n[years_run+1,1,]
plot(sim)

params@w_full[1]

plot(x=params@w_full, y=sim@n_d[years_run,], type="b",log="xy",xlab = "Size (g)", ylab="Abundance Of Dead Fish")
plot(x=params@w, y=sim@n_d[years_run+1,], type="b",log="xy",xlab = "Size (g)", ylab="Abundance Of Dead Fish")
plot(x=params@w, y=sim@n[years_run,1,]*getFMort(sim@params, effort=1)[1,], type="b",log="xy", xlab = "Size (g)", ylab="Abundance Of New Discards")
plot(x=params@w, y=sim@n[years_run+1,1,]*getFMort(sim@params, effort=1)[1,], type="b",log="xy", xlab = "Size (g)", ylab="Abundance Of New Discards")

persp3D(x = (1:(years_run+1)), y = params@w, z = sim@n_d, xlab="Time", ylab="Size (g)", zlab="Abundance Of Dead Fish",
        lighting = TRUE) 
plotrgl(smooth = TRUE)
persp3D(x = (1:(years_run+1)), y = log(params@w), z = log(0.01+sim@n_d), xlab="Time", ylab="log(Size (g))", zlab="log(Abundance Of Dead Fish)",
        lighting = TRUE) 
plotrgl(smooth = TRUE)

#########
params_without_discarding <- MizerParams(params_data[6,], interaction = inter, no_w = 200, min_landing_weight = 0,
                                         disintegration_rate = disirate, fraction_discarded = 1, predation_multiplier = multip,
                                         intrinsic_annihilation = intri)
sim_without_discarding <- project(params_without_discarding, effort = 1, t_max = years_run, dt = 0.1, t_save = 1)
nNoDisc <- sim_without_discarding@n[years_run+1,1,]
# 
plot(sim_without_discarding)

plot(x=params@w, y=(nDisc-nNoDisc), type="b",log="xy",xlab = "Size (g)", ylab="Increase In Abundance From Discarding")

plot(x=params@w, y=nNoDisc, type="l",log="xy",xlab = "Size (g)", ylab="Whiting Abundance")
lines(x=params@w, y=nDisc, col ="red")

sum(nDisc - nNoDisc <1)

plot(x=params@w[(1:30)], y=sim@n_d[years_run,(1:30)], type="b",log="y",xlab = "Size (g)", ylab="log(Abundance Of Dead Fish)")

# plot change in abundance due to discarding
# plot change in biomass due to discarding
nDisc[1]

plot(x=params@w, y=log(nDisc/nNoDisc), type="b",log="x",xlab = "Size (g)", ylab="Increase In Abundance From Discarding")

