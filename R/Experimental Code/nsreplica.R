library(mizer)
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
#load("Fmat.RData")
#Fmat <- Fmat[,11:(dim(Fmat)[2]-1)]
load("Landings.RData")
#params_data$catchability <- as.numeric(colMeans((f_history)[19:29,]))
L <- t(landings)
L[is.na(L)] <- 0
log_landings <- log10(10^(-10)+L[18:(dim(L)[1]-1),])
#####
params_data$sel_func <- "sigmoid_length"
params_data$l25 <- c(7.6, 9.8, 8.7, 10.1, 11.5, 19.8, 16.4, 19.8, 11.5,
                     19.1, 13.2, 35.3)
params_data$l50 <- c(8.1, 11.8, 12.2, 20.8, 17.0, 29.0, 25.8, 29.0, 17.0,
                     24.3, 22.9, 43.6)
params_data$a <- c(0.007, 0.001, 0.009, 0.002, 0.010, 0.006, 0.008, 0.004,
                   0.007, 0.005, 0.005, 0.007)
params_data$b <- c(3.014, 3.320, 2.941, 3.429, 2.986, 3.080, 3.019, 3.198,
                   3.101, 3.160, 3.173, 3.075)
f_history <- read.csv("./vignettes/NS_f_history.csv", row.names=1)
f_history <- as(f_history, "matrix")
params_data$catchability <- as.numeric(f_history["1990",])
params <- MizerParams(params_data, inter, kappa = 9.27e10)
relative_effort <- sweep(f_history,2,f_history["1990",],"/")
relative_effort[as.character(1988:1992),]
initial_effort <- matrix(relative_effort[1,],byrow=TRUE, nrow=100,
                         ncol=ncol(relative_effort), dimnames = list(1867:1966))
relative_effort <- rbind(initial_effort,relative_effort)
sim <- project(params, effort=relative_effort, dt = 0.5, t_save = 1)