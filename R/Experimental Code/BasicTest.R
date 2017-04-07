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


params_data <- read.csv("./vignettes/NS_species_params.csv")

params_data$sel_func <- "sigmoid_length"
params_data$l25 <- c(7.6, 9.8, 8.7, 10.1, 11.5, 19.8, 16.4, 19.8, 11.5,
                     19.1, 13.2, 35.3)
params_data$l50 <- c(8.1, 11.8, 12.2, 20.8, 17.0, 29.0, 25.8, 29.0, 17.0,
                     24.3, 22.9, 43.6)
params_data$a <- c(0.007, 0.001, 0.009, 0.002, 0.010, 0.006, 0.008, 0.004,
                   0.007, 0.005, 0.005, 0.007)
params_data$b <- c(3.014, 3.320, 2.941, 3.429, 2.986, 3.080, 3.019, 3.198,
                   3.101, 3.160, 3.173, 3.075)

inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
params <- MizerParams(params_data, interaction = inter, no_w = 200, min_landing_weight = 1100, disintegration_rate = 0.2)

params@w[1]

sim <- project(params, effort = 1, t_max = 10, dt = 0.1, t_save = 1)

plot(sim)

plot(x=params@w_full, y=sim@n_d[11,], type="b",log="xy")

plot(x=params@w_full, y=sim@n_d[10,], type="b",log="xy")

plot(x=params@w, y=sim@n[11,2,], type="b",log="xy")

sum(sim@n_d[11,]==0)

plot(sim@n_d[1,])

ggg[100:120]

params_data

params@discard_fraction[1,]

params@disintegration


#####################

no_sp <- 3
min_w <- 0.1
max_w <- 5000
no_w <- 200
min_w_pp <- 1e-8
no_w_pp <- 20
species_names <- c("Cod","Haddock","Whiting")
test_params <- MizerParams(no_sp, min_w = min_w, max_w = max_w, no_w = no_w, min_w_pp = min_w_pp, no_w_pp = no_w_pp, species_names = species_names)

length(test_params@w_full)
no_w+no_w_pp
no_w
no_w_pp

max(diff(log(test_params@w_full)))

min(diff(log(test_params@w_full)))


test_params@w_full[1+length(test_params@w_full)-length(test_params@w)]

test_params@w[1]


#expect_that(dim(test_params@pred_kernel), equals(c(no_sp,no_w,no_w+no_w_pp)))

dim(test_params@pred_kernel)

c(no_sp,no_w,no_w+no_w_pp)

c(no_sp,no_w,length(test_params@w_full))

#expect_that(dim(test_params@pred_kernel), equals(c(no_sp,no_w,length(test_params@w_full))))


##############

data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
sim <- project(params, effort=1, t_max=40, dt = 1, t_save = 1)
# Just check that it doesn't crash
plot(sim)
#expect_that(plot(sim),!throws_error())    
