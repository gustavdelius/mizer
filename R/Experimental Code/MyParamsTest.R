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


no_sp <- 3
min_w <- 0.1
max_w <- 5000
no_w <- 200
min_w_pp <- 1e-8
no_w_pp <- 20
species_names <- c("Cod","Haddock","Whiting")
test_params <- MizerParams(no_sp, min_w = min_w, max_w = max_w, no_w = no_w, min_w_pp = min_w_pp, no_w_pp = no_w_pp, species_names = species_names)



# The two tests below are not appropriate for our new version of mizer with fft, so instead I added a new test for the new w_full
## expect_that(length(test_params@w_full), equals(no_w+no_w_pp))
## expect_that(length(test_params@dw_full), equals(no_w+no_w_pp))
# Check that that log of w_full is evenly spaced


max(diff(log(test_params@w_full))) 

min(diff(log(test_params@w_full)))

max(diff(log(test_params@w_full))) == min(diff(log(test_params@w_full)))

max(diff(log(test_params@w_full)))-min(diff(log(test_params@w_full)))<10^(-10)
