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

data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
# Randomize selectivity and catchability for proper test
params@catchability[] <- runif(prod(dim(params@catchability)), min=0, max=1)
params@selectivity[] <- runif(prod(dim(params@selectivity)), min=0, max=1)
no_sp <- dim(params@catchability)[2]
no_w <- length(params@w)
no_w_full <- length(params@w_full)
n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
n_full <- abs(rnorm(no_w_full))
# Two methods:
# Params + pred_rate
# Params + n + n_pp
pred_rate <- getPredRate(params,n,n_full)
m21 <- getM2(params,pred_rate=pred_rate)
m22 <- getM2(params,n,n_full)
# Test dims
#expect_equal(dim(m21), c(no_sp,no_w))
#expect_equal(dim(m21), c(no_sp,no_w))

#2 expect_equal(m21, m22)

plot(m21[1,])
lines(m22[1,])

