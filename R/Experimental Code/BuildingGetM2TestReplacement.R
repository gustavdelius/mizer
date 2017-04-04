


source('./R/MizerParams-classNewSlots.R')
#source('./R/MizerParams-class.r')
source('./R/MizerSim-class.R')
source('./R/project_methods.R')
source('./R/selectivity_funcs.R')
source('./R/summary_methods.R')
source('./R/wrapper_functions.R')
source('./R/plots.R')
source('./R/project.R')
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)

source('./R/project_methodsFFT2usingSlots.R')


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
    ###expect_equal(dim(m21), c(no_sp,no_w))
    ###expect_equal(dim(m21), c(no_sp,no_w))
    ## Replace this test with something appropriate for fft## expect_equal(m21, m22)
    
    ##############
    
    object <- params
    noSpecies <- dim(object@interaction)[1]
    muVals <- matrix(0, nrow = noSpecies, ncol = length(params@w))
    w <- params@w
    x <- log(w)
    x <- x - x[1]
    dx <- x[2]-x[1]
    feeding_level <- getFeedingLevel(object, n=n, n_pp=n_full)
    no_P <- length(object@smatM[1,])
    muIntermediate <- matrix(0, nrow = noSpecies, ncol = length(params@w))
    for (j in 1:noSpecies){
        f <- (1-feeding_level[j,])*object@search_vol[j,]*n[j,]*w
        f <- c(f[1:length(x)], rep(0, no_P-length(x)))
        mortalityIntegral <- dx*Re(fft((object@fsmatM[j,])*fft(f), inverse=TRUE)/no_P)
        muIntermediate[j, ] <- c(mortalityIntegral[(no_P-1):no_P], mortalityIntegral[1:(length(x)-1-1)])
    }
    for (i in 1:noSpecies){
        for (j in 1:noSpecies){
            muVals[i, ] <- muVals[i, ]+object@interaction[j,i]*muIntermediate[j, ]
        }
        
    }
    m22[1,]==muVals[1,]


