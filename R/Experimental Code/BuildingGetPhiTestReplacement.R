


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

data(NS_species_params_gears)
data(inter)
params <- MizerParams(NS_species_params_gears, inter)
no_sp <- dim(params@catchability)[2]
no_w <- length(params@w)
no_w_full <- length(params@w_full)
n <- abs(array(rnorm(no_w * no_sp), dim = c(no_sp, no_w)))
n_full <- abs(rnorm(no_w_full))
pp <- getPhiPrey(params,n,n_full)
# test dim
#expect_that(dim(pp), equals(c(no_sp,no_w)))
# Test numbers are right
# Hideous - doing it by hand for first predator - should be the same
n <- abs(matrix(rnorm(no_sp*no_w),nrow = no_sp,ncol = no_w))
n_pp <- abs(rnorm(no_w_full))
neff1 <- rep(0,length(params@w))
for (i in 1:no_w)
    neff1[i] <- neff1[i] + sum(params@interaction[1,] * n[,i] * params@w[i] * params@dw[i])
w_offset <- no_w_full - no_w
pks <- rep(NA,length(params@w))
pkpp <- rep(NA,length(params@w))
for (i in 1:length(params@w)){
    pks[i] <- sum(params@pred_kernel[1,i,(w_offset+1):no_w_full] * neff1)
    pkpp[i] <- sum(n_pp * params@w_full * params@dw_full * params@pred_kernel[1,i,])
}
pp1 <- pks + pkpp


ppp <- getPhiPrey(params,n,n_pp)
pp <- getPhiPrey(params,n,n_pp)

max(abs(pp1-ppp[1,]))<10^(-8)





###########################
ppp <- getPhiPrey(params,n,n_pp)

object <- params


###########################


w0 <- object@w[1]
Beta <- log(object@species_params$beta)
sigma <- object@species_params$sigma
wFull <- object@w_full
xFull <- log(wFull)
xFull <- xFull - xFull[1]
dx <- xFull[2]-xFull[1]
idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
fullEnergyMat <- matrix(0, nrow = dim(n)[1], ncol=length(object@w))
fishEaten <- rep(0, length.out = length(wFull))
for(i in 1:(dim(n)[1])) {
    fishEaten[idx_sp] <- (object@interaction %*% n)[i, ]
    f2 <- wFull*wFull*(n_pp + fishEaten)
    fullEnergy <- dx*Re(fft(object@fsmat[i, ]*fft(f2), inverse=TRUE)/length(object@smat[i,]))
    fullEnergyMat[i,] <- fullEnergy[idx_sp]
}

##############################

source('./R/project_methodsFFT2usingSlots.R')
pppfft <- getPhiPrey(params,n,n_pp)

pppfft[1,]==fullEnergyMat[1,] 

all.equal(pppfft,fullEnergyMat)

