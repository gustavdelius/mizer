s_params <- set_scaling_model()
sim <- project(s_params, t_max=5, effort = 0)
plotSpectra(sim)

gg<-getEGrowth(s_params, s_params@initial_n, s_params@initial_n_pp )
#ggnn <- gg*s_params@initial_n
#ggnn_back <- cbind(RDI,ggnn[,1:(dim(ggnn)[2]-1)])
#getZ
#parts for resource spectrum

length(diff(gg[4,]))
head(diff(gg[4,]))
head(gg[4,])

# remove last element, add -ve term on end
ggnn <- gg*s_params@initial_n
RR <- getRDD(s_params, s_params@initial_n, s_params@initial_n_pp)
FF <- ggnn
ZZ <- getZ(s_params, s_params@initial_n, s_params@initial_n_pp, effort = 0)
for ( i in (1:dim(ggnn)[1])){
X <- ggnn[i,] 
FF[i,] <- -(X - c(RR[i],X[1:(length(X)-1)]))/s_params@dw - ZZ[i,]*s_params@initial_n[i,]
}
# better to make the above using cbind to add in reproductive flux

# sort out F for the background resource too

#tmp <- (sim@params@rr_pp * sim@params@cc_pp / (sim@params@rr_pp + m2_background))
#n_pp <- tmp - (tmp - n_pp) * exp(-(sim@params@rr_pp+m2_background)*dt)

length(s_params@rr_pp)
m2_background <- getM2Background(s_params, n=s_params@initial_n, n_pp=s_params@initial_n_pp)

FF_resource <- s_params@rr_pp*(s_params@cc_pp-s_params@initial_n_pp)-m2_background*s_params@initial_n_pp

# make into a function of n and n_pp

# test does it go to zero at the steady state of my simulation

fish_rates <- function(params, n, n_pp){
    gg<-getEGrowth(s_params, n, n_pp )
    ggnn <- gg*n
    RR <- getRDD(params, n, n_pp)
    FF <- ggnn
    ZZ <- getZ(params, n, n_pp, effort = 0)
    for ( i in (1:dim(ggnn)[1])){
        X <- ggnn[i,] 
        FF[i,] <- -(X - c(RR[i],X[1:(length(X)-1)]))/params@dw - ZZ[i,]*n[i,]
    }
    return(FF)
}

resource_rates <- function(params, n, n_pp){
    m2_background <- getM2Background(params, n=n, n_pp=n_pp)
    
    FF_resource <- params@rr_pp*(params@cc_pp-n_pp)-m2_background*n_pp
    return(FF_resource)
}


my_n <- s_params@initial_n
my_n_pp <- s_params@initial_n_pp
FF <- fish_rates(s_params, my_n, my_n_pp)
FF_resource <- resource_rates(s_params, my_n, my_n_pp)
dataout <- c(FF_resource,as.vector(FF))
