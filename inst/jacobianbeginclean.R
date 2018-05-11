
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



s_params <- set_scaling_model()
sim <- project(s_params, t_max=5, effort = 0)
plot(sim)

################

my_n <- sim@n[dim(sim@n)[1],,]
my_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]
FF <- fish_rates(s_params, my_n, my_n_pp)
FF_resource <- resource_rates(s_params, my_n, my_n_pp)
dataout <- c(FF_resource,as.vector(FF))
max(abs(dataout))


