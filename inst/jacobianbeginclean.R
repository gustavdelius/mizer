
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

#

s_params <- set_scaling_model()

nstart <- s_params@initial_n
nstart[1,] <- nstart[1,]

sim <- project(s_params, t_max=15, effort = 0, initial_n = nstart, t_save = 1)
plot(sim)

###########################

library("rootSolve")

model <- function(x) {
   F1 <- x[1] + x[2] + x[3]^2 -12
   F2 <- x[1]^2 - x[2] + x[3] -2
   F3 <- 2*x[1] - x[2]^2 + x[3] -1
   c(F1 = F1, F2 = F2, F3 = F3)
   }


ss <- multiroot(f = model, start = c(1, 1, 1))

ss$root

ss$f.root

################

results <- 1:dim(sim@n)[1]

for (T in (1:dim(sim@n)[1])){
    
    my_n <- sim@n[T,,]
    #my_n[1, ] <- 1000*my_n[1, ]
    my_n_pp <- sim@n_pp[T,]
    FF <- fish_rates(s_params, my_n, my_n_pp)
    FF_resource <- resource_rates(s_params, my_n, my_n_pp)
    dataout <- c(FF_resource,as.vector(FF))
    results[T] <- max(abs(dataout))
}
plot(results,type="l")
########################

used_params <- s_params
change_function <- function(data_in){
  my_n_pp <- data_in[1:length(used_params@cc_pp)]
  #n <- data_in[(length(used_params@cc_pp)+1):length(data_in)]
 # matrix(vec,nrow = 7,ncol = 7)
  my_n <- array(data_in[(length(used_params@cc_pp)+1):length(data_in)], dim = dim(used_params@initial_n))
  
FF <- fish_rates(used_params, my_n, my_n_pp)
FF_resource <- resource_rates(used_params, my_n, my_n_pp)
dataout <- c(FF_resource,as.vector(FF))
return(dataout)
}

##########################

trial_data_in <- c(my_n_pp, as.vector(my_n))

ss_go <- multiroot(f = change_function, start = trial_data_in)

