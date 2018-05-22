

fish_rates <- function(params, n, n_pp){
  gg<-getEGrowth(params, n, n_pp )
  ggnn <- gg*n
  #ggnn <- gg*n*params@w^(params@lambda)
  RR <- getRDD(params, n, n_pp)
  FF <- ggnn
  ZZ <- getZ(params, n, n_pp, effort = 0)
  for ( i in (1:dim(ggnn)[1])){
    X <- ggnn[i,] 
    if (i>1){
      X[params@species_params$w_min_idx[i]-1] <- RR[i]
      FF[i,] <- -(X - c(0,X[1:(length(X)-1)]))/params@dw - ZZ[i,]*n[i,]
    } else {
      FF[i,] <- -(X - c(RR[1],X[1:(length(X)-1)]))/params@dw - ZZ[i,]*n[i,]
    }
  }
  return(FF)
}

resource_rates <- function(params, n, n_pp){
  m2_background <- getM2Background(params, n=n, n_pp=n_pp)
  
  FF_resource <- params@rr_pp*(params@cc_pp-n_pp)-m2_background*n_pp
  return(FF_resource)
}

#
s_params <- set_scaling_model(no_sp = 4, no_w = 50)

nstart <- s_params@initial_n
nstart[1,] <- nstart[1,]

sim <- project(s_params, t_max=2, effort = 0, initial_n = nstart, t_save = 1)
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

#trial_data_in <- c(my_n_pp, as.vector(my_n))
trial_data_in <- c(s_params@initial_n_pp, s_params@initial_n)
ss_go <- multiroot(f = change_function, start = trial_data_in)

root_n_pp <- ss_go$root[1:length(used_params@cc_pp)]
root_n <- array(ss_go$root[(length(used_params@cc_pp)+1):length(ss_go$root)], dim = dim(used_params@initial_n))

sim <- project(used_params, t_max=15, effort = 0, initial_n = root_n, t_save = 1, initial_n_pp = root_n_pp)
plot(sim)

# numerically evaluate the Jacobian

mod <- function (t=0,y, parms=NULL,...) {return(change_function(y))}

jacobian_numeric <- jacobian.full(y = trial_data_in, func = mod)

# #74 In this code we have 4 scale invariant species. We use default stock recruitment 
# relationships. In this case we get newton raphson to hit on a steady state where all 
# 4 species coexist, but not near a power law. However, I think the state is weakly unstable, 
# and starting at our pre-computed initial condition allows newton raphson to find the steady state.
# However running size spectrum dynamics before inputting to newton raphson can lead 
# one to drift out of the newton-raphson basin of attraction of the coexistant steady state.

################## BB 

library(BB)

require("setRNG")
setRNG(list(kind="Wichmann-Hill", normal.kind="Box-Muller", seed=1236))

#ss <- multiroot(f = model, start = c(1, 1, 1))

ss_BB <- dfsane(par = c(1, 1, 1), fn = model, control = list(maxit=15000), method = 3)

ss_BB <- BBsolve(par = c(1, 1, 1), fn = model)


#ss_go_BB <- dfsane(par = trial_data_in, fn = change_function)
ss_go_BB <- BBsolve(par = trial_data_in, fn = change_function)


#ans <- dfsane(par=p0, fn=expo3)

root_n_pp_BB <- ss_go_BB$par[1:length(used_params@cc_pp)]
root_n_BB <- array(ss_go_BB$par[(length(used_params@cc_pp)+1):length(ss_go_BB$par)], dim = dim(used_params@initial_n))

sim <- project(used_params, t_max=15, effort = 0, initial_n = root_n_BB, t_save = 1, initial_n_pp = root_n_pp_BB)
plot(sim)

#

# got BB package working, but it wont seem to converge in this instance. Next I should go 
# back to Newton Raphson, and take out the do loop, other simplifications are: 
# only considering occupieable species-weight bins, and just considering aggregate abundancs, 
# and abundance at egg size

