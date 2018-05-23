library("rootSolve")


fish_rates <- function(params, n, n_pp){
  gg<-getEGrowth(params, n, n_pp )
  ggnn <- gg*n
  #ggnn <- gg*n*params@w^(params@lambda)
  RR <- getRDD(params, n, n_pp)
  
  shift <- cbind(rep(0, dim(ggnn)[1]), ggnn[,1:(dim(ggnn)[2]-1)])
  
  #for (i in (1:dim(ggnn)[1])){
  #  shift[i, params@species_params$w_min_idx[i]] <- RR[i]  
  #}
  
  no_sp <- dim(ggnn)[1]
  indx <- (params@species_params$w_min_idx-1) * no_sp + (1:no_sp)
  shift[indx] <- RR
  
  
  
  FF <- sweep(-(ggnn - shift), 2, params@dw, "/") - n*getZ(params, n, n_pp, effort = 0)
  return(FF)
}


resource_rates <- function(params, n, n_pp){
  m2_background <- getM2Background(params, n=n, n_pp=n_pp)
  
  FF_resource <- params@rr_pp*(params@cc_pp-n_pp)-m2_background*n_pp
  return(FF_resource)
}

#
#s_params <- set_scaling_model(no_sp = 4, no_w = 50)
s_params <- set_scaling_model(no_sp = 5)


nstart <- s_params@initial_n
nstart[1,] <- nstart[1,]

#sim <- project(s_params, t_max=2, effort = 0, initial_n = nstart, t_save = 1)
#plot(sim)

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

