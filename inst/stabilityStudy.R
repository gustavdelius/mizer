library("BB")
params <- set_scaling_model(no_sp = 8)


sim <- project(params, t_max=15, effort = 0, t_save = 1)
plot(sim)


###################
library("rootSolve")

s_params <- params
used_params <- params


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


fish_ratesB <- function(params, nB, n_ppB){
  rate_n <- fish_rates(params, sweep(nB, 2, params@w^(params@lambda), "/"), n_ppB*(params@w_full^(-params@lambda))) 
  return(sweep(rate_n, 2, params@w^(params@lambda), "*"))
}  
resource_ratesB <- function(params, nB, n_ppB){
  rate_n_pp <-resource_rates(params, sweep(nB, 2, params@w^(params@lambda), "/"), n_ppB*(params@w_full^(-params@lambda)))
  return(rate_n_pp*(params@w_full^(params@lambda)))
}  



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


used_params <- s_params
change_functionB <- function(data_in){
  my_n_pp <- data_in[1:length(used_params@cc_pp)]
  #n <- data_in[(length(used_params@cc_pp)+1):length(data_in)]
  # matrix(vec,nrow = 7,ncol = 7)
  my_n <- array(data_in[(length(used_params@cc_pp)+1):length(data_in)], dim = dim(used_params@initial_n))
  
  FF <- fish_ratesB(used_params, my_n, my_n_pp)
  FF_resource <- resource_ratesB(used_params, my_n, my_n_pp)
  dataout <- c(FF_resource,as.vector(FF))
  return(dataout)
}


####################


output_ss <- c(params@initial_n_pp, params@initial_n)

mod <- function (t=0,y, parms=NULL,...) {return(change_function(y))}
jacobian_numeric <- jacobian.full(y = output_ss, func = mod)
EE <- eigen(jacobian_numeric)

real_parts <- sapply(EE$values, function(x) Re(x))
max(real_parts)


######################

nB <- sweep(params@initial_n, 2, used_params@w^(used_params@lambda), "*")*(used_params@initial_n>0)
n_ppB <- params@initial_n_pp*(used_params@w_full^(used_params@lambda))*(used_params@initial_n_pp>0)

################

output_ssB <- c(n_ppB, nB)

modB <- function (t=0,y, parms=NULL,...) {return(change_functionB(y))}
jacobian_numericB <- jacobian.full(y = output_ssB, func = modB)
EEB <- eigen(jacobian_numericB)

real_partsB <- sapply(EEB$values, function(x) Re(x))
max(real_partsB)

