#biomass2density_fish <- function(params, n)
#density2biomass

# explore two species case more carefully

#fish_rates_B(nB,)

#fish_ratesB <- function(params, nB, n_ppB){
#  return(fish_rates(params, sweep(nB, 2, params@w^(params@lambda), "/"), n_ppB*(params@w_full^(-params@lambda))))
#}  


################################

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


fish_ratesB <- function(params, nB, n_ppB){
  rate_n <- fish_rates(params, sweep(nB, 2, params@w^(params@lambda), "/"), n_ppB*(params@w_full^(-params@lambda))) 
  return(sweep(rate_n, 2, params@w^(params@lambda), "*"))
}  
resource_ratesB <- function(params, nB, n_ppB){
  rate_n_pp <-resource_rates(params, sweep(nB, 2, params@w^(params@lambda), "/"), n_ppB*(params@w_full^(-params@lambda)))
  return(rate_n_pp*(params@w_full^(params@lambda)))
}  


#
#s_params <- set_scaling_model(no_sp = 2, no_w = 100)
s_params <- set_scaling_model(no_sp = 4, no_w = 100)


sim <- project(s_params, t_max=15, effort = 0, initial_n = s_params@initial_n, t_save = 1, initial_n_pp = s_params@initial_n_pp)
plot(sim)


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


##########################

#trial_data_in <- c(my_n_pp, as.vector(my_n))


trial_data_in <- c(s_params@initial_n_pp, s_params@initial_n)

trial_data_inB <- c(s_params@initial_n_pp*(s_params@w_full^(s_params@lambda)), 
                    sweep(s_params@initial_n, 2, s_params@w^(s_params@lambda), "*"))


#ss_go <- multiroot(f = change_function, start = trial_data_in)
#root_n_pp <- ss_go$root[1:length(used_params@cc_pp)]
#root_n <- array(ss_go$root[(length(used_params@cc_pp)+1):length(ss_go$root)], dim = dim(used_params@initial_n))

ss_goB <- multiroot(f = change_functionB, start = trial_data_inB, rtol = 0, atol = 1e-15, ctol = 0)


root_n_ppB <- ss_goB$root[1:length(used_params@cc_pp)]
root_nB <- array(ss_goB$root[(length(used_params@cc_pp)+1):length(ss_goB$root)], dim = dim(used_params@initial_n))

#root_n <- sweep(root_nB, 2, used_params@w^(used_params@lambda), "/")
#root_n_pp <- root_n_ppB*(used_params@w_full^(-used_params@lambda))

root_n <- sweep(root_nB, 2, used_params@w^(used_params@lambda), "/")*(used_params@initial_n>0)
root_n_pp <- root_n_ppB*(used_params@w_full^(-used_params@lambda))*(used_params@initial_n_pp>0)


sim <- project(used_params, t_max=15, effort = 0, initial_n = root_n, t_save = 1, initial_n_pp = root_n_pp)
plot(sim)
########################## other stuff from 2speciesattractive.R removed below


output_ss <- c(root_n_pp, root_n)

mod <- function (t=0,y, parms=NULL,...) {return(change_function(y))}
jacobian_numeric <- jacobian.full(y = output_ss, func = mod)
EE <- eigen(jacobian_numeric)

real_parts <- sapply(EE$values, function(x) Re(x))
max(real_parts)

######## why does the jacobian (biomass) code below break ?

output_ssB <- c(root_n_ppB, root_nB)

modB <- function (t=0,y, parms=NULL,...) {return(change_functionB(y))}
jacobian_numericB <- jacobian.full(y = output_ssB, func = modB)
EEB <- eigen(jacobian_numericB)

real_partsB <- sapply(EEB$values, function(x) Re(x))
max(real_partsB)

plot(used_params@w,root_n[2,],log="xy")

################

run2fixed <- function(nI,n_ppI){
  
  
  sim <- project(used_params, t_max=150, effort = 0, initial_n = nI, t_save = 1, initial_n_pp = n_ppI)
  n <- sim@n[dim(sim@n)[1],,]
  
  n_pp <- sim@n_pp[dim(sim@n_pp)[1],]
  
  
  trial_data_inB <- c(n_pp*(s_params@w_full^(s_params@lambda)), 
                      sweep(n, 2, s_params@w^(s_params@lambda), "*"))
  
  ss_goB <- multiroot(f = change_functionB, start = trial_data_inB, rtol = 0, atol = 1e-15, ctol = 0)
  
  
  
  root_n_ppB <- ss_goB$root[1:length(used_params@cc_pp)]
  root_nB <- array(ss_goB$root[(length(used_params@cc_pp)+1):length(ss_goB$root)], dim = dim(used_params@initial_n))
  
  
  
  root_n <- sweep(root_nB, 2, used_params@w^(used_params@lambda), "/")*(used_params@initial_n>0)
  root_n_pp <- root_n_ppB*(used_params@w_full^(-used_params@lambda))*(used_params@initial_n_pp>0)
  
  return(list(n=root_n,n_pp=root_n_pp))
  
}

RF <- run2fixed(s_params@initial_n,s_params@initial_n_pp)

plot(s_params@w,RF$n[1,],log="xy", type="l",ylim=c(min(RF$n[RF$n>0]),max(RF$n)))
lines(s_params@w,RF$n[2,])

random_ini <- function(L){
  nn <- s_params@initial_n
  for (i in 1:dim(nn)[1]){
    nn[i,] <- nn[i,]*10^(runif(1,-L,L))
  }
  return(nn)
}

#n_start <- random_ini(3)

# another steady state where species 2 is not v abundant
# rr <- c(0.9552858, -2.4299018)

for (i in 1:20){
  rr <- runif(dim(s_params@initial_n)[1],-3,3)
  print(rr)
  n_start <- s_params@initial_n
  for (i in 1:dim(n_start)[1]){
    n_start[i,] <- n_start[i,]*10^(rr[i])
  }
  RF <- run2fixed(n_start,s_params@initial_n_pp)
  plot(s_params@w,RF$n[1,],log="xy", type="l",ylim=c(min(RF$n[RF$n>0]),max(RF$n)))
  lines(s_params@w,RF$n[2,])
}

# another steady state N_1 (w) *10^(0.9552858)
# N_2(w) *10^(-2.4299018)

n_start <- s_params@initial_n
for (i in 1:dim(n_start)[1]){
  X <- runif(dim(n_start)[2],0 ,10^(-7))
  H <- n_start[i,]
  H[H>0] <- X[H>0]
  n_start[i,] <- H
  #  n_start[i,] <- runif(0,10^6,dim(n_start)[2])
}
RF <- run2fixed(n_start,s_params@initial_n_pp)
plot(s_params@w,RF$n[1,],log="xy", type="l",ylim=c(min(RF$n[RF$n>0]),max(RF$n)))
lines(s_params@w,RF$n[2,], col="red")
sum(RF$n<0)
plot(project(used_params, t_max=15, effort = 0, initial_n = RF$n, t_save = 1, initial_n_pp = RF$n_pp))


#74 rewrote newton raphson solver in terms of `quasi-biomass`, and am running this for 
# the two species scale invariant model. Next is to see if this runs in BB

#74 I fixed the rescaled version, however for some reason the solver is now varyingnull points

#74 done a few random pertubations of two species case, then running mizer (150 yrs), then 
# newton raphson, so far my few observations support the hypothesis that 
# there is just one attractive interior steady state

#74 used 70 weight bins to speed up the code, and found another steady state
