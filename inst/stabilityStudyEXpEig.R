library("BB")
params <- set_scaling_model(no_sp = 5)


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
max(Re(EE$values))

lead <- which(Re(EE$values)==max(Re(EE$values)))
Re(EE$values[lead])
Im(EE$values[lead])

min(Mod(EE$values))

# check these eigen things each satisfy the eigen equation

eig_error <- EE$values
#for (i in (1:length(EE$values))){
  for (i in (1:100)){
    
  ans2 <- (jacobian_numeric %*% EE$vectors[,i] - EE$values[i]*EE$vectors[,i])
eig_error[i] <- max(Mod(ans2))
}
max(Mod(eig_error[1:100]))

i <- lead[1]
ans_lead <- (jacobian_numeric %*% EE$vectors[,i] - EE$values[i]*EE$vectors[,i])
max(Mod(ans_lead))


############

tangent_vector <- EE$vectors[,i]
fixed_pt <- c(params@initial_n_pp, params@initial_n)
new_start <- fixed_pt+(10^(-5))*Re(tangent_vector)

#ss_go <- multiroot(f = change_function, start = new_start)

# too slow to check this, and we should use gustav's faster algorithm to see if there are 
# many fixed points 

tang_n_pp <- Re(tangent_vector[1:length(params@cc_pp)])
#n <- data_in[(length(used_params@cc_pp)+1):length(data_in)]
# matrix(vec,nrow = 7,ncol = 7)
tang_n <- Re(array(tangent_vector[(length(params@cc_pp)+1):length(tangent_vector)], dim = dim(params@initial_n)))

plot(params@w_full,tang_n_pp,log="x")

plot(params@w,tang_n[1,],log="x",ylim=c(min(tang_n),max(tang_n)),type="l")
for (j in 1:dim(tang_n)[1]){
  lines(params@w,tang_n[j,])
}


####################################

plot(params@w_full,tang_n_pp,log="xy")


plot(params@w_full,tang_n_pp*params@w_full^(params@lambda),log="x")

HH <- sweep(tang_n, 2, params@w^(params@lambda), "*")

plot(params@w,tang_n[1,]*params@w^(params@lambda),log="x",ylim=c(min(HH),max(HH)),type="l")
for (j in 2:(dim(tang_n)[1])){
  lines(params@w,tang_n[j,]*params@w^(params@lambda),col=gray(1-j/5))
}

plot(params@w,tang_n[1,]*params@w^(params@lambda),log="x",ylim=c(min(HH),max(HH)),type="l")
for (j in 2:(dim(tang_n)[1])){
  lines(params@w,tang_n[j,]*params@w^(params@lambda),col=gray((j-1)/5))
}


for (j in 1:5){
  print(params@w[max((1:length(tang_n[j,]))[tang_n[j,]!=0])])
}

for (j in 1:5){
  print(c(params@species_params$w_min[j],params@w[min((1:length(tang_n[j,]))[tang_n[j,]!=0])],tang_n[j,params@species_params$w_min_idx[j]]*params@species_params$w_min[j]^params@lambda
))
}

ss <- tang_n[1,]
ss[] <- 0
for (j in 1:5){
  ss <- ss + tang_n[j,]
}
plot(params@w,ss,log="x")


plot(sapply(EE$values, function(x) Mod(x))) 

EE$values[lead[1]]-Conj(EE$values[lead[2]])

sum(EE$vectors[lead[1]]-Conj(EE$vectors[lead[2]]))


tangent_vectorA <- EE$vectors[,lead[1]]

tang_n_ppI <- Im(tangent_vectorA[1:length(params@cc_pp)])
#n <- data_in[(length(used_params@cc_pp)+1):length(data_in)]
# matrix(vec,nrow = 7,ncol = 7)
tang_nI <- Im(array(tangent_vectorA[(length(params@cc_pp)+1):length(tangent_vectorA)], dim = dim(params@initial_n)))



plot(params@w_full,tang_n_ppI*params@w_full^(params@lambda),log="x")

HH <- sweep(tang_nI, 2, params@w^(params@lambda), "*")

plot(params@w,tang_nI[1,]*params@w^(params@lambda),log="x",ylim=c(min(HH),max(HH)),type="l")
for (j in 2:(dim(tang_nI)[1])){
  lines(params@w,tang_nI[j,]*params@w^(params@lambda),col=gray((j-1)/5))
}

t <- 80
a <- Re(EE$values[lead[1]])
b <- Im(EE$values[lead[1]])
alpha <- Re(EE$vectors[,lead[1]])
beta <- Im(EE$vectors[,lead[1]])
osc <- exp(a*t)*(alpha*cos(b*t) - beta*sin(b*t))

################

tang_n_ppS <- osc[1:length(params@cc_pp)]
#n <- data_in[(length(used_params@cc_pp)+1):length(data_in)]
# matrix(vec,nrow = 7,ncol = 7)
tang_nS <- array(osc[(length(params@cc_pp)+1):length(osc)], dim = dim(params@initial_n))

#plot(params@w_full,tang_n_ppS*params@w_full^(params@lambda),log="x")

HH <- sweep(tang_nS, 2, params@w^(params@lambda), "*")

plot(params@w,tang_nS[1,]*params@w^(params@lambda),log="x",ylim=c(min(HH),max(HH)),type="l")
for (j in 2:(dim(tang_nS)[1])){
  lines(params@w,tang_nS[j,]*params@w^(params@lambda),col=gray((j-1)/5))
}


################



# plot leading eigenvector


# dot product of normalized eigenvectors (any close to parrallel ?)



conditional <- kappa(jacobian_numeric)
conditional

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


####################

library("RSpectra")

EEnew <- eigs(jacobian_numeric, 1, which = "LR")$value
EEnew
EEnewB <- eigs(jacobian_numericB, 1, which = "LR")$value
EEnewB

eigs(matrix(1:9,ncol=3,nrow=3), 1, which = "LR")$value

eigen(matrix(1:9,ncol=3,nrow=3))$values


