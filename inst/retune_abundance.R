


#@ retune_abundance <- function(params,A){
params <- set_scaling_model()
A <- rep(1,length(params@species_params$w_inf))
A[4:7] <- NA



# get a list of the species that we can tune abundance mulpliers of
all_background <- (1:length(A))[is.na(A)]
largest_background <- which.max(params@species_params$w_inf[all_background])
A2 <- A
# we are assuming that the the abundance multiplier of the largest backgroud 
# species should be held fixed at 1, if though it was initially an NA
A2[largest_background] <- 1


# we make a list L of species we will vary the abundance parameters of
# (everything but largest background)
#L <- all_background[all_background!=largest_background]
L <- (1:length(A))[is.na(A2)]
idx_start <- sum(params@w<=min(params@species_params$w_mat))
idx_stop <- sum(params@w<=max(params@species_params$w_inf))
RR <- matrix(0,nrow = length(L), ncol = length(L))
QQ <- (1:length(L))

Lcomp <- (1:length(A))[!is.na(A2)]

old_n <- params@initial_n
no_sp <- length(params@species_params$w_inf)

#A3 <- A2
#A3[is.na(A3)] <- 1
#for (i in 1:no_sp){
#  old_n[i,] <- A3[i]*params@initial_n[i,]
#}
#cc <- params@kappa*params@w^(-params@lambda)
cc <- colSums(params@initial_n[all_background,])
#rho <- colSums(params@initial_n[Lcomp,])
rho <- colSums(params@initial_n[Lcomp,])

den <- cc^2
den[den==0] <- 10^(-50)
for (i in (1:length(L))){
  QQ[i] <- sum((params@initial_n[L[i],]*(cc-rho)*params@dw/(den))[idx_start:(idx_stop-1)])
  for (j in (1:length(L))){
    RR[i,j] <- sum((params@initial_n[L[i],]*params@initial_n[L[j],]*params@dw/(den))[idx_start:(idx_stop-1)])
  }
}

A2[L] <- solve(RR,QQ)

#@  return(A2)}
new_n <- params@initial_n
for (i in 1:no_sp){
  new_n[i,] <- A2[i]*params@initial_n[i,]
}
A2
plot(params@w,new_n[1,],log="xy", ylim=c(10^(-3),max(new_n)))
for (i in (2:no_sp)){
  lines(params@w,new_n[i,])
}

plot(params@w[idx_start:(idx_stop-1)],(colSums(new_n)/(params@kappa*params@w^(-params@lambda)))[idx_start:(idx_stop-1)],log="x")

plot(params@w,colSums(new_n),log="xy")
A2

# #20 #42  Used colSums(params@initial_n[all_background,]) for cc now by only summing over the background species, and 
# adjusted denominator cc^2 in RR and QQ, so division by 0 does not occur. 
# Why does this code give -ve abundance multipliers. 
# but now the code gives -ve abundances, and does not return pow law.