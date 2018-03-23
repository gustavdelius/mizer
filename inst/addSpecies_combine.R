
############# retune_abundance code for later use




retune_abundance <- function(params,A){
#params <- set_scaling_model()
#A <- rep(1,length(params@species_params$w_inf))
#A[4:7] <- NA



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

A3 <- A2
A3[is.na(A3)] <- 1

for (i in 1:no_sp){
  old_n[i,] <- A3[i]*params@initial_n[i,]
}
#cc <- params@kappa*params@w^(-params@lambda)
cc <- colSums(old_n)
#rho <- colSums(params@initial_n[Lcomp,])
rho <- colSums(old_n[Lcomp,])


for (i in (1:length(L))){
  QQ[i] <- sum((params@initial_n[L[i],]*(cc-rho)*params@dw/(cc^2))[idx_start:(idx_stop-1)])
  for (j in (1:length(L))){
    RR[i,j] <- sum((params@initial_n[L[i],]*params@initial_n[L[j],]*params@dw/(cc^2))[idx_start:(idx_stop-1)])
  }
}

A2[L] <- solve(RR,QQ)

  return(A2)}
params <- set_scaling_model()
A <- rep(1,length(params@species_params$w_inf))
A[4:7] <- NA
A2 <- retune_abundance(params,A)
A2
#new_n <- params@initial_n
#for (i in 1:no_sp){
#  new_n[i,] <- A2[i]*params@initial_n[i,]
#}
#A2
#plot(params@w,new_n[1,],log="xy", ylim=c(10^(-3),max(new_n)))
#for (i in (2:no_sp)){
#  lines(params@w,new_n[i,])
#}
#plot(params@w[idx_start:(idx_stop-1)],(colSums(new_n)/(params@kappa*params@w^(-params@lambda)))[idx_start:(idx_stop-1)],log="x")



############################
########## add_species code below, ready for wrap, inputs are params & species_params


#! note we are assuming theta_ij=1, for all i,j

#! the code currently breaks under the default rfac=inf case, because it 
# relies on some values of rmax being passed through

#params <- set_scaling_model(rfac = 10^10)

params <- set_scaling_model()

if (length(params@species_params$r_max)==0){
  params@species_params$r_max <- rep(10^(50),length(params@species_params$w_inf))
}

sim <- project(params, t_max=5, effort = 0)
plot(sim)

species_params <- params@species_params[10,]
species_params$beta <- 50

if (length(species_params$r_max)==0){
  species_params$r_max <- 10^(50)
}


#####################

combi_species_params <- data.frame(
  species = 1:(length(params@species_params$species)+length(species_params$species)),
  w_min = c(params@species_params$w_min,species_params$w_min),
  w_inf = c(params@species_params$w_inf,species_params$w_inf),
  w_mat = c(params@species_params$w_mat,species_params$w_mat),
  h = c(params@species_params$h,species_params$h),
  ks = c(params@species_params$ks,species_params$ks),
  beta = c(params@species_params$beta,species_params$beta),
  sigma = c(params@species_params$sigma,species_params$sigma),
  z0 = 0,
  alpha = c(params@species_params$alpha,species_params$alpha),
  erepro = c(params@species_params$erepro,species_params$erepro),
  sel_func = "knife_edge",
  # not used but required
  knife_edge_size = c(params@species_params$knife_edge_size,species_params$knife_edge_size),
  gear = c(params@species_params$gear,species_params$gear),
  r_max = c(params@species_params$r_max,species_params$r_max)
)

combi_params <-
  MizerParams(
    combi_species_params,
    p = params@p,
    n = params@n,
    q = params@q,
    lambda = params@lambda,
    f0 = params@f0,
    kappa = params@kappa,
    min_w = min(params@w),
    max_w = max(params@w),
    no_w = length(params@w),
    min_w_pp = min(params@w_full),
    w_pp_cutoff = max(params@w_full),
    r_pp = (params@rr_pp/(params@w_full^(params@p-1)))[1]
  )


combi_params@initial_n_pp <- params@initial_n_pp
combi_params@cc_pp <- params@cc_pp

new_sp <- length(params@species_params$species)+1
for (i in (1:(new_sp-1))){
  combi_params@initial_n[i,] <- params@initial_n[i,]
  combi_params@species_params$erepro[i] <- params@species_params$erepro[i] 
  combi_params@psi[i,] <- params@psi[i,]
  combi_params@mu_b[i,] <- params@mu_b[i,]
}
# other important info to pass through correspond to 
# the parts of params that got modified after set_scaling made it initially.
combi_params@species_params$erepro[new_sp] <- params@species_params$erepro[(new_sp-1)]
#! maybe we do not have to change psi[new_sp,]
combi_params@psi[new_sp,] <- (combi_params@w / combi_params@species_params$w_inf[new_sp]) ^ (1 - combi_params@n)
combi_params@psi[new_sp, combi_params@w < (combi_params@species_params$w_mat[new_sp] - 1e-10)] <- 0
combi_params@psi[new_sp, combi_params@w > (combi_params@species_params$w_inf[new_sp] - 1e-10)] <- 1
combi_params@mu_b[new_sp,] <- params@mu_b[(new_sp-1),]
#! what about params@srr ? do I have to pass this through when rmax is off ?
#! do I have to set rmax off if it is off in two inputs ?
# use rest of info to fill out n_new_sp properly
# combi_params@initial_n[new_sp,] <- params@initial_n[(new_sp-1),]
################################
combi_params@interaction[new_sp,new_sp] <- 0
mumu <- getZ(combi_params, combi_params@initial_n, combi_params@initial_n_pp,0)[new_sp,]
gg <- getEGrowth(combi_params, combi_params@initial_n, combi_params@initial_n_pp)[new_sp,]

w_inf_idx <- length(combi_params@w[combi_params@w<=combi_params@species_params$w_inf[new_sp]])
integrand <- params@dw[combi_params@species_params$w_min_idx[new_sp]:w_inf_idx]*mumu[
  combi_params@species_params$w_min_idx[new_sp]:w_inf_idx]/gg[
    combi_params@species_params$w_min_idx[new_sp]:w_inf_idx]
combi_params@initial_n[new_sp,] <- 0
#! Here we are setting the abundance multiplier to 1. More useful to set it 
# to the value it should be if we were adding in another scale free species.
mult <- 1.5*10^(11)
combi_params@initial_n[new_sp,combi_params@species_params$w_min_idx[new_sp]:w_inf_idx] <- 
  mult*exp(-cumsum(integrand))/gg[combi_params@species_params$w_min_idx[new_sp]:w_inf_idx]
#! not really sure why the above seems to create a NAN in the last weight box.
combi_params@initial_n[is.nan(combi_params@initial_n)] <- 0

combi_params@interaction[new_sp,new_sp] <- 1
################## reset erepro
#! this version readjusts erepro in the case where rmax=inf, need to think more 
# about what to do in general case, and look at wrapper functions.
gg0 <- gg[combi_params@species_params$w_min_idx[new_sp]]
mumu0 <- mumu[combi_params@species_params$w_min_idx[new_sp]]
rdi <- getRDI(combi_params, combi_params@initial_n, combi_params@initial_n_pp)[new_sp]
DW <- combi_params@dw[combi_params@species_params$w_min_idx[new_sp]]
#combi_params@species_params$erepro[new_sp] <- combi_params@species_params$erepro[new_sp]*(
#  combi_params@initial_n[new_sp,combi_params@species_params$w_min_idx[new_sp]]*(gg0+DW*mumu0))/rdi

H <- rdi/combi_params@species_params$erepro[new_sp]
X <- combi_params@initial_n[new_sp,combi_params@species_params$w_min_idx[new_sp]]*(gg0+DW*mumu0)

combi_params@species_params$erepro[new_sp] <- 
  combi_params@species_params$r_max[new_sp]*X/(H*combi_params@species_params$r_max[new_sp]+H*X)
##################
sim <- project(combi_params, t_max=5, effort = 0)
plot(sim)

###########
A <- rep(NA,length(combi_params@species_params$w_inf))
A[length(combi_params@species_params$w_inf)] <- 1
AA <- retune_abundance(combi_params,A)
new_n <- combi_params@initial_n
for (i in 1:length(combi_params@species_params$w_inf)){
new_n[i,] <- AA[i]*combi_params@initial_n[i,]
}

AA
plot(params@w,new_n[1,],log="xy", ylim=c(10^(-3),max(new_n)))
for (i in (2:no_sp)){
  lines(params@w,new_n[i,])
}
plot(params@w,colSums(new_n),log="xy")

# #20 #42 Started addSpecies_combine.R where I plugged together add_species 
# and retune abundance, although 
# retune abundance always seems to return 1's
