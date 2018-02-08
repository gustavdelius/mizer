library(progress)
library(mizer)

params_data_NS <- read.csv("./vignettes/NS_species_params.csv")
params <- MizerParams(params_data_NS)
no_sp <- dim(params@interaction)[1]
#mu0 <- 10
mu0 <- 100
#mu0 <- ((1-f0) * sqrt(2*pi) * kappa * gamma * sigma *
#          (beta^(n-1)) * exp(sigma^2 * (n-1)^2 / 2))[ii]
hbar <- params@species_params$alpha*params@species_params$h*params@f0-params@species_params$ks
# set these abundance multipliers later
sol_mult <- (1:no_sp)*10^(1)

pow <- (mu0/hbar/(1-params@n))

n_output <- params@psi
for (rep_idx in (1:no_sp)){
n_mult <- (1 - (params@w/params@species_params$w_inf[rep_idx])^(1-params@n))^(pow[rep_idx]-1) * (1 - (params@species_params$w_mat[rep_idx]/params@species_params$w_inf[rep_idx])^(1-params@n))^(-pow[rep_idx])
n_mult[params@w < params@species_params$w_mat[rep_idx]] <- 1
n_mult[params@w >= params@species_params$w_inf[rep_idx]] <- 0
n_exact <- params@w  # Just to get array with correct dimensions and names
n_exact <- ((params@species_params$w_min[rep_idx]/params@w)^(mu0/hbar[rep_idx]) / (hbar[rep_idx] * params@w^params@n)) * n_mult
n_exact[params@w < params@species_params$w_min[rep_idx]] <- 0
n_output[rep_idx,] <- sol_mult[rep_idx]*n_exact
}

################### setup plankton, background death and psi

plankton_vec <- (params@kappa*params@w^(-params@lambda))-colSums(n_output)
## better check sol_mult is small enough that plankton_vec has no negative entries
#plankton_vec[plankton_vec<0] <- 0
#plankton_vec[min(which(plankton_vec==0)):length(plankton_vec)] <- 0
params@cc_pp[sum(params@w_full<=w[1]):length(params@cc_pp)] <- plankton_vec
initial_n_pp <- params@cc_pp
# The cc_pp factor needs to be higher than the desired steady state in
# order to compensate for predation mortality
m2_background <- getM2Background(params, n_output, initial_n_pp)
params@cc_pp <- (1+m2_background/params@rr_pp) * initial_n_pp

# ----
#' ### Setup background death and steplike psi

m2 <- getM2(params, n_output, initial_n_pp)

for (i in 1:no_sp) {
  params@psi[i, ] <- (params@w/params@species_params$w_inf[i])^(1-params@n)
  params@psi[i, params@w < (params@species_params$w_mat[i]-1e-10)] <- 0
  params@psi[i, params@w > (params@species_params$w_inf[i]-1e-10)] <- 1
  params@mu_b[i, ] <- mu0 * params@w^(params@n-1) - m2[i,]
}

## better check mu0 is large enough that mu_b has no negative entries


#' ### Set erepro to meet boundary condition

rdi <- getRDI(params, n_output, initial_n_pp)
gg <- getEGrowth(params, n_output, initial_n_pp)
mumu <- getZ(params, n_output, initial_n_pp, effort = 0)
erepro_final <- rdi
for (i in (1:no_sp)){
  #  erepro_final[i] <- erepro*(gg[i,params@species_params$w_min_idx[i]]*n_output[i,params@species_params$w_min_idx[i]])/
  #    rdi[i]
  gg0 <- gg[i,params@species_params$w_min_idx[i]]
  mumu0 <- mumu[i,params@species_params$w_min_idx[i]]
  DW <- params@dw[params@species_params$w_min_idx[i]]
  erepro_final[i] <- params@species_params$erepro[i]*(n_output[i,params@species_params$w_min_idx[i]]*(gg0+DW*mumu0))/rdi[i]
}



# turn on density dependence 

# add density dependence
nn <- n_output
nn[nn==0] <- 1
params@chi <- 0.5
params@ddd <- nn^(params@chi)

# run simulation
n_output2 <- n_output
n_output2[no_sp,] <- 10*n_output[no_sp,]
t_max <- 400
sim <- project(params, t_max=t_max ,dt=0.1, t_save=t_max/100, effort = 0, 
               initial_n = n_output2, initial_n_pp = initial_n_pp)
plot(sim)

plot(params@w,plankton_vec,log="xy")
plot(params@w,n_output[1,],log="xy")
plot(params@w,n_exact,log="xy")

sum(params@mu_b<0)
sum(plankton_vec<0)

plot(params@w,sim@n[dim(sim@n)[1],1,],log="xy")
for (i in 2:no_sp){
  lines(params@w,sim@n[dim(sim@n)[1],i,])
}
# plotSpectra(sim) fails since min(sim@n[sim@n>0]) is so small
