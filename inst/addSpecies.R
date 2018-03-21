#! note we are assuming theta_ij=1, for all i,j

#! the code currently breaks under the default rfac=inf case, because it 
# relies on some values of rmax being passed through

params <- set_scaling_model(rfac = 10^10)
sim <- project(params, t_max=5, effort = 0)
plot(sim)

species_params <- params@species_params[10,]

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
combi_params@initial_n[new_sp,new_sp] <- 0
mumu <- getZ(combi_params, combi_params@initial_n, combi_params@initial_n_pp,0)[new_sp,]
gg <- getEGrowth(combi_params, combi_params@initial_n, combi_params@initial_n_pp)[new_sp,]

w_inf_idx <- length(combi_params@w[combi_params@w<=combi_params@species_params$w_inf[new_sp]])
integrand <- params@dw[combi_params@species_params$w_min_idx[new_sp]:w_inf_idx]*mumu[
  combi_params@species_params$w_min_idx[new_sp]:w_inf_idx]/gg[
    combi_params@species_params$w_min_idx[new_sp]:w_inf_idx]
combi_params@initial_n[new_sp,] <- 0
#! Here we are setting the abundance multiplier to 1. More useful to set it 
# to the value it should be if we were adding in another scale free species.
combi_params@initial_n[new_sp,combi_params@species_params$w_min_idx[new_sp]:w_inf_idx] <- 
  exp(-cumsum(integrand))/gg[combi_params@species_params$w_min_idx[new_sp]:w_inf_idx]
#! not really sure why the above seems to create a NAN in the last weight box.
combi_params@initial_n[is.nan(combi_params@initial_n)] <- 0

combi_params@initial_n[new_sp,new_sp] <- 1
################## reset erepro
#! this version readjusts erepro in the case where rmax=inf, need to think more 
# about what to do in general case, and look at wrapper functions.
gg0 <- gg[combi_params@species_params$w_min_idx[new_sp]]
mumu0 <- mumu[combi_params@species_params$w_min_idx[new_sp]]
rdi <- getRDI(combi_params, combi_params@initial_n, combi_params@initial_n_pp)[new_sp]
DW <- combi_params@dw[combi_params@species_params$w_min_idx[new_sp]]
combi_params@species_params$erepro[new_sp] <- combi_params@species_params$erepro[new_sp]*(
  combi_params@initial_n[new_sp,combi_params@species_params$w_min_idx[new_sp]]*(gg0+DW*mumu0))/rdi

##################
sim <- project(combi_params, t_max=5, effort = 0)
plot(sim)

# #20 Created a new branch called adsp, which is a copy of the 
# scaling branch, except that I have copied the params@p slots etc., 
# over from the density_dev branch. Started writing add_species. First 
# goal is to try and get this params combining code working. 
# Next step is to use getM2 etc., to generate the initial_n for 
# the new species.

# #20 Have got params code working without crashing. Next I need to: (1) 
# add code to pass forward the initial conditions for params, (2) write 
# getEgrowth and getZ based solver, to get the initial conditions for the 
# new species, (3) test the code by recombining species [1:(no_sp-1)] and 
# no_sp, from the scale free trait based model. (4) test code for adding 
# red mullet.

# #20 Got params combined code working basically.I think we should write 
# this code, so that if the newly added species is the next largest one 
# in the set_scaling model then it will work smoothly.
#The next steps are 
# (1) to put a new solver into the code to get the un-multiplied form 
# of the solution using old params data. To do this step, we can use 
# combi_params at the end of code, and mess with its theta matrix and solve. 
# Rather than not using an abundance multiplier for the new species, 
# we could it an abundance multiplier corresponding 
# to the abundance multiplier which would be associated with such a 
# species in the set_scaling model. (2) Test code when
# combining with scale free case. (3) Test code with mullet. (4) 
# make it as procedure, with help

# #20 and #42 . Added solver for new species, that determines its 
# abundance and erepro under the assumption of non-interaction. Marked specific issues with #! in code.
# Issues are (1)  maybe we do not have to change psi[new_sp,], (2) do we need to change
# params@srr when rmax off if it is off in two inputs ? The code currently breaks under the default rfac=inf case,
# because it relies on some values of rmax being passed through (3) Here we are setting the
# abundance multiplier to 1. More useful to set it to the value it should be if we 
# were adding in another scale free species. (4) Not really sure why the solution made 
# seems to have a NAN in the last weight box. (5) this version readjusts erepro in the case where rmax=inf, need to think more 
# about what to do in general case, and look at wrapper functions.





