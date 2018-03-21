params <- set_scaling_model()
sim <- project(params, t_max=5, effort = 0)
plot(sim)

species_params <- params@species_params[11,]

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
combi_params@psi[new_sp,] <- (combi_params@w / combi_params@species_params$w_inf[new_sp]) ^ (1 - combi_params@n)
combi_params@psi[new_sp, combi_params@w < (combi_params@species_params$w_mat[new_sp] - 1e-10)] <- 0
combi_params@psi[new_sp, combi_params@w > (combi_params@species_params$w_inf[new_sp] - 1e-10)] <- 1
combi_params@mu_b[new_sp,] <- params@mu_b[(new_sp-1),]
# what about params@srr ? do I have to pass this through when rmax is off ?
# use rest of info to fill out n_new_sp properly
combi_params@initial_n[new_sp,] <- params@initial_n[(new_sp-1),]



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
