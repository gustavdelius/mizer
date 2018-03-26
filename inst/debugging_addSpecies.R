
# #20 #42 run addSpecies_combine.R first, then run this for debugging after hake is added. 
# the code below line 5 is lifed from add_species, to spot bug after we add have

A
params <- combi_params


combi_species_params <- rbind(params@species_params, species_params)
if (length(combi_species_params$r_max) == 0) {
  combi_species_params$r_max <- 10 ^ (50)
}
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
    r_pp = (params@rr_pp / (params@w_full ^ (params@p - 1)))[1]
  )


combi_params@initial_n_pp <- params@initial_n_pp
combi_params@cc_pp <- params@cc_pp

new_sp <- length(params@species_params$species) + 1

combi_params@initial_n[1:(new_sp - 1), ] <- params@initial_n
combi_params@species_params$erepro[1:(new_sp - 1)] <-
  params@species_params$erepro
combi_params@psi[1:(new_sp - 1), ] <- params@psi
combi_params@mu_b[1:(new_sp - 1), ] <- params@mu_b

# other important info to pass through correspond to
# the parts of params that got modified after set_scaling made it initially.
combi_params@species_params$erepro[new_sp] <- 0.1
#! maybe we do not have to change psi[new_sp,]
combi_params@psi[new_sp, ] <-
  (combi_params@w / combi_params@species_params$w_inf[new_sp]) ^ (1 - combi_params@n)
combi_params@psi[new_sp, combi_params@w < (combi_params@species_params$w_mat[new_sp] - 1e-10)] <-
  0
combi_params@psi[new_sp, combi_params@w > (combi_params@species_params$w_inf[new_sp] - 1e-10)] <-
  1
combi_params@mu_b[new_sp, ] <- params@mu_b[(new_sp - 1), ]
#! what about params@srr ? do I have to pass this through when rmax is off ?
#! do I have to set rmax off if it is off in two inputs ?
combi_params@srr <- params@srr
# use rest of info to fill out n_new_sp properly
# combi_params@initial_n[new_sp,] <- params@initial_n[(new_sp-1),]
################################
combi_params@interaction[new_sp, new_sp] <- 0
sum(is.na(combi_params@initial_n))
mumu <-
  getZ(combi_params,
       combi_params@initial_n,
       combi_params@initial_n_pp,
       effort = 0)[new_sp, ]
sum(is.na(mumu))