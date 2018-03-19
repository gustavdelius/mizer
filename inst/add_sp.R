params <- set_scaling_model()
sim <- project(params, t_max=5, effort = 0)
plot(sim)

species_params <- params@species_params

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
  gear = c(params@species_params$gear,species_params$gear)
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
    r_pp = params@rr_pp
  )
