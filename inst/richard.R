plankton_logistic <- function(params, n, n_pp, B, rates, dt, ...) {
    i <- 0.1 * params@w_full^(-params@lambda)
    f <- params@rr_pp * n_pp * (1 - n_pp / params@cc_pp) + i - 
        rates$plankton_mort 
    return(n_pp + dt * f)
}

species_params <- data.frame(
    species = "Anchovy",
    w_min = 0.0003,
    w_mat = 10,
    w_inf = 10,
    erepro = 0.5,
    alpha = 0.1,
    ks = 0,
    gamma = 2000,
    ppmr_min = 100,
    ppmr_max = 30000,
    pred_kernel_type = "box",
    h = Inf)

lambda <- 2

params <- set_multispecies_model(
    species_params,
    lambda = lambda,
    kappa = 100 * 10^(-3*lambda),
    w_pp_cutoff = 0.1,
    q = 0.8)

plotlySpectra(params)

sim <- project(params)
