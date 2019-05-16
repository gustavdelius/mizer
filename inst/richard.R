plankton_logistic <- function(params, n, n_pp, B, rates, dt, ...) {
    i <- 0.1 * params@w_full^(-params@lambda)
    f <- params@rr_pp * n_pp * (1 - n_pp / params@cc_pp) + i - 
        rates$plankton_mort 
    return(n_pp + dt * f)
}

norm_box_pred_kernel <- function(ppmr, ppmr_min, ppmr_max) {
    assert_that(ppmr_min < ppmr_max)
    phi <- rep(1, length(ppmr))
    phi[ppmr > ppmr_max] <- 0
    phi[ppmr < ppmr_min] <- 0
    # Do not allow feeding at own size
    phi[1] <- 0
    # normalise in log space
    logppmr <- log(ppmr)
    dl <- logppmr[2] - logppmr[1]
    N <- sum(phi) * dl
    phi <- phi / N
    return(phi)
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
    pred_kernel_type = "norm_box",
    h = Inf,
    r_max = Inf)

lambda <- 2

params <- set_multispecies_model(
    species_params,
    lambda = lambda,
    kappa = 100 * 10^(-3*lambda),
    w_pp_cutoff = 0.1,
    q = 0.8)

plotlySpectra(params)

sim <- project(params)
plotlySpectra(sim)
plotlyBiomass(sim)

sim <- project(sim)
plotlyBiomass(sim)
