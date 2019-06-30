plankton_logistic <- function(params, 
                              n = params@initial_n, 
                              n_pp = params@initial_n_pp, 
                              B = params@initial_B, 
                              rates = getRates(params), 
                              dt = 0.1, ...) {
    i <- 0.1 * params@w_full^(-params@lambda) * exp(-6.9*(lambda - 1))
    f <- params@rr_pp * n_pp * (1 - n_pp / params@cc_pp) + i - 
        rates$plankton_mort * n_pp 
    f[is.na(f)] <- 0
    return(n_pp + dt * f)
}

norm_box_pred_kernel <- function(ppmr, ppmr_min, ppmr_max) {
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

lambda <- 2
n <- 0.7

species_params <- data.frame(
    species = "Anchovy",
    w_min = 0.0003,
    w_mat = 10,
    m = 0.2 + n,
    w_inf = 66.5,
    erepro = 0.5,
    alpha = 0.1,
    ks = 0,
    gamma = 1920,# * 75000,
    ppmr_min = 100,
    ppmr_max = 30000,
    pred_kernel_type = "norm_box",
    h = Inf,
    r_max = Inf)


params <- set_multispecies_model(
    no_w = 200,
    species_params,
    lambda = lambda,
    kappa = 1 * exp(-6.9*(lambda - 1)),
    w_pp_cutoff = 0.1,
    q = 0.8,
    plankton_dynamics = plankton_logistic)
r0 <- 1
params@rr_pp[] <- r0 * params@w_full^(0.85 - 1)

# Larval mortality
mu_l <- 0
w_l <- 0.03
rho_l <- 5
mu_b <- mu_l / (1 + (params@w / w_l)^rho_l)
# Senescent mortality
mu_s <- 0.001
w_s <- 9
rho_s <- 5
mu_b[params@w >= w_s] <- (mu_s * (params@w / w_s)^rho_s)[params@w >= w_s]
params@mu_b[] <- mu_b

# no cannibalism
params@interaction[] <- 0

# initial anchovy abundance
params@initial_n[] <- 0.001 * params@w^(-1.8)

sim <- project(params, t_max = 150, dt = 0.05)
plotlySpectra(sim)
plotlyBiomass(sim)

sim <- project(sim, t_max = 200, dt = 0.01)
plotlyBiomass(sim)
plotlySpectra(sim, wlim = c(1e-10, NA), power = 1, ylim = c(1e-3, NA))

sim2 <- project(sim, dt = 0.01)
plotlyBiomass(sim2)
plotlySpectra(sim2, wlim = c(1e-10, NA), power = 1, ylim = c(1e-3, NA))

last <- dim(sim@n)[1]

plotGrowthCurves(sim, species = "Anchovy")
growth <- getEGrowth(params, sim@n[last, , ], sim@n_pp[last, ])
plot(params@w, growth, type = "l", log = "xy")

ws <- getGrowthCurves(sim, species = "Anchovy", max_age = 1000)
years <- colnames(ws)
years[1] <- 1
plot(years, ws[1,], type = "l", log = "x", 
     ylab = "weight [g]", xlab = "year")

pm <- getPlanktonMort(params, sim@n[last, , ], sim@n_pp[last, ])
plankton_biomass <- sim@n_pp[last, ] * params@w_full^2
plot(params@w_full, pm, type = "l", log = "xy",
     ylim = c(1e-1,100), col = "blue",
     main = "Cause of plankton collapse",
     xlab = "Size [g]",
     ylab = "Rate [1/year]")
lines(params@w_full, params@rr_pp, lty = "longdash")
lines(params@w_full, plankton_biomass / max(plankton_biomass) * 100,
      lty = "dotted", col = "green", lwd = 2)

plot(params@w_full, sim@n_pp[last, ] * params@w_full^2, type = "l", log = "xy")
