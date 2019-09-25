# Set up model ####

setAnchovyMort <- 
  function(params,
           # background mortality
           mu_0 = 1,
           w_0 = params@species_params$w_min,
           rho_b = -0.25,
           # Larval mortality
           mu_l = 0,
           w_l = 0.03,
           rho_l = 5,
           # Senescent mortality
           mu_s = 0.001,
           w_s = 0.5,
           rho_s = 1) {
    mu_b <- rep(0, length(params@w))
    mu_b[params@w <= w_s] <- (mu_0 * (params@w / w_0)^rho_b)[params@w < w_s]
    if (mu_0 > 0) {
      mu_s <- min(mu_b[params@w <= w_s])
    }
    mu_b[params@w >= w_s] <- (mu_s * (params@w / w_s)^rho_s)[params@w >= w_s]
    # Add larval mortality
    mu_b <- mu_b + mu_l / (1 + (params@w / w_l)^rho_l)
    
    params@mu_b[] <- mu_b
    return(params)
  }

plankton_logistic <- function(params, 
                              n = params@initial_n, 
                              n_pp = params@initial_n_pp, 
                              B = params@initial_B, 
                              rates = getRates(params), 
                              dt = 0.1, ...) {
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

setModel <- function(
  lambda = 2,
  q = 0.85,
  a0 = 100,
  r0 = 10,
  gamma = 1920 #1920,# * 75000
  ) {
  kappa = a0 * exp(-6.9*(lambda - 1))
  n = 2/3 # irrelevant value
  
  species_params <- data.frame(
    species = "Anchovy",
    w_min = 0.0003,
    w_mat = 10,
    m = 0.2 + n,
    w_inf = 66.5,
    erepro = 0.5,
    alpha = 0.1,
    ks = 0,
    gamma = gamma,
    ppmr_min = 100,
    ppmr_max = 30000,
    pred_kernel_type = "norm_box",
    h = Inf,
    r_max = Inf)
  
  
  params <- set_multispecies_model(
    no_w = 200,
    species_params,
    lambda = lambda,
    kappa = kappa,
    w_pp_cutoff = 0.1,
    q = q,
    plankton_dynamics = plankton_logistic)
  
  params@rr_pp[] <- r0 * params@w_full^(0.85 - 1)
  return(params)
}

# Old case ####
i0 <- 0.1
params <- setModel(a0 = 1, r0 = 1, q = 0.8) %>% 
  setAnchovyMort(mu_l = 0, mu_0 = 0)
# no cannibalism
params@interaction[] <- 0
# initial anchovy abundance
params@initial_n[] <- 0.001 * params@w^(-1.8)
params@initial_n_pp[] <- params@cc_pp

sim <- project(params, t_max = 40, dt = 0.01)
plotlyBiomass(sim)
plotlySpectra(sim, wlim = c(1e-10, NA), power = 2, ylim = c(1e-5, NA))

# New case ####
i0 <- 10
params <- setModel(a0 = 100, r0 = 10, q = 0.85, gamma = 750) %>% 
  setAnchovyMort(mu_l = 0)
i <- i0 * params@w_full^(-params@lambda) * exp(-6.9*(params@lambda - 1))
# no cannibalism
params@interaction[] <- 0
# initial anchovy abundance
params@initial_n[] <- 0.001 * params@w^(-1.8)
params@initial_n_pp[] <- params@cc_pp

sim <- project(params, t_max = 20, dt = 0.001)
plotlyBiomass(sim)
plotlySpectra(sim, wlim = c(1e-10, NA), power = 2, ylim = c(1e-5, NA))

sim2 <- project(sim, t_max = 10, dt = 0.001)
plotlyBiomass(sim2)
plotlySpectra(sim2, wlim = c(1e-10, NA), power = 2, ylim = c(1e-5, NA))
sim2 <- project(sim2, t_max = 10, dt = 0.001)

last <- dim(sim2@n)[1]
params@initial_n[] <- sim2@n[last, , ]
params@initial_n_pp[] = sim2@n_pp[last, ]

# Growth curves ####
growth <- getEGrowth(params)
plot(params@w, growth, type = "l", log = "xy")

ws <- getGrowthCurves(sim2, species = "Anchovy", max_age = 1000)
years <- colnames(ws)
years[1] <- 1
plot(years, ws[1,], type = "l", log = "x", 
     ylab = "weight [g]", xlab = "year")

# Plankton mortality ####
pm <- getPlanktonMort(params)
plankton_biomass <- params@initial_n_pp * params@w_full^2
plot(params@w_full, pm, type = "l", log = "xy",
     ylim = c(1e-1,100), col = "blue",
     main = "Cause of plankton collapse",
     xlab = "Size [g]",
     ylab = "Rate [1/year]")
lines(params@w_full, params@rr_pp, lty = "longdash")
lines(params@w_full, plankton_biomass / max(plankton_biomass) * 100,
      lty = "dotted", col = "green", lwd = 2)

plot(params@w_full, params@initial_n_pp * params@w_full^2, type = "l", log = "xy")


# Switch on cannibalism ####
params@interaction[] <- 1
# No larval mortality
params <- setAnchovyMort(params, mu_l = 0)
# initial anchovy abundance
params@initial_n[] <- 0.001 * params@w^(-1.8)
params@initial_n_pp[] <- params@cc_pp

# Simulate ####
sim <- project(params, t_max = 20, dt = 0.001)
plotlyBiomass(sim)
plotlySpectra(sim, wlim = c(1e-10, NA), power = 2, ylim = c(1e-5, NA))

last = dim(sim@n)[1]
mu <- getMort(params, n = sim@n[last, , ], n_pp = sim@n_pp[last, ], effort = 0)[1, ]
plot(params@w, mu, type = "l", log = "x")

sim2 <- project(sim, t_max = 10, dt = 0.001)
plotlyBiomass(sim2)
plotlySpectra(sim2, wlim = c(1e-10, NA), power = 2, ylim = c(1e-5, NA),
              time_range = 2)
sim2 <- project(sim2, t_max = 10, dt = 0.001)

sim <- project(params, dt = 0.01)
plotBiomass(sim)
plotSpectra(sim)
sim <- project(sim, dt = 0.01)

# Switch on larval mortality ####
params <- setAnchovyMort(params, mu_l = 21)
# no cannibalism
params@interaction[] <- 0
# initial anchovy abundance
params@initial_n[] <- 0.001 * params@w^(-1.8)
params@initial_n_pp[] <- params@cc_pp

# Simulate
sim <- project(params, dt = 0.001)
plotBiomass(sim)
plotSpectra(sim)
sim <- project(sim, dt = 0.01)

params <- setAnchovyMort(params, mu_l = 0.01)
sim <- project(params, dt = 0.01)
plotBiomass(sim)
plotSpectra(sim)

params@initial_n[] <- 0.001 * params@w^(-1.8)
params <- setAnchovyMort(params, mu_l = 0.2)
sim <- project(params, dt = 0.01)
plotBiomass(sim)
plotSpectra(sim)
sim <- project(sim, dt = 0.01)
plotBiomass(sim)
plotSpectra(sim)


params <- setAnchovyMort(params, mu_l = 1)
sim <- project(params, dt = 0.01)
plotBiomass(sim)


