# Our first attempt at setting up a mizer model for the Blanes system
species_params <- read.csv("inst/blanes/species_params.csv")

species_params$r_max <- Inf
species_params$ks <- 0

theta <- read.csv("inst/blanes/theta.csv", header = FALSE)
theta <- t(as.matrix(theta))

f0 <- 0.6
q <- 0.88
n <- 0.7
p <- 0.7
lambda <- 2 + q - n
kappa <- 1

no_sp <- length(species_params$species)

params <- MizerParams(species_params,
                      f0 = f0, kappa = kappa, p = p, n = n, q = q, 
                      lambda = lambda, z0pre = 2)

# Set up resources
resource_params <- list("detritus_external" = 0,
                        "carrion_external" = 0,
                        "background_external" = 0)
resource_dynamics <- 
    list("detritus" = function(params, n, n_pp, B, rates, dt, ...) B["detritus"],
         "carrion" = function(params, n, n_pp, B, rates, dt, ...) B["carrion"],
         "background" = function(params, n, n_pp, B, rates, dt, ...) B["background"])
params <- setResourceDynamics(params,
                              resource_dynamics = resource_dynamics,
                              resource_params = resource_params)
params@initial_B <- c(detritus = 1, carrion = 1, background = 1)

# Start with no interaction
interaction0 <- theta
interaction0[] <- 0
interaction_p <- params@species_params$interaction_p
params@species_params$interaction_p[] <- 0
params <- setInteraction(params, interaction = interaction0)

# Make everything feed on detritus
rate <- f0 / (1 - f0) * params@species_params$h
params@species_params$rho_detritus <- rate * species_params$detQ
params@species_params$rho_carrion <- rate * species_params$scavQ
params@species_params$rho_background <- 
    rate * !(species_params$detQ | species_params$scavQ)
params <- setResourceEncounter(params)

# Set initial solution in fixed background
mumu <- getMort(params, effort = 0)  # Death rate
gg <- getEGrowth(params)  # Growth rate
for (i in 1:no_sp) {
    mu <- mumu[i, ]
    g <- gg[i, ]
    w_inf_idx <- sum(params@w <= species_params$w_inf[i])
    idx <- params@w_min_idx[i]:(w_inf_idx - 1)
    idxs <- params@w_min_idx[i]:(w_inf_idx)
    # Steady state solution of the upwind-difference scheme used in project
    params@initial_n[i, idxs] <- 
        c(1, cumprod(g[idx] / ((g + mu * params@dw)[idx + 1])))
    # rescale abundance to agree with observation
    if (is.na(species_params$cutoff_size[i])) {
        stop(paste("No cutoff size is given for ",
                   species_params$species[i]))
    }
    sel <- params@w >= species_params$cutoff_size[i]
    if (!is.na(species_params$abundance_observed[i])) {
        abundance <- sum(params@initial_n[i, sel] * params@dw[sel])
        params@initial_n[i, ] <- params@initial_n[i, ] / abundance *
            species_params$abundance_observed[i]
    }
    if (!is.na(species_params$biomass_observed[i])) {
        biomass <- sum(params@initial_n[i, sel] * params@w[sel] * params@dw[sel])
        params@initial_n[i, ] <- params@initial_n[i, ] / biomass *
            species_params$biomass_observed[i]
    }
    
}

# Retune the values of erepro, so that we are at steady state.
rdi <- getRDI(params)
for (i in 1:no_sp) {
    gg0 <- gg[i, params@w_min_idx[i]]
    mumu0 <- mumu[i, params@w_min_idx[i]]
    DW <- params@dw[params@w_min_idx[i]]
    if (rdi[i] == 0) {
        stop("No reproduction")
    }
    params@species_params$erepro[i] <- params@species_params$erepro[i] *
            (params@initial_n[i, params@w_min_idx[i]] *
                 (gg0 + DW * mumu0)) / rdi[i]
}

params <- steady(params, t_max = 100)
plotlySpectra(params)

# Slowly switch on interaction
epsilon <- 1
params@species_params$interaction_p[] <- interaction_p * epsilon
params <- setInteraction(params, interaction = theta * epsilon)

params <- steady(params, t_max = 100)
plotlySpectra(params)

# Rescale abundance again
mumu <- getMort(params, effort = 0)  # Death rate
gg <- getEGrowth(params)  # Growth rate
for (i in 1:no_sp) {
    if (is.na(species_params$cutoff_size[i])) {
        stop(paste("No cutoff size is given for ",
                   species_params$species[i]))
    }
    sel <- params@w >= species_params$cutoff_size[i]
    if (!is.na(species_params$abundance_observed[i])) {
        abundance <- sum(params@initial_n[i, sel] * params@dw[sel])
        params@initial_n[i, ] <- params@initial_n[i, ] / abundance *
            species_params$abundance_observed[i]
    }
    if (!is.na(species_params$biomass_observed[i])) {
        biomass <- sum(params@initial_n[i, sel] * params@w[sel] * params@dw[sel])
        params@initial_n[i, ] <- params@initial_n[i, ] / biomass *
            species_params$biomass_observed[i]
    }
    
}

params <- steady(params, t_max = 100)
plotlySpectra(params)

saveRDS(params, file = "inst/tuning/params.rds")

# Slowly switch off background
params@species_params$rho_detritus[species_params$rho_detritus == 0] <- 
    (1 - epsilon) * (f0 / (1 - f0) * params@species_params$h)[species_params$rho_detritus == 0]
params <- setResourceEncounter(params)

plotlySpectra(params)
# Run to steady state with constant reproduction
(rdd <- getRDD(params))
oldsrr <- params@srr
params@srr <- function(rdi, species_params) {rdd}
sim <- project(params, t_max = 100)
plotBiomass(sim)
params@srr <- oldsrr
plotlySpectra(sim)

