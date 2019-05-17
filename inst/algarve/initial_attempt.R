# Our first attempt at setting up a mizer model for the lagoon
species_params <- read.csv("inst/algarve/parameters.csv")[1:24, ]
names(species_params)[names(species_params) == "Species"] <- "species"
species_params$beta[is.na(species_params$beta)] <- 100
species_params$sigma[is.na(species_params$sigma)] <- 2
species_params$w_mat[species_params$w_mat == 0] <- NA
species_params$r_max <- Inf

species_params$cutoff_size <- species_params$minimum.sampling.size..g.
species_params$cutoff_size[24] <- species_params$w_inf[24] / 2
species_params$biomass_observed <- species_params$Biomass.in.habitat.area..gDWm.2.

interaction <- read.csv("inst/algarve/interation_matrix.csv", row.names = 1)
interaction <- as.matrix(interaction)

resource_params <- list("detritus_external" = 0)
resource_dynamics <- list("detritus" = function(params, n, n_pp, B, rates, dt, ...) B["detritus"])

f0 <- 0.6
kappa <- 1
n <- 2 / 3
q <- 3 / 4
lambda = 2 + q - n

no_sp <- length(species_params$species)

params <- MizerParams(species_params,
                      resource_dynamics = resource_dynamics,
                      resource_params = resource_params,
                      f0 = f0, kappa = kappa, n = n, q = q, lambda = lambda)

params@initial_B[] <- 1

# Start with no interaction
interaction0 <- interaction
interaction0[] <- 0
interaction_p <- params@species_params$interaction_p
params@species_params$interaction_p[] <- 0
params <- setInteraction(params, interaction = interaction0)

# Make everything feed on detritus
params@species_params$rho_detritus <- f0 / (1 - f0) * params@species_params$h
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
    if (!is.na(species_params$biomass_observed[i])) {
        biomass <- sum(params@initial_n[i, sel] * params@w[sel] * params@dw[sel])
        params@initial_n[i, ] <- params@initial_n[i, ] / biomass *
            species_params$biomass_observed[i]
    } else {
        stop("Missing biomass")
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
params <- setInteraction(params, interaction = interaction * epsilon)

params <- steady(params, t_max = 100)
plotlySpectra(params)

saveRDS(params, file = "inst/algarve/params.rds")

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

