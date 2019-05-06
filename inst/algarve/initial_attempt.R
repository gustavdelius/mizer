# Our first attempt at setting up a mizer model for the lagoon
species_params <- read.csv("inst/algarve/parameters.csv")[1:24, ]
names(species_params)[names(species_params) == "Species"] <- "species"
<<<<<<< HEAD
species_params$beta[is.na(species_params$beta)] <- 100
species_params$sigma[is.na(species_params$sigma)] <- 2
species_params$w_mat[species_params$w_mat == 0] <- NA
species_params$r_max <- Inf

species_params$ks <- 0
=======
>>>>>>> bd6456e706404d0e13937c821fe8a68d4ef7a5ae

interaction <- read.csv("inst/algarve/interation_matrix.csv", row.names = 1)
interaction <- as.matrix(interaction)


resource_params <- list("detritus_external" = 0)

resource_dynamics <- list("detritus" = function(params, n, n_pp, B, rates, dt, ...) B["detritus"])

params <- MizerParams(species_params,
                      interaction = interaction,
                      resource_dynamics = resource_dynamics,
                      resource_params = resource_params)

params@initial_B[] <- 10^8
plotlySpectra(params, plankton = FALSE)

# Run to steady state with constant reproduction
(rdd <- getRDD(params))
params <- steady(params, t_max = 100)
plotlySpectra(params)

# Play with allometric death rate
mu_b0 <- rep(1, nrow(params@species_params))
mu_b <- outer(mu_b0, params@w^(-0.2))
params <- setBMort(params, mu_b = mu_b)
params <- steady(params, t_max = 100)
plotlySpectra(params)

# Calculate total biomass
biomass <- rowSums(sweep(params@initial_n, 2, params@w * params@dw, "*"))

# rescale abundances
factor <- 1e-14
params@initial_n <- params@initial_n * factor
params@initial_n_pp <- params@initial_n_pp * factor
params@initial_B <- params@initial_B * factor
params@cc_pp <- params@cc_pp * factor
params@species_params$gamma <- params@species_params$gamma / factor
params <- setSearchVolume(params)
params@rho <- params@rho / factor
params <- steady(params)

