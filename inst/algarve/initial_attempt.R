# Our first attempt at setting up a mizer model for the lagoon
species_params <- read.csv("inst/algarve/parameters.csv")[1:24, ]
names(species_params)[names(species_params) == "Species"] <- "species"

interaction <- read.csv("inst/algarve/interation_matrix.csv", row.names = 1)
interaction <- as.matrix(interaction)

resource_dynamics <- list("detritus" = function(params, n, n_pp, B, rates, dt, ...) B["detritus"])

params <- MizerParams(species_params,
                      interaction = interaction,
                      rho = rho,
                      resource_dynamics = resource_dynamics)

params@initial_B <- 10^8
names(params@initial_B) <- "detritus"
plotSpectra(params, plankton = FALSE)

# Run to steady state with constant reproduction
rdd <- getRDD(params)
params@srr <- function(rdi, species_params) {rdd}
sim <- project(params, t_max = 100)
plotBiomass(sim)
