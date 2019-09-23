library(tidyverse)
options("mizer_new" = TRUE)

# Load species params ####

theta <- read.csv("inst/blanes/theta.csv", header = FALSE)
theta <- t(as.matrix(theta))

sp <- read.csv("inst/blanes/species_params.csv")
no_sp <- nrow(sp)
sp$R_max <- Inf
sp$r_max <- NULL
sp$alpha <- 0.4

no_w <- 200
min_w <- min(sp$w_min)
# Extend background species a bit beyond largest foreground species
max_w <- max(sp$w_inf) * 4

min_w_mat = min_w * 10^5 / 20

# Create background ####
# also the smallest eggs should have something to eat
dx <- log10(max_w / min_w) / (no_w - 1)
ppmr <- 10^(seq(from = 0, by = dx, length.out = 3 * no_w))
phis <- get_phi(sp, ppmr)
max_ppmr <- apply(phis, 1, function(x) ppmr[max(which(x != 0)) + 1])
min_w_pp <- min(sp$w_min / max_ppmr)

p <- set_scaling_model(
    min_w_pp = NA,
    no_sp = 14, 
    min_w_inf = min_w * 10^5, 
    max_w_inf = max_w,
    min_egg = min_w,
    min_w_mat = min_w_mat,
    w_pp_cutoff = min_w_mat,
    no_w = no_w,
    n = 0.7,
    q = 0.8,
    h = 30,
    sigma = 2,
    r_pp = 2000,
    bmort_prop = 0.01)

p <- p %>% 
    steady() %>% 
    markBackground()

plotSpectra(p, power = 2, total = TRUE, 
            ylim = c(10^-6, NA), wlim = c(10^-8, NA))

# Add species as planktivores ####
sp$interaction_p <- 1
ps <- addSpecies(p, sp, interaction = theta) %>% 
    steady()

plotlySpectra(ps, power = 2, total = TRUE, 
              ylim = c(10^-6, NA), wlim = c(10^-8, NA))

# Attempt to fix very low growth rates ####
fix <- which(ps@species_params$h[15:39] < 80)
factor <- 80 / ps@species_params$h[15:39][fix]
sp$k_vb[fix] <- sp$k_vb[fix] * factor

pss <- addSpecies(p, sp, interaction = theta) 
pss <- rescaleAbundance(pss, 0.1)
pss <- steady(pss)
plotlySpectra(pss, power = 2, total = TRUE, 
              ylim = c(10^-6, NA), wlim = c(10^-8, NA))
pst <- tuneParams(pss)
saveRDS(pst, file = "inst/blanes/pst.rds")

pst <- readRDS("inst/blanes/pst.rds")
# Add resource components ####
# Set up params with constant resource biomass
resource_dynamics <- 
    list("detritus" = function(params, n, n_pp, B, rates, dt, ...) B["detritus"],
         "carrion" = function(params, n, n_pp, B, rates, dt, ...) B["carrion"])
resource_params <- list("detritus_external" = 0, #detritus_cons,
                          "detritus_proportion" = 0,
                          "carrion_external" = 0) #carrion_cons)
pst <- setResourceDynamics(pst, resource_dynamics = resource_dynamics,
                           resource_params = resource_params)

# default for rho
detQ <- pst@species_params$detQ
detQ[1:14] <- 0
scavQ <- pst@species_params$scavQ
scavQ[1:14] <- 0
pst@species_params$rho_detritus <- pst@species_params$h * 
    pst@f0 / (1 - pst@f0) * detQ
pst@species_params$rho_carrion <- pst@species_params$h * 
    pst@f0 / (1 - pst@f0) * scavQ
pst <- setResourceEncounter(pst)

# Switch off plankton for detritivores and scavengers ####
pst@species_params$interaction_p[detQ == 1 | scavQ == 1] <- 0

pst@initial_B <- c(detritus = 1, carrion = 1)

pstt <- steady(pst)
plotlySpectra(psttt, power = 2, total = TRUE, 
            ylim = c(10^-6, NA), wlim = c(10^-8, NA))

# Switch off plankton for everyone else by hand ####
psttt <- tuneParams(pstt)

psttt <- removeSpecies(psttt, 1:14)
psttt@cc_pp[] <- 0
psttt@initial_n_pp[] <- 0

psttt <- tuneParams(psttt)

# Reduce background mortality

psttt@mu_b <- psttt@mu_b * 1.1
p1 <- steady(psttt)
plotlySpectra(p1, power = 2, total = TRUE, 
              ylim = c(10^-6, NA), wlim = c(10^-8, NA))
p1 <- tuneParams(p1)

p2 <- setBMort(p1) %>% steady()
p2 <- tuneParams(p2)
