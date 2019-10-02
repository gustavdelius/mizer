sp <- params@species_params
sp$ks <- NULL
lambda <- 2.12
min_egg <- min(sp$w_min)
max_w_inf <- max(sp$w_inf)
p <- markBackground(
    set_scaling_model(min_w_pp = 1e-10,
                      no_sp = 10, no_w = 400, min_w_inf = 2, max_w_inf = max_w_inf,
                      min_egg = min_egg, min_w_mat = 2 / 10^0.8, 
                      w_pp_cutoff = 2 / 10^0.8, perfect = FALSE,
                      lambda = lambda, beta = 60, sigma = 2,
                      knife_edge_size = Inf,
                      n = 0.7, p = 0.7)
)
p <- steady(p)
plotSpectra(p, power = 2)

ps <- addSpecies(p, sp)

sp$beta[1] <- 60
