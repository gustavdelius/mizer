source("./inst/gustav/tragedy_fns.R")

# Find optimal knife_edge_size for the trait_based model

knife_edge_size <- 1000
effort <- 0.5
t_max <- 100
params_knife <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = knife_edge_size)
sim <- project(params_knife, effort = effort, t_max = t_max)
plot(sim)
final_n <- sim@n[dim(sim@n)[1], , ] 
final_n_pp <- sim@n_pp[dim(sim@n_pp)[1], ]

#' Calculate yield from given mesh size
yield <- function(size, effort, t_max, initial_n, initial_n_pp) {
    params <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = size)
    sim <<- project(params, effort = effort, t_max = t_max, initial_n = initial_n, initial_n_pp = initial_n_pp)
    n_total <- apply(sim@n, c(1,3), sum)
    biomass_total <- sweep(n_total,2,sim@params@w * sim@params@dw, "*")
    yield_total <- rowSums(biomass_total[ ,sim@params@w > size]) * effort
    plot(yield_total, type="l", ylim=c(0.99*min(yield_total),1.01*max(yield_total)))
    y <- yield_total[length(yield_total)]
    print(y)
    return(list(y=y, final_n=sim@n[dim(sim@n)[1], , ],final_n_pp=sim@n_pp[dim(sim@n_pp)[1], ]))
}

t_max<-50
l <- yield(1000, effort, t_max, initial_n=final_n, initial_n_pp=final_n_pp)
for (size in seq(from = 190, to = 150, by = -10)) {
    l <- yield(size, effort, t_max, initial_n=l$final_n, initial_n_pp=l$final_n_pp)
}
