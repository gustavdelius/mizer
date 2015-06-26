source("./inst/gustav/tragedy_fns.R")

# Run the trait_based model
knife_edge_size <- 1000
effort <- 0.5
t_max <- 50
params_knife <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = knife_edge_size)
sim <- project(params_knife, effort = effort, t_max = t_max)
plot(sim)
# Plot yield
plot_yield(sim)

# Rerun with different knife_edge_size to find that size that optimises the yield
knife_edge_size <- 0.13
t_max <- 10
params_knife <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = knife_edge_size)
sim <- project(params_knife, effort = effort, t_max = t_max, initial_n = sim@n[dim(sim@n)[1], , ], initial_n_pp = sim@n_pp[dim(sim@n_pp)[1], ])
c1<-plot_yield(sim); c1

plotSpectra(sim)

# Simulation with exploiting fisher
tragedy_size <- 0.9*knife_edge_size
n_fishers <- 5
params_tragedy<- set_tragedy_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = knife_edge_size,
                                  n_fishers = n_fishers, tragedy_size = tragedy_size)
sim_tragedy <- project(params_tragedy, effort = effort, t_max = t_max, initial_n = sim@n[dim(sim@n)[1], , ], initial_n_pp = sim@n_pp[dim(sim@n_pp)[1], ])
e1<-plot_yield(sim_tragedy); e1
# Note that his yield higher than before, so the fisher will decrease the net size

# Eventually all other fishers lower their mesh size as well
knife_edge_size <- tragedy_size
params_knife <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = knife_edge_size)
sim <- project(params_knife, effort = effort, t_max = t_max, initial_n = sim@n[dim(sim@n)[1], , ], initial_n_pp = sim@n_pp[dim(sim@n_pp)[1], ])
c2<-plot_yield(sim);c2
# So now the situation is worse for everyone, but not by that much

# Now the exploiting fisher lowers his mesh size a second time
tragedy_size <- 0.9*knife_edge_size
params_tragedy<- set_tragedy_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = knife_edge_size,
                                   n_fishers = n_fishers, tragedy_size = tragedy_size)
sim_tragedy <- project(params_tragedy, effort = effort, t_max = t_max, initial_n = sim@n[dim(sim@n)[1], , ], initial_n_pp = sim@n_pp[dim(sim@n_pp)[1], ])
e2<-plot_yield(sim_tragedy);e2
# Which again improves his yield compared to the community yield

# The other fishers follow suit
knife_edge_size <- tragedy_size
params_knife <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = knife_edge_size)
sim <- project(params_knife, effort = effort, t_max = t_max, initial_n = sim@n[dim(sim@n)[1], , ], initial_n_pp = sim@n_pp[dim(sim@n_pp)[1], ])
c3<-plot_yield(sim);c3

# Now the exploiting fisher lowers his mesh size a third time
tragedy_size <- 0.9*knife_edge_size
params_tragedy<- set_tragedy_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = knife_edge_size,
                                   n_fishers = n_fishers, tragedy_size = tragedy_size)
sim_tragedy <- project(params_tragedy, effort = effort, t_max = t_max, initial_n = sim@n[dim(sim@n)[1], , ], initial_n_pp = sim@n_pp[dim(sim@n_pp)[1], ])
e3<-plot_yield(sim_tragedy);e3
# His yield no longer improves. So the mesh size will not decrease any further

# Thus the total yield has only decreased from
c1
# to
c3
# which is only a decrease of 
(c1-c3)/c1*100
# percent

