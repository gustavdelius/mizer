# Calculate average psi

# Set up a trait based model with 20 species
params <- set_trait_model(no_sp = 40, min_w_inf = 10, max_w_inf = 1e5)
sim <- project(params, t_max=75, effort = 0)
plot(sim)

# average over the first 20 species
b<-colSums(sim@n[76, 1:20, ])
a<-colSums(sim@n[76, 1:20, ]*sim@params@psi[1:20, ])
plot(a/b~sim@params@w, log="x")

# Draw in maturity of smallest and largest species
abline(v=sim@params@species_params$w_mat[1])
abline(v=sim@params@species_params$w_mat[20])
# and maximum size of smallest species
abline(v=sim@params@species_params$w_inf[1])
