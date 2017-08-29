library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library(mizer)

params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")

#get north sea parameters
#NS_params <- MizerParams(params_data, interaction = inter)
NS_params <- MizerParams(params_data)

# make a system with no_traits plus whiting

no_traits <- 10

params <- set_trait_model(no_sp = no_traits+1, min_w_inf = 10, max_w_inf = 1e5)

#params@species_params[no_traits+1,] <- NS_params@species_params[6,]


sim <- project(params, t_max=76, effort = 0)
plot(sim)

params@species_params$species
# add mackrel


params2 <- params
params2@species_params$w_inf[no_traits+1] <- NS_params@species_params$w_inf[6]
params2@species_params$w_mat[no_traits+1] <- NS_params@species_params$w_mat[6]
params2@species_params$h[no_traits+1] <- NS_params@species_params$h[6]
params2@species_params$gamma[no_traits+1] <- NS_params@species_params$gamma[6]
params2@species_params$ks[no_traits+1] <- NS_params@species_params$ks[6]
params2@species_params$beta[no_traits+1] <- NS_params@species_params$beta[6]
params2@species_params$sigma[no_traits+1] <- NS_params@species_params$sigma[6]
params2@species_params$z0[no_traits+1] <- NS_params@species_params$z0[6]
params2@species_params$alpha[no_traits+1] <- NS_params@species_params$alpha[6]
params2@species_params$sel_func[no_traits+1] <- NS_params@species_params$sel_func[6]
params2@species_params$knife_edge_size[no_traits+1] <- NS_params@species_params$knife_edge_size[6]
params2@species_params$gear[no_traits+1] <- NS_params@species_params$gear[6]
params2@species_params$erepro[no_traits+1] <- NS_params@species_params$erepro[6]
params2@species_params$k[no_traits+1] <- NS_params@species_params$k[6]
params2@species_params$w_min[no_traits+1] <- NS_params@species_params$w_min[6]
params2@species_params$w_min_idx[no_traits+1] <- NS_params@species_params$w_min_idx[6]
params2@species_params$r_max[no_traits+1] <- NS_params@species_params$r_max[6]

sim <- project(params2, t_max=7, effort = 0)
plot(sim)
