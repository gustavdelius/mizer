library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library(mizer)
params_data <- read.csv("./vignettes/NS_species_params.csv")
inter <- read.csv("./vignettes/inter.csv", row.names=1)
inter <- as(inter, "matrix")
NS_params <- MizerParams(params_data)


no_traits <- 10
params <- set_trait_model(no_sp = no_traits+1, min_w_inf = 10, max_w_inf = 1e6)
#plot(project(params, t_max=16, effort = 0))
kappa_def <- 0.005
k0_def <- 50
#scalec <- 10^11
#scalec <- 54000*10^9
scalec <- 10^13
paramsK <- set_trait_model(no_sp = no_traits+1,
                           min_w_inf = 10,
                           max_w_inf = 1e6,kappa=scalec*kappa_def,
                           k0=scalec*k0_def)
#plot(project(paramsK, t_max=16, effort = 0))

#############

params2 <- paramsK
NS_params <- MizerParams(params_data,min_w_pp=min(params@w_full), min_w=min(params@w),max_w=max(params@w),no_w=100)
#all(params@w_full==NS_params@w_full)


#rp <- 5
rp <- no_traits+1
params2@species_params$w_inf[rp] <- NS_params@species_params$w_inf[6]
params2@species_params$w_mat[rp] <- NS_params@species_params$w_mat[6]
params2@species_params$h[rp] <- NS_params@species_params$h[6]
params2@species_params$gamma[rp] <- NS_params@species_params$gamma[6]
params2@species_params$ks[rp] <- NS_params@species_params$ks[6]
params2@species_params$beta[rp] <- NS_params@species_params$beta[6]
params2@species_params$sigma[rp] <- NS_params@species_params$sigma[6]
params2@species_params$z0[rp] <- NS_params@species_params$z0[6]
params2@species_params$alpha[rp] <- NS_params@species_params$alpha[6]
params2@species_params$sel_func[rp] <- NS_params@species_params$sel_func[6]
params2@species_params$knife_edge_size[rp] <- NS_params@species_params$knife_edge_size[6]
params2@species_params$erepro[rp] <- NS_params@species_params$erepro[6]
params2@species_params$k[rp] <- NS_params@species_params$k[6]
params2@species_params$w_min[rp] <- NS_params@species_params$w_min[6]
params2@species_params$w_min_idx[rp] <- NS_params@species_params$w_min_idx[6]
params2@species_params$r_max[rp] <- NS_params@species_params$r_max[6]
params2@psi[rp,] <- NS_params@psi[6,] 
params2@intake_max[rp,] <- NS_params@intake_max[6,] 
params2@search_vol[rp,] <- NS_params@search_vol[6,] 
params2@activity[rp,] <- NS_params@activity[6,] 
params2@std_metab[rp,] <- NS_params@std_metab[6,] 
params2@ft_pred_kernel_e[rp,] <- NS_params@ft_pred_kernel_e[6,] 

sim2 <- project(params2, t_max=100, effort = 0)
plot(sim2)

plot(params2@w,sim2@n[50,rp,],log="xy")
