params_data <- read.csv(system.file("extdata", "speciesNCME_edited2.csv", package = "mizer"))
effort <- 1.4

no_sp <- dim(params_data)[1]

l25 <- c(1.0e+29,     1.9e+00,     4.0e+00,     5.0e+00,     2.9e+00,     3.2e+01,     4.9e+00,     4.9e+01 )
l50 <-  c(1.1e+29,     2.0e+00,     5.0e+00,     6.0e+00,     3.0e+00,     3.6e+01,     8.0e+00,     5.1e+01) 
names(l25) <- as.character(params_data$species)
names(l50) <- as.character(params_data$species)

###

no_pts <- 10
#foreground_multiplier <- 10^(-(1:no_pts))
foreground_multiplier <- (1:no_pts)*0.1/no_pts

deviation_from_power_law <- 1:no_pts
for (kk in (1:no_pts)){
    
#foreground_mult <- 10^(-6) # in usage in app, we used 0.1
    foreground_mult <- foreground_multiplier[kk]
    
# 
SSBvals <- c(0.0186762374072649,
                                                                     0.000720757184972784,
                                                                     7.97858473029264e-05,
                                                                     2.51194830974027e-05,
                                                                     7.70690253424067e-06,
                                                                     0.000645875410082031,
                                                                     0.00023397243034225,
                                                                     0.000131377499657936)*foreground_mult

###

p <- setBackground(
    set_scaling_model(min_w_pp = 1e-12,
                      no_sp = 10, no_w = 400, min_w_inf = 2, max_w_inf = 6e5,
                      min_egg = 1e-4, min_w_mat = 2 / 10^0.6, 
                      lambda = 2.12,
                      knife_edge_size = Inf)
)


all_efforts <- c(0, 1.4, 1.1)
names(all_efforts) <- c("knife_edge_gear", "sigmoid_gear", "sigmoid_gear_Anchovy")
effort <- all_efforts[1:2]
for (i in (1:no_sp)) {
    if (params_data$species[i] == "Anchovy") {
        effort <- c(effort, all_efforts[3])
    }
    a_m <- params_data$a2[i]
    b_m <- params_data$b2[i]
    L_inf_m <- params_data$Linf[i]
    L_mat <- params_data$Lmat[i]
    if (params_data$species[i] == "Anchovy") {
        gear <- "sigmoid_gear_Anchovy"
    } else {
        gear <- "sigmoid_gear"
    }
    species_params <- data.frame(
        species = as.character(params_data$species[i]),
        w_min = params_data$Wegg[i],
        w_inf = params_data$w_inf[i],
        w_mat = params_data$w_mat[i],
        beta = params_data$beta[i],
        sigma = log(params_data$sigma[i]),
        z0 = 0,
        alpha = 0.6,
        erepro = 0.1, # unknown, determined later
        sel_func = "sigmoid_length",
        gear = gear,
        l25 = l25[i],
        l50 = l50[i],
        k = 0,
        k_vb = params_data$k_vb[i],
        a = a_m,
        b = b_m
    )
    
    p <- addSpecies(p, species_params, effort = effort, rfac=Inf, SSB = SSBvals[i])
}

# Run to steady state
p <- steady(p, effort = effort, t_max = 500,  tol = 1e-2)
# sim <- project(p, t_max = 15, t_save = 0.1, effort = effort)
# plot(sim)
# fg <- (7:14)
# 
# resulting_SSB <- SSBvals
# 
# counter <- 0
# for (i in fg){
#     counter <- counter + 1
#     resulting_SSB[counter] <- sum(sim@n[dim(sim@n)[1],i,] *sim@params@w * sim@params@dw * sim@params@psi[i, ])
# }
# 
# absolute_error <- sum((SSBvals-resulting_SSB)^2)
# relative_error <- sum(((SSBvals-resulting_SSB)/SSBvals)^2)
# SSBvals
# resulting_SSB

#community_abundance <- colSums(sim@n[dim(sim@n)[1],,])+sim@n_pp[dim(sim@n_pp)[1],(sim@params@w_full>=min(sim@params@w))]

community_abundance <- colSums(p@initial_n)+p@initial_n_pp[p@w_full>=min(p@w)]

#plot(sim@params@w,community_abundance,log="xy")

#power_law <- sim@params@kappa*sim@params@w^(-sim@params@lambda)

power_law <- p@kappa*p@w^(-p@lambda)


power_law_deviation <- sum(((community_abundance-power_law)/power_law)^2)

deviation_from_power_law[kk] <- power_law_deviation
}
plot(foreground_multiplier, deviation_from_power_law, type="l")


#############


# 
# SSB_table_custom <- function(sim){
#     no_sp <- dim(sim@params@psi)[1]
#     trueres <- matrix(0,nrow = no_sp,ncol=3)
#     for (i in (1:no_sp)){
#         trueres[i,3] <- sum(sim@n[dim(sim@n)[1],i,] *sim@params@w * sim@params@dw * sim@params@psi[i, ])
#     }
#     trueres[,1] <- names(sim@params@psi[,1])
#     trueres[,2] <- trueres[,1]
#     trueres[fg,2] <- SSBvals
#     colnames(trueres) <- c("Species","Target SSB","Final SSB")
#     rownames(trueres) <- names(sim@params@psi[,1])
#     return(trueres)
# }
# 
# SSB_tab <- SSB_table_custom(sim)
# absolute_error <- sum((SSB_tab[fg,2]-SSB_tab[fg,3])^2)
# relative_error <- sum(((SSB_tab[fg,2]-SSB_tab[fg,3])/SSB_tab[fg,2])^2)





# declare SSB in here, using same values as at the start of app_functional, and rerun, 
# investigating the behaivours as we alter the fraction of SSB taken up by the foreground species
# and investigate the deviations between the target, and actual growth rates. 

# we have computed the relative error in resulting SSB in this old humboldt system, and this is ready to wrap in a loop 
# to investigate the accuracy of the model, for different levels of abundance of the foreground species. 
# We can plot how the relative error in the SSB, as well as the distance from the steady state change 
# with the amount of foreground species. 

# To do next:
#
# next, add computations of deviation from the powerlaw, (could also see sum of squared error of abundances, 
# over all weights and species, and relative version of this)
# next write loop and make plots,
# next add a global foreground abundance multiplier to the app, and investigate it there too,
# next try and make a version of humboldt that acts like an inert version of the app 