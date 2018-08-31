

# assuming SSB_observed[k] is the value for the kth foreground species, which corresponds with the kth 
# entry in the following list
#foreground_indices <- (1:length(p@A))[!is.na(p@A)]
# run humboldt first

reweight_for_SSB <- function(p,SSB){
    new_n <- p@initial_n
    for (i in (1:length(p@species_params$species))){
        if (!is.na(SSB[i])) {
            current_SSB <- sum(p@initial_n[i,] * p@w * p@dw * p@psi[i, ])
            new_n[i,] <- p@initial_n[i,]*SSB[i]/current_SSB
        }
    }
    p@initial_n <- new_n
    return(p)
}

get_current_SSB <- function(p){
    res <- p@species_params$w_mat
    for (i in (1:length(p@species_params$species))){
        res[i] <- sum(p@initial_n[i,] * p@w * p@dw * p@psi[i, ])
    }
    return(res)
}


fix_SSB <- function(p,SSB,trials=10, effort = 0, t_max = 100, t_save = 1,
                    tol = 1e-2, dt =0.1) {
    q <- p
    for (t in (1:trials)){
        q <- reweight_for_SSB(q,SSB)
        # q <- steady(q, effort = c(knife_edge_gear = 0, sigmoid_gear = input$effort, sigmoid_gear_Anchovy = input$Anchovy_effort), 
        #            t_max = 100, tol = 1e-2,
        #            shiny_progress = progress)
        q <- steady(q, effort = effort, 
                    t_max = t_max, tol = tol, t_per = t_max, t_save = t_save, dt = dt)
        print(get_current_SSB(q))
    }
    return(q)
}



fake_obs_SSB <- get_current_SSB(p)
fake_obs_SSB[is.na(p@A)] <- NA
foreground_indices <- (1:length(p@A))[!is.na(p@A)]
# double SSB of 1st foreground species, to make an obs_SSB for testing
#fake_obs_SSB[foreground_indices[1]] <- 2* fake_obs_SSB[foreground_indices[1]]

fake_obs_SSB[foreground_indices] <- 2* fake_obs_SSB[foreground_indices]


fake_obs_SSB
get_current_SSB(p)
p_fixed <- fix_SSB(p=p,SSB=fake_obs_SSB,trials=1, effort = effort, t_max = 10, 
                   t_save = 0.1, tol = 1e-2, dt = 0.01)
get_current_SSB(p_fixed)    

mult_rdd <- rep(1, 14)
mult_rdd[foreground_indices[1]] <- 2
p_m <- steady(p, mult_rdd = mult_rdd, t_max = 10,t_save = 0.1, tol = 1e-2,
              t_per = 10)

# p <- steady(p, effort = effort, t_max = 500,  tol = 1e-3)
# 
# 
# humboldt_params@species_params$SSB[] <- NA
# # this is a total hack, and I should have added the SSB into the underlying object, 
# # and I should not reply upon the ordering of the species.
# humboldt_params@species_params$SSB[(no_sp-7): no_sp] <- c(0.0186762374072649,
#                                                           0.000720757184972784,
#                                                           7.97858473029264e-05,
#                                                           2.51194830974027e-05,
#                                                           7.70690253424067e-06,
#                                                           0.000645875410082031,
#                                                           0.00023397243034225,
#                                                           0.000131377499657936)*0.1