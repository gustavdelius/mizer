params_data <- read.csv(system.file("extdata", "speciesNCME_edited2.csv", package = "mizer"))
effort <- 1.4

no_sp <- dim(params_data)[1]

l25 <- c(1.0e+29,     1.9e+00,     4.0e+00,     5.0e+00,     2.9e+00,     3.2e+01,     4.9e+00,     4.9e+01 )
l50 <-  c(1.1e+29,     2.0e+00,     5.0e+00,     6.0e+00,     3.0e+00,     3.6e+01,     8.0e+00,     5.1e+01) 
names(l25) <- as.character(params_data$species)
names(l50) <- as.character(params_data$species)

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
    
    p <- addSpecies(p, species_params, effort = effort, rfac=Inf)
}

# Run to steady state
p <- steady(p, effort = effort, t_max = 500,  tol = 1e-3)
sim <- project(p, t_max = 15, t_save = 0.1, effort = effort)
plot(sim)

############## varying effort demonstration
r <- p
effort2 <- 2*effort
T <- 5
for (t in 1:T){
    X <- effort2*(t/T)+effort*(1-t/T)
    r <- steady(r, effort = X, t_max = 500,  tol = 1e-2)
}
sim <- project(r, t_max = 15, t_save = 0.1, effort = effort)
plot(sim)



################ Trying out the morph idea with code that takes the current homboldt system
# P at steady state, with its interaction matrix full of ones, and then considers another 
# system Q, which is just the same as P, except the interaction matrix of Q is full of psudorandomly generated 
# real numbers between 0 and 1. Then we use the morph idea, gradually changing from P to Q, by 
# gradually changing the interaction matrix, and keeping track of the steady state as we do this.
# Next I want to design an algorithm that does the deformation with variable step size, as follows: 
# The progress from P to Q can be kept track of by a fraction 0<=L<=1, meaning we are looking 
# for a steady state of LQ+(1-L)P, and the step size s is the difference between the next and 
# current L values consider, which generating the sucsession of sequences. In the algorithms here 
# I keep the step size constant, but what we could do, if we are able to sucsessfully find 
# a steady state over the last step of our procedure, then we can multiply the step size by a `ramp up` 
# constant bigger than 1, wheras if the procedure fails to find the next steady state in the 
# sequence, that means that the step size is too large, and so we multiply the step size by a 
# `ramp down` constant, less than one. Also, we keep count of how many failures we have, 
# and if they exceed a given threshold we halt. We could also include newton raphson, then 
# I guess as long as the steady state remains present and changes continuously throughout deformation,
# from P to Q the morph procedure should always 
# work in principal, in the ideal continuous case with unlimited run time.


new_theta <- p@interaction
new_theta[]<-((1:length( p@interaction) %% 17)+1)/17 # psuedo random interaction matrix, for testing
q <- p
q@interaction <- new_theta
r <- p
T <- 50
for (t in 1:T){
    rr <- r
    rr@interaction <- q@interaction*(t/T)+ p@interaction*(1-t/T)
    r <- steady(rr, effort = effort, t_max = 500,  tol = 1e-2) 
}
sim <- project(r, t_max = 15, t_save = 0.1, effort = effort)
plot(sim)

##########################


#new_theta <- p@interaction
#new_theta[]<-runif(length( p@interaction))



# To illustraite this morph idea, in a simple way, I got the humboldt system at steady state P, 
# and I let Q be the same system, but with twice the fishing effort, and then I gradually 
# changed the system from P to Q, and was able to use steady to sucsessfully track the steady state. 
# 