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
# current L values consider, which generating the sucsession of mizer systems. In the algorithms here 
# I keep the step size constant, but what we could do, if we are able to sucsessfully find 
# a steady state over the last step of our procedure, then we can multiply the step size by a `ramp up` 
# constant bigger than 1 (but not making the step size too big to overshoot the destination L=1), wheras if the procedure fails to find the next steady state in the 
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


#########################

new_min_idx <- function(w,w_min){
    res <- w_min
    for (i in 1:length(res)){
        res[i] <- sum(w<=w_min[i])
    }
    return(res)
}


morph <- function(P,Q,T=50){
    r <- P
    for (t in 1:T){
        L <- t/T
        new_sp_data <- data.frame(
            species = Q@species_params$species,
            w_min = L*Q@species_params$w_min+(1-L)*P@species_params$w_min,
            w_inf = L*Q@species_params$w_inf+(1-L)*P@species_params$w_inf,
            w_mat = L*Q@species_params$w_mat+(1-L)*P@species_params$w_mat,
            w_min_idx = new_min_idx(Q@w,(L*Q@species_params$w_min+(1-L)*P@species_params$w_min)),
            h = L*Q@species_params$h+(1-L)*P@species_params$h,
            ks = L*Q@species_params$ks+(1-L)*P@species_params$ks,
            beta = L*Q@species_params$beta+(1-L)*P@species_params$beta,
            sigma = L*Q@species_params$sigma+(1-L)*P@species_params$sigma,
            z0 = L*Q@species_params$z0+(1-L)*P@species_params$z0,
            alpha = L*Q@species_params$alpha+(1-L)*P@species_params$alpha,
            erepro = L*Q@species_params$erepro+(1-L)*P@species_params$erepro,
            sel_func = Q@species_params$sel_func,
            # not used but required
            knife_edge_size = L*Q@species_params$knife_edge_size+(1-L)*P@species_params$knife_edge_size,
            gear = Q@species_params$gear,
            m = L*Q@species_params$m+(1-L)*P@species_params$m,
            w25 = L*Q@species_params$w25+(1-L)*P@species_params$w25,
            k = L*Q@species_params$k+(1-L)*P@species_params$k,
            gamma = L*Q@species_params$gamma+(1-L)*P@species_params$gamma,
            r_max = L*Q@species_params$r_max+(1-L)*P@species_params$r_max,
            linetype = Q@species_params$linetype,
            linecolour = Q@species_params$linecolour,
            l25 = L*Q@species_params$l25+(1-L)*P@species_params$l25,
            l50 = L*Q@species_params$l50+(1-L)*P@species_params$l50,
            k_vb = L*Q@species_params$k_vb+(1-L)*P@species_params$k_vb,
            a = L*Q@species_params$a+(1-L)*P@species_params$a,
            b = L*Q@species_params$b+(1-L)*P@species_params$b
        )
        
        rr <- multispeciesParams(object = new_sp_data,
                                 interaction = Q@interaction*L+ P@interaction*(1-L),
                                 min_w = min(Q@w)*L+(1-L)*min(P@w),
                                 max_w = max(Q@w)*L+(1-L)*max(P@w),
                                 no_w = length(P@w), # this better be the same in p and q
                                 min_w_pp = min(Q@w_full)*L+(1-L)*min(P@w_full),
                                 n = L*Q@n+(1-L)*P@n,
                                 p = L*Q@p+(1-L)*P@p,
                                 q = L*Q@p+(1-L)*P@q,
                                 r_pp = L*(Q@rr_pp / (Q@w_full ^ (Q@p - 1)))[1]+(1-L)*(P@rr_pp / (P@w_full ^ (P@p - 1)))[1],
                                 kappa = L*Q@kappa + (1-L)*P@kappa,
                                 lambda = L*Q@lambda + (1-L)*P@lambda,
                                 w_pp_cutoff = max(Q@w), # not sure about how to set this, it would be nice if it were passed along
                                 f0 = L*Q@f0+(1-L)*P@f0 # ,
                                 # z0pre = 0.6, z0exp = n - 1  # these arnt mentioned in addSpecies either
        )
        
        rr@initial_n <- r@initial_n
        rr@initial_n_pp <- r@initial_n_pp
        
        
        r <- steady(rr, effort = effort, t_max = 500,  tol = 1e-2) 
    }
    return(r)
}
########### testing it ############

PP <- p
QQ <- PP
QQ@interaction[] <- ((1:length(PP@interaction) %% 17)+1)/17 # psuedo random interaction matrix, for testing
RR <- morph(PP,QQ,50)
sim <- project(RR, t_max = 15, t_save = 0.1, effort = effort)
plot(sim)

# wrote a draft of the morph procedure so far. I want to debug it, and introduce variable step size.

# take apart and reconstruct and reconstruct a params object
# write morph procedure, but initiallly just changing theta
# gradually introduce other global variables into morph, checking it runs as we add flexability
# also allow dataframes to change
# introduce variable step size by using the following code structure for morph instead

r <- p
T <- 50
for (t in 1:T){
    rr <- r
    rr@interaction <- q@interaction*(t/T)+ p@interaction*(1-t/T)
    r <- steady(rr, effort = effort, t_max = 500,  tol = 1e-2) 
}
sim <- project(r, t_max = 15, t_save = 0.1, effort = effort)
plot(sim)


r <- p
new_theta <- r@interaction
new_theta[]<-((1:length(r@interaction) %% 17)+1)/17 # psuedo random interaction matrix, for testing
q <- p
q@interaction <- new_theta

max_fails <- 100
fail_count <- 0
T <- 50
dL <- 1/50
L <- 0
ramp_up <- 1.1
ramp_down <- 0.99
while (L<1){
    L_try <- L + dL
    rr <- r
    rr@interaction <- q@interaction*(t/T)+ p@interaction*(1-t/T)
    r_try <- steady(rr, effort = effort, t_max = 500,  tol = 1e-2)
    if (steady was successful) {
        r <- r_try
        L <- L_try
        dL <- max(c(dL*ramp_up,1-L))
    } else {
        fail_count <- fail_count+1
        if (fail_count>max_fails){
            break
        }
        dL <- dL*ramp_down
    }
}
# r is our objective

# written procedure structure for the morph code with variable step size, but have to figure 
# out how to detect and get a numerical trigger from when steady() fails to converge.