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
#sim <- project(p, t_max = 15, t_save = 0.1, effort = effort)
#plot(sim)

#dim(sim@n)

nn <- p@initial_n
nn[1, ] <- 10*nn[1,]
t_max <- 15
t_save <- 0.1
sim <- project(p, t_max = t_max, t_save = t_save, effort = effort, initial_n = nn)
plot(sim)

#sim@n[,1,]

library("plot3D")

multi <- function(x,y){x*y}

biomass_spectra_dynamics <- sweep(sim@n[,1,],2,p@w,multi)
plot(p@w,biomass_spectra_dynamics[1,],log="xy")
#persp(biomass_spectra_dynamics)
my_times <- seq(from = 0, to = t_max, by = t_save)
#contour2D(z=biomass_spectra_dynamics,y=log(p@w),x=my_times)
#persp(z=log(biomass_spectra_dynamics),y=log(p@w),x=my_times,
#      zlim = c(log(min(biomass_spectra_dynamics[biomass_spectra_dynamics>0])),max(log(biomass_spectra_dynamics))))
#persp3D(z=log(biomass_spectra_dynamics),y=log(p@w),x=my_times,
#      zlim = c(log(min(biomass_spectra_dynamics[biomass_spectra_dynamics>0])),max(log(biomass_spectra_dynamics))))

#L <- log(min(biomass_spectra_dynamics[biomass_spectra_dynamics>0]))

truncated_log_biomass_spectra_dynamics <- log(biomass_spectra_dynamics)
truncated_log_biomass_spectra_dynamics[biomass_spectra_dynamics==0] <- min(log(biomass_spectra_dynamics)[biomass_spectra_dynamics>0])-1
contour2D(z=truncated_log_biomass_spectra_dynamics,y=log(p@w),x=my_times,xlab="time", ylab="log(weight)")

#made perp3D plot of biomass density of species 1 in the humboldt system changing over time
persp3D(z=truncated_log_biomass_spectra_dynamics,y=log(p@w),x=my_times,xlab="time", ylab="log(weight)", zlab = "log(biomass_density",zlim=c(min(truncated_log_biomass_spectra_dynamics),-0.1),phi = 40, theta = 40)


# finding new ways to plot size spectrum dynamics
#z = log(biomass_density), for species 1 in humboldt.R (lines 1..60) system, 
#starting from a state like the steady state, but with the initial abundance of
#species 1 multiplied by 10.  
# after this we should look at ways of visualizing biomass flow, and influences 
# of species on one anothers death and growth using digraphs with edge thickness 
# corresponding to interaction levels. Or one could just plot the k strongest links. 
# For example, if E_ij is the biomass of species i that is consumed by a member 
# of species j at maturity weight, and k >0 then we could consider the directed graph 
# where there is a link from i to j if E_i.j is one of the k largest entries in the set 
# E_1.j, E_2.j,..,E_no_sp.j  (i.e. if species i is amoungst the k most influencial species upon
# the amount of biomass obtained by mature members of species j). This should give a way to 
# visualize which species effect one another, perhaps a bit like protien-protien interaction networks. 
# To follow the analogy further, I guess possitive biomass flow is like an `activiator` link 
# but somehow we may also want to visualize the most important `inhibitter` links that 
# correspond to death effects, but I think it might be cleaner to visualize the 
# death and biomass flow digraphs seperately.

########### now for similar plots using community abundance

NC <- sim@n[,1,]
NC[] <- 0
for (i in 1:dim(sim@n)[2]){
    NC <- NC + sim@n[,i,]
}

biomass_spectra_dynamics_aggregate <- sweep(NC,2,p@w,multi)

truncated_log_biomass_spectra_dynamics_aggregate <- log(biomass_spectra_dynamics_aggregate)
truncated_log_biomass_spectra_dynamics_aggregate[biomass_spectra_dynamics_aggregate==0] <- min(log(biomass_spectra_dynamics_aggregate)[biomass_spectra_dynamics_aggregate>0])-1

# I made similar plots for community abundance. It would be good if we could colour the 
# different parts of the surface according to the species which has the largest biomass 
# density at a particular time.
contour2D(z=truncated_log_biomass_spectra_dynamics_aggregate,y=log(p@w),x=my_times,xlab="time", ylab="log(weight)")

#made perp3D plot of biomass density of the aggregation of all species over time
persp3D(z=truncated_log_biomass_spectra_dynamics_aggregate,y=log(p@w),x=my_times,xlab="time", ylab="log(weight)", zlab = "log(biomass_density",zlim=c(min(truncated_log_biomass_spectra_dynamics),-0.1),phi = 35, theta = 130)

