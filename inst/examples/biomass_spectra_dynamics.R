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
library(plot3D)
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
library(rgl)
plot3d(z=truncated_log_biomass_spectra_dynamics_aggregate,y=log(p@w),x=my_times)


open3d()
persp3d(z=truncated_log_biomass_spectra_dynamics_aggregate,y=log(p@w),x=my_times, col = "blue")
#grid3d(c("x", "y+", "z"))

############## same as above, but using N_i(w)*w^2 instead of N_i(w)*w

biomass_spectra_dynamics_aggregate2 <- sweep(NC,2,p@w^2,multi)

truncated_log_biomass_spectra_dynamics_aggregate2 <- log(biomass_spectra_dynamics_aggregate2)
truncated_log_biomass_spectra_dynamics_aggregate2[biomass_spectra_dynamics_aggregate2==0] <- min(log(biomass_spectra_dynamics_aggregate2)[biomass_spectra_dynamics_aggregate2>0])-1
persp3d(z=truncated_log_biomass_spectra_dynamics_aggregate2,y=log(p@w),x=my_times, col = "blue")

############### next get multiple surfaces on one plot
# https://stackoverflow.com/questions/39948720/plot-multiple-surfaces
# #Create Plotly object
# plot_ly(showscale = FALSE) %>%
#     
#     #Volcano surface    
#     add_surface(z = volcano) %>%
#     
#     #First rectangle
#     add_surface(x = c(10, 60),
#                 y = c(10, 50),
#                 z = matrix(160, nrow = 2, ncol = 2)) %>%
#     
#     #Second rectangle
#     add_surface(x = c(10, 60),
#                 y = c(10, 50),
#                 z = matrix(180, nrow = 2, ncol = 2))

library(plotly)
plot_ly(showscale = FALSE)
    add_surface(z=truncated_log_biomass_spectra_dynamics_aggregate2,y=log(p@w),x=my_times, col = "blue")

    ######################
    
   # install.packages("magrittr")
    library(magrittr)
    
    z <- c(
        c(8.83,8.89,8.81,8.87,8.9,8.87),
        c(8.89,8.94,8.85,8.94,8.96,8.92),
        c(8.84,8.9,8.82,8.92,8.93,8.91),
        c(8.79,8.85,8.79,8.9,8.94,8.92),
        c(8.79,8.88,8.81,8.9,8.95,8.92),
        c(8.8,8.82,8.78,8.91,8.94,8.92),
        c(8.75,8.78,8.77,8.91,8.95,8.92),
        c(8.8,8.8,8.77,8.91,8.95,8.94),
        c(8.74,8.81,8.76,8.93,8.98,8.99),
        c(8.89,8.99,8.92,9.1,9.13,9.11),
        c(8.97,8.97,8.91,9.09,9.11,9.11),
        c(9.04,9.08,9.05,9.25,9.28,9.27),
        c(9,9.01,9,9.2,9.23,9.2),
        c(8.99,8.99,8.98,9.18,9.2,9.19),
        c(8.93,8.97,8.97,9.18,9.2,9.18)
    )
    dim(z) <- c(15,6)
    z2 <- z + 1
    z3 <- z - 1
    
     plot_ly(showscale = FALSE) %>%
        add_surface(z = ~z) %>%
        add_surface(z = ~z2, opacity = 0.98) %>%
        add_surface(z = ~z3, opacity = 0.98)
    
    # Create a shareable link to your chart
    # Set up API credentials: https://plot.ly/r/getting-started
    chart_link = api_create(p, filename="surface-3")
    chart_link
    
    ##############
    # f
    # get weight-time-abundance data from wim
    
    multi <- function(x,y){x*y}
    
 #   species_sheet_plot <- function(sim,weight_exp=1,t_save=0.1,t_max=15){
    weight_exp<- 1
    
        my_times <- seq(from = 0, to = t_max, by = t_save)
        
        biomass_spectra_dynamics_gen <- sweep(sim@n,3,sim@params@w^weight_exp,multi)
        
        truncated_log_biomass_spectra_dynamics_gen <- log(biomass_spectra_dynamics_gen)
        truncated_log_biomass_spectra_dynamics_gen[biomass_spectra_dynamics_gen==0] <- min(log(biomass_spectra_dynamics_gen)[biomass_spectra_dynamics_gen>0])-1
        
        #x <- my_times
        #y <- log(p@w)
        z1 <- truncated_log_biomass_spectra_dynamics_gen[,1,]
        z2 <- truncated_log_biomass_spectra_dynamics_gen[,2,]
        z3 <- truncated_log_biomass_spectra_dynamics_gen[,3,]
        z4 <- truncated_log_biomass_spectra_dynamics_gen[,4,]
        z5 <- truncated_log_biomass_spectra_dynamics_gen[,5,]
        z6 <- truncated_log_biomass_spectra_dynamics_gen[,6,]
        z7 <- truncated_log_biomass_spectra_dynamics_gen[,7,]
        z8 <- truncated_log_biomass_spectra_dynamics_gen[,8,]
        z9 <- truncated_log_biomass_spectra_dynamics_gen[,9,]
        z10 <- truncated_log_biomass_spectra_dynamics_gen[,10,]
        z11 <- truncated_log_biomass_spectra_dynamics_gen[,11,]
        z12 <- truncated_log_biomass_spectra_dynamics_gen[,12,]
        z13 <- truncated_log_biomass_spectra_dynamics_gen[,13,]
        z14 <- truncated_log_biomass_spectra_dynamics_gen[,14,]
        
        plot_ly(showscale = FALSE) %>%
            add_surface(z = ~z1, x = ~my_times, y = ~log(p@w), opacity = 0.5) %>%
            add_surface(z = ~z2, x = ~my_times, y = ~log(p@w), opacity = 0.5)  %>%
            add_surface(z = ~z3, x = ~my_times, y = ~log(p@w), opacity = 0.5)  %>%
            add_surface(z = ~z4, x = ~my_times, y = ~log(p@w), opacity = 0.5)  %>%
            add_surface(z = ~z5, x = ~my_times, y = ~log(p@w), opacity = 0.5)  %>%
            add_surface(z = ~z6, x = ~my_times, y = ~log(p@w), opacity = 0.5)  %>%
            add_surface(z = ~z7, x = ~my_times, y = ~log(p@w), opacity = 0.98,color="red")  %>%
            add_surface(z = ~z8, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="blue")  %>%
            add_surface(z = ~z9, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="green")  %>%
            add_surface(z = ~z10, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="orange")  %>%
            add_surface(z = ~z11, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="yellow")  %>%
            add_surface(z = ~z12, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="purple")  %>%
            add_surface(z = ~z13, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="pink")  %>%
            add_surface(z = ~z14, x = ~my_times, y = ~log(p@w), opacity = 0.98) 
            
        
    #     ####################
    #     
    #     plot_ly(showscale = FALSE) %>%
    #         add_surface(z = ~z) %>%
    #         add_surface(z = ~z2, opacity = 0.98) %>%
    #         add_surface(z = ~z3, opacity = 0.98)
    #     
    #         persp3D(z=truncated_log_biomass_spectra_dynamics,y=log(p@w),x=my_times,xlab="time", ylab="log(weight)", zlab = "log(biomass_density",zlim=c(min(truncated_log_biomass_spectra_dynamics),-0.1),phi = 40, theta = 40)
    #     
    #     
    #     
    #     #################
    #     
    #     
    #     z3 <- z - 1
    #     
    #     plot_ly(showscale = FALSE) %>%
    #         add_surface(z = ~z) %>%
    #         add_surface(z = ~z2, opacity = 0.98) %>%
    #         add_surface(z = ~z3, opacity = 0.98)
    #     
    # }

        ########################### multiply by w^2 instead to get flatter plots #####
        
        #   species_sheet_plot <- function(sim,weight_exp=1,t_save=0.1,t_max=15){
        weight_exp<- 2
        
        my_times <- seq(from = 0, to = t_max, by = t_save)
        
        biomass_spectra_dynamics_gen <- sweep(sim@n,3,sim@params@w^weight_exp,multi)
        
        truncated_log_biomass_spectra_dynamics_gen <- log(biomass_spectra_dynamics_gen)
        truncated_log_biomass_spectra_dynamics_gen[biomass_spectra_dynamics_gen==0] <- min(log(biomass_spectra_dynamics_gen)[biomass_spectra_dynamics_gen>0])-1
        
        #x <- my_times
        #y <- log(p@w)
        z1 <- truncated_log_biomass_spectra_dynamics_gen[,1,]
        z2 <- truncated_log_biomass_spectra_dynamics_gen[,2,]
        z3 <- truncated_log_biomass_spectra_dynamics_gen[,3,]
        z4 <- truncated_log_biomass_spectra_dynamics_gen[,4,]
        z5 <- truncated_log_biomass_spectra_dynamics_gen[,5,]
        z6 <- truncated_log_biomass_spectra_dynamics_gen[,6,]
        z7 <- truncated_log_biomass_spectra_dynamics_gen[,7,]
        z8 <- truncated_log_biomass_spectra_dynamics_gen[,8,]
        z9 <- truncated_log_biomass_spectra_dynamics_gen[,9,]
        z10 <- truncated_log_biomass_spectra_dynamics_gen[,10,]
        z11 <- truncated_log_biomass_spectra_dynamics_gen[,11,]
        z12 <- truncated_log_biomass_spectra_dynamics_gen[,12,]
        z13 <- truncated_log_biomass_spectra_dynamics_gen[,13,]
        z14 <- truncated_log_biomass_spectra_dynamics_gen[,14,]
        
        plot_ly(showscale = FALSE) %>%
            add_surface(z = ~z1, x = ~my_times, y = ~log(p@w), opacity = 0.5) %>%
            add_surface(z = ~z2, x = ~my_times, y = ~log(p@w), opacity = 0.5)  %>%
            add_surface(z = ~z3, x = ~my_times, y = ~log(p@w), opacity = 0.5)  %>%
            add_surface(z = ~z4, x = ~my_times, y = ~log(p@w), opacity = 0.5)  %>%
            add_surface(z = ~z5, x = ~my_times, y = ~log(p@w), opacity = 0.5)  %>%
            add_surface(z = ~z6, x = ~my_times, y = ~log(p@w), opacity = 0.5)  %>%
            add_surface(z = ~z7, x = ~my_times, y = ~log(p@w), opacity = 0.98,color="red")  %>%
            add_surface(z = ~z8, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="blue")  %>%
            add_surface(z = ~z9, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="green")  %>%
            add_surface(z = ~z10, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="orange")  %>%
            add_surface(z = ~z11, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="yellow")  %>%
            add_surface(z = ~z12, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="purple")  %>%
            add_surface(z = ~z13, x = ~my_times, y = ~log(p@w), opacity = 0.98, color="pink")  %>%
            add_surface(z = ~z14, x = ~my_times, y = ~log(p@w), opacity = 0.98) 
        
        