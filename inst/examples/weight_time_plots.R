############## humboldt current setup
######### another copy of the preamble from humboldt.R

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
p <- steady(p, effort = effort, t_max = 500,  tol = 1e-2)

############## alter linecolors and line types
p@species_params$linetype[is.na(p@species_params$linetype)] <- "solid"
p@species_params$linecolour <- sapply((1:length(p@species_params$w_mat)),
                                      function(x) rgb(x/length(p@species_params$w_mat),.74,.18))



############################ sim from perturbed initial condition


# nn <- p@initial_n
# 
# for (i in 1:no_sp){
#     nn[i, ] <- runif(length(nn[i, ]))
# }
nn <- p@initial_n
nn[1, ] <- 2*nn[1,]

t_max <- 15
t_save <- 0.1
sim <- project(p, t_max = t_max, t_save = t_save, effort = effort, initial_n = nn)
plot(sim)

#! ?? line colors have been altered, so why don't they appear in plot(sim) ?

##############################################

########## weight-time plot showing most abundant species

MMa <- matrix(0,nrow = dim(sim@n)[1], ncol = dim(sim@n)[3])

for (t in (1:dim(sim@n)[1])){
    for (w in 1:dim(sim@n)[3]){
        MMa[t,w] <- which.max(sim@n[t,,w])
    }
    
}

################## working finished code for ggplot2 and making weight-time plot

library(reshape2)
library(ggplot2)


t_max <- 15
t_save <- 0.1
times <- seq(0,t_max,length.out = dim(sim@n)[1])
weights <- sim@params@w
m <- MMa
df <- melt(m)
df$value <- as.factor(df$value)
colnames(df) <- c("time",  "log10_weight",  "value")
df$time <- times[df$time]
df$log10_weight <- log10(weights[df$log10_weight])
ggplot(df) + geom_raster(aes(fill=value,x=time,y=log10_weight))+
    scale_fill_manual(values = sim@params@species_params$linecolour)
#! # again, why are these different ?
ggplot(df) + geom_raster(aes(fill=value,x=time,y=log10_weight))+
    scale_fill_manual(values = sapply((1:max(as.numeric(df$value))),function(x) rgb(x/max(as.numeric(df$value)),.74,.18))
    )

################## get dominant predator values at different weight-time 
# using species 2, but actually prey identity j is irrelevant since theta_ij=1

j <- 2

# time by weight                  
MM <- matrix(0,nrow = dim(sim@n)[1], ncol = dim(sim@n)[3])

for (t in (1:dim(sim@n)[1])){
    n <- sim@n[t,,]
    mumu <- sim@params@psi
    for (i in 1:no_sp){
        m <- n
        m[,] <- 0
        m[i,] <- n[i,]
        full_pr <- getPredRate(sim@params,m,sim@n_pp[t,])
        mumu[i,] <- full_pr[j,(dim(full_pr)[2]-length(sim@params@w)+1):dim(full_pr)[2]]
    }
    for (w in 1:dim(sim@n)[3]){
        MM[t,w] <- which.max(mumu[,w])
    }
    
}

################## use ggplot2 to plot dominant predator

library(reshape2)
library(ggplot2)


t_max <- 15
t_save <- 0.1
times <- seq(0,t_max,length.out = dim(sim@n)[1])
weights <- sim@params@w
m <- MM
df <- melt(m)
df$value <- as.factor(df$value)
colnames(df) <- c("time",  "log10_weight",  "value")
df$time <- times[df$time]
df$log10_weight <- log10(weights[df$log10_weight])

#! why does this plot look all green ?
ggplot(df) + geom_raster(aes(fill=value,x=time,y=log10_weight))+
    scale_fill_manual(values = sim@params@species_params$linecolour)


ggplot(df) + geom_raster(aes(fill=value,x=time,y=log10_weight))+
    scale_fill_manual(values = sapply((1:max(as.numeric(df$value))),function(x) rgb(x/max(as.numeric(df$value)),.74,.18))
    )



