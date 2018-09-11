############## use scale invariant setup
no_sp <- 2
p <- set_scaling_model(no_sp = no_sp)

############## alter linecolors and line types
p@species_params$linetype <- "solid"
p@species_params$linecolour <- sapply((1:length(p@species_params$w_mat)),
                                      function(x) rgb(x/length(p@species_params$w_mat),.74,.18))
###### find steady state

effort <- 0
t_max <- 30
t_save <- 0.1

sim <- project(p, t_max = t_max, t_save = t_save, effort = effort)
plot(sim)


p@initial_n <- sim@n[dim(sim@n)[1],,]
p@initial_n_pp <- sim@n_pp[dim(sim@n_pp)[1],]



############################ sim from perturbed initial condition


# nn <- p@initial_n
# 
# for (i in 1:no_sp){
#     nn[i, ] <- runif(length(nn[i, ]))
# }
nn <- p@initial_n
nn[1, ] <- nn[1,]*1/2

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



