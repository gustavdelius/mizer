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

params@species_params[no_traits+1,] <- NS_params@species_params[6,]

params@species_params[no_traits+1,]

NS_params@species_params[6,]

#for (i in (1:dim(params@species_params)[2])){
for (i in (2:10)){
  
  params@species_params[no_traits+1,i] <- NS_params@species_params[6,i]
  
}
i <- 12
params@species_params[no_traits+1,i] <- NS_params@species_params[6,i]
i <- 14
params@species_params[no_traits+1,i] <- NS_params@species_params[6,i]
i <- 15
params@species_params[no_traits+1,i] <- NS_params@species_params[6,i]
i <- 16
params@species_params[no_traits+1,i] <- NS_params@species_params[6,i]
i <- 17
params@species_params[no_traits+1,i] <- NS_params@species_params[6,i]
i <- 18
params@species_params[no_traits+1,i] <- NS_params@species_params[6,i]

sim <- project(params, t_max=76, effort = 0)
plot(sim)

# add mackrel