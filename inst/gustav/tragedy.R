# Source the mizer code without loading it as a package
source('./R/MizerParams-class.r')
source('./R/MizerSim-class.r')
source('./R/project_methods.r')
source('./R/selectivity_funcs.r')
source('./R/summary_methods.r')
source('./R/wrapper_functions.R')
source('./R/plots.r')
source('./R/project.r')

# Load required packages
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)

#' Selectivity function for investigating the tragedy of the commons
#'
#' One fisher uses a smaller mesh size than all the other fishers
#'
#' @param w The size of the individual.
#' @param knife_edge_size The size at which the knife-edge operates.
#' @param n_fishers Number of fishers
#' @param tragedy_size The size at which the exploiting fisher starts fishing
#' @export
tragedy_edge <- function(w, knife_edge_size, n_fishers, tragedy_size){
  sel <- rep(0, length(w))
  sel[w >= tragedy_size] <- 1/n_fishers
  sel[w >= knife_edge_size] <- 1
  return(sel)
} 

#' Modified version of set_trait_model which sets harvesting so that one fisher uses small net
set_tragedy_model <- function(no_sp = 10,
                            min_w_inf = 10,
                            max_w_inf = 1e5,
                            no_w = 100,
                            min_w = 0.001,
                            max_w = max_w_inf * 1.1,
                            min_w_pp = 1e-10,
                            no_w_pp = round(no_w)*0.3,
                            w_pp_cutoff = 1,
                            k0 = 50, # recruitment adjustment parameter
                            n = 2/3,
                            p = 0.75,
                            q = 0.9, 
                            eta = 0.25,
                            r_pp = 4,
                            kappa = 0.005,
                            lambda = 2+q-n,
                            alpha = 0.6,
                            ks = 4,
                            z0pre = 0.6,
                            h = 30,
                            beta = 100,
                            sigma = 1.3,
                            f0 = 0.5,
                            gamma = NA,
                            knife_edge_size = 1000,
                            gear_names = "knife_edge_gear",
                            tragedy_size = 100,
                            n_fishers = 5,
                            ...){
  # If not supplied, calculate gamma using equation 2.1 in A&P 2010
  if(is.na(gamma)){
    alpha_e <- sqrt(2*pi) * sigma * beta^(lambda-2) * exp((lambda-2)^2 * sigma^2 / 2) # see A&P 2009
    gamma <- h * f0 / (alpha_e * kappa * (1-f0)) # see A&P 2009 
  }
  w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
  w_mat <- w_inf * eta
  
  # Check gears
  if (length(knife_edge_size) > no_sp){
    stop("There cannot be more gears than species in the model")
  }
  if ((length(knife_edge_size) > 1) & (length(knife_edge_size) != no_sp)){
    warning("Number of gears is less than number of species so gear information is being recycled. Is this what you want?")
  }
  if ((length(gear_names) != 1) & (length(gear_names) != no_sp)){
    stop("Length of gear_names argument must equal the number of species.")
  }
  
  # Make the species parameters data.frame
  trait_params_df <- data.frame(
    species = 1:no_sp,
    w_inf = w_inf,
    w_mat = w_mat,
    h = h, # max food intake
    gamma = gamma, # vol. search rate,
    ks = ks,# standard metabolism coefficient,
    beta = beta,
    sigma = sigma,
    z0 = z0pre * w_inf^(n-1), # background mortality
    alpha = alpha,
    #r_max = r_max,
    sel_func = "tragedy_edge",
    knife_edge_size = knife_edge_size,
    tragedy_size = tragedy_size,
    n_fishers = n_fishers,
    gear = gear_names,
    erepro = 1 # not used but included out of necessity
  )
  # Make the MizerParams
  trait_params <- MizerParams(trait_params_df, min_w = min_w, max_w=max_w, no_w = no_w, min_w_pp = min_w_pp, w_pp_cutoff = w_pp_cutoff, n = n, p=p, q=q, r_pp=r_pp, kappa=kappa, lambda = lambda) 
  # Sort out maximum recruitment - see A&P 2009
  # Get max flux at recruitment boundary, R_max
  # R -> | -> g0 N0
  # R is egg flux, in numbers per time
  # Actual flux at recruitment boundary = RDD = NDD * g0 (where g0 is growth rate)
  # So in our BH SRR we need R_max comparable to RDI (to get RDD)
  # R_max = N0_max * g0 (g0 is the average growth rate of smallest size, i.e. at f0 = 0.5)
  # N0 given by Appendix A of A&P 2010 - see Ken's email 12/08/13
  # Taken from Ken's code 12/08/13 - equation in paper is wrong!
  alpha_p <- f0 * h * beta^(2 * n - q - 1) * exp((2 * n * (q - 1) - q^2 + 1) * sigma^2 / 2)
  alpha_rec <- alpha_p / (alpha * h * f0 - ks)
  # Calculating dw using Ken's code - see Ken's email 12/08/13
  tmpA <- w_inf[1]
  tmpB <- (log10(w_inf[length(w_inf)]) - log10(w_inf[1])) / (no_sp - 1) # Difference between logged w_infs, fine
  dw_winf <- tmpB * tmpA *10^(tmpB*((1:no_sp)-1)) # ?
  N0_max <- k0 * w_inf^(n*2-q-3+alpha_rec) * dw_winf  # Why * dw_winf, not / ? Ken confirms * in email
  # No need to include (1 - psi) in growth equation because allocation to reproduction at this size = 0, so 1 - psi = 1
  g0 <- (alpha * f0 * h * trait_params@w[1]^n - ks * trait_params@w[1]^p)
  r_max <- N0_max * g0
  
  trait_params@species_params$r_max <- r_max
  
  return(trait_params)
}


#' Plot the yield
ploty <- function(sim, sim_tragedy) {
  y <- data.frame(time=1:50,equal=rowSums(getYield(sim)),exploiting=rowSums(getYield(sim_tragedy)))
  ym <- melt(y[20:50, ], id.vars="time")
  qplot(time, value, data=ym, color=variable, geom="line")
}

params_knife <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = 1000)
sim <- project(params_knife, effort = 0.5, t_max = 50)
plot(sim)

# Simulation with exploiting fisher
params_tragedy<- set_tragedy_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = 1000,
                                  n_fishers = 5, tragedy_size = 100)
sim_tragedy <- project(params_tragedy, effort = 0.5, t_max = 50)
plot(sim_tragedy)

# Plot the yield
ploty(sim, sim_tragedy)

# Not reached steady state yet, so run for longer
sim <- project(params_knife, effort = 0.5, t_max = 50, initial_n = sim@n[51,,])
sim_tragedy <- project(params_tragedy, effort = 0.5, t_max = 50, initial_n = sim_tragedy@n[51,,])
ploty(sim, sim_tragedy)

# After running the above two sections until steady state, now change the
# mesh size
params_knife <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = 0.13)
params_tragedy<- set_tragedy_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = 0.13,
                                   n_fishers = 5, tragedy_size = 0.12)
sim <- project(params_knife, effort = 0.5, t_max = 50, initial_n = sim@n[51,,])
sim_tragedy <- project(params_tragedy, effort = 0.5, t_max = 50, initial_n = sim_tragedy@n[51,,])
ploty(sim, sim_tragedy)



params_knife <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5, knife_edge_size = 0.12)
sim <- project(params_knife, effort = 0.5, t_max = 50, initial_n = sim@n[51,,])
ys <- rowSums(getYield(sim))
ys[50]
