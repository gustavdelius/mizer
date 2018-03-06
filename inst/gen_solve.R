#' ---
#' title: "A Stable Community"
#' author: "Richard Southwell"
#' output: pdf_document
#' ---

library(progress)
library(mizer)

#' We use the analytic solution to our trait based model, with variable egg size, to set up a 
#' community thats aggregation produces a background power law size spectrum that is self stableizing.
#' The stablity of this fixed point is promoted by our inclusion of density dependence.
#' 

#' ## Without density dependence
# ----
#' ### Set parameters 
#' Global Parameters
f0 <- 0.6
alpha <- 0.4
r_pp <- 1e-1
n <- 2/3
p <- n
q <- 3/4
lambda <- 2+q-n
kappa <- 7e10
#' Species Parameters
erepro <- 0.1 # Will be changed later
beta <- 100
sigma <- 1.3
h <- 30
ks <- 4

# ----
#' ### Set grid points and characteristic sizes 


dist_sp <- 0.2
no_sp <- 2/dist_sp + 1
species <- 1:no_sp
x_min <- seq(-4, by = dist_sp, length.out = no_sp)
w_min <- 10^x_min
w_inf <- 10^(x_min+5)
w_mat <- 10^(x_min+4.4)  # This is about a quarter of w_inf
min_w <- min(w_min)
max_w <- max(w_inf)
no_w <- log10(max_w/min_w)*100+1
min_w_pp <- 1e-8

# ----
#' ### Build Params Object 

species_params <- data.frame(
  species = 1:no_sp,
  w_min = w_min,
  w_inf = w_inf,
  w_mat = w_mat,
  h = h,
  ks = ks,
  beta = beta,
  sigma = sigma,
  z0 = 0,
  alpha = alpha,
  erepro = erepro,
  sel_func = "knife_edge", # not used but required
  knife_edge_size = 1000
)

params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
                      kappa = kappa, min_w = min_w, max_w = max_w, no_w = no_w, 
                      min_w_pp = min_w_pp, w_pp_cutoff = max_w, r_pp = r_pp,
                      chi = 0)
#' gamma is determined by mizerparams. Note that density dependence is currently off
gamma <- params@species_params$gamma[1]
w <- params@w
##################### solution maker : we have to specify effort and abundance_multipliers


# make a pretend solution, and extract growth rate
n_pretend <- params@psi
n_pretend[,] <- 0
n_pp_pretend <- params@kappa*params@w_full^(-params@lambda)
gg <- getEGrowth(params,n_pretend,n_pp_pretend)
# at this point we can compare and check that gg=(1-psi)*hbar*w^(n)

# make another pretend solution, and extract death rate
effort <- 0
n_pretend2 <- n_pretend
n_pretend2[1,] <- params@kappa*params@w^(-params@lambda)
n_pp_pretend2 <- n_pp_pretend
n_pp_pretend2[] <- 0
mumu <- getZ(params, n_pretend2, n_pp_pretend2, effort)
# we model by making the background abundances etc,. to realize these growth and death rates
# we could have fishing effort etc. on if we like, or we could fill to a mu0*w^(n-1) type curve

# now we construct the solutions
no_sp <- dim(params@psi)[1]
abundance_multipliers <- rep(10^(-8),no_sp)
w_inf_idx <- rep(0,no_sp)
n_real <- n_pretend
for (i in (1:no_sp)){
  w_inf_idx[i] <- length(params@w[params@w<=params@species_params$w_inf[i]])
  integrand <- mumu[i,params@species_params$w_min_idx[i]:w_inf_idx[i]]/gg[i,params@species_params$w_min_idx[i]:w_inf_idx[i]]
  n_real[i,params@species_params$w_min_idx[i]:w_inf_idx[i]] <- abundance_multipliers[i]*
    exp(-cumsum(integrand))/gg[i,params@species_params$w_min_idx[i]:w_inf_idx[i]]
}

# plot and check above matches with stable_community.R step by step

##################

#' ### Setup plankton

plankton_vec <- (kappa*w^(-lambda))-colSums(n_real)
negwarn <- function() warning("negative entries")
if (length(plankton_vec[plankton_vec<0])>0){
  negwarn()
}
params@cc_pp[sum(params@w_full<=params@w[1]):length(params@cc_pp)] <- plankton_vec
initial_n_pp <- params@cc_pp
# The cc_pp factor needs to be higher than the desired steady state in
# order to compensate for predation mortality
m2_background <- getM2Background(params, n_real, initial_n_pp)
params@cc_pp <- (1+m2_background/params@rr_pp) * initial_n_pp

################# setup background death

mu_impartial <- getZ(params, n_real, initial_n_pp)

for (i in 1:no_sp) {
  params@mu_b[i, ] <- mumu[i,] - mu_impartial[i,]
}

### mmet reproduction boundary condition

rdi <- getRDI(params, n_real, initial_n_pp)
gg_new <- getEGrowth(params, n_real, initial_n_pp)
# # #
effort <- 0
mumu <- getZ(params, n_output, initial_n_pp, effort = effort)
erepro_final <- rdi
for (i in (1:no_sp)){
  #  erepro_final[i] <- erepro*(gg[i,params@species_params$w_min_idx[i]]*n_output[i,params@species_params$w_min_idx[i]])/
  #    rdi[i]
  gg0 <- gg[i,params@species_params$w_min_idx[i]]
  mumu0 <- mumu[i,params@species_params$w_min_idx[i]]
  DW <- params@dw[params@species_params$w_min_idx[i]]
  erepro_final[i] <- erepro*(n_output[i,params@species_params$w_min_idx[i]]*(gg0+DW*mumu0))/rdi[i]
}

params@species_params$erepro <- erepro_final

