#' ---
#' title: "Assembling a community"
#' output: html_document
#' ---
library(mizer)

# ----
#' ## Set parameters
chi <- 0.05
f0 <- 0.6
alpha <- 0.4
r_pp <- 1  # Not used because we commented out the plankton dynamics
n <- 2/3
p <- n
q <- 0.75
lambda <- 2+q-n
erepro <- 0.1 # Will be overwritten later
R <- 1e10  # The rate of reproduction

beta <- 100
sigma <- 1.3
h <- 30
ks <- 4
kappa <- 7e10

# ----
#' ## Calculate single-species solution
#' 
w_min <- 1e-2
w_inf <- 1e3
w_mat <- 10^2.4  # About 251, a quarter of w_inf
min_w_pp <- 1e-7  # Only have to make sure the smallest fish are perfectly fed
# Chose number of gridpoints so that w_mat and w_inf lie on gridpoints
no_w <- log10(w_inf/w_min)*100+1  

species_params <- data.frame(
    species = "Single",
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
                      kappa = kappa, min_w = w_min, max_w = w_inf, no_w = no_w, 
                      min_w_pp = min_w_pp, w_pp_cutoff = w_inf, r_pp = r_pp)

gamma <- params@species_params$gamma[1]
w <- params@w

# ----
#' ### Exact solution
#' 
#' We set the background death to $\mu_0 w^{n-1}$ and choose
#' $$ 
#' \mu_0 = (1-f_0) \sqrt{2\pi} \kappa \gamma \sigma 
#'            \beta^{n-1} \exp(\sigma^2 (n-1)^2 / 2)
#' $$
#' which is the death rate that is produced by predation if the predators
#' follow the same power law as the plankton. 
#' We could equally well have chosen any other constant for $\mu_0$.
mu0 <- (1-f0) * sqrt(2*pi) * kappa * gamma * sigma *
    (beta^(n-1)) * exp(sigma^2 * (n-1)^2 / 2)
params@mu_b[1, ] <- mu0 * w^(n-1)
#' The rate at which energy is available for growth and reproduction
#' is $\bar{h} w^n$, with $\bar{h} = \alpha h f_0 - k_s$.
hbar <- alpha * h * f0 - ks
#' n_exact is calculated using the analytic expression for the solution,
#' $$
#' N(w) = \frac{R}{\bar{h} w^n} (w_0/w)^{\mu_0/\bar{h}}
#'  \times \begin{cases} 1 & w<w_{mat}\\
#'  \left(1-(w/W_{inf})^{1-n}\right)^{\mu_0/\bar{h}/(1-n)-1}
#'  \left(1-(w_{mat}/W_{inf})^{1-n}\right)^{-\mu_0/\bar{h}/(1-n)} & w>w_{mat}
#'  \end{cases}
#' $$
pow <- mu0/hbar/(1-n)
n_mult <- (1 - (w/w_inf)^(1-n))^(pow-1) * (1 - (w_mat/w_inf)^(1-n))^(-pow)
n_mult[w < w_mat] <- 1
n_exact <- params@psi  # Just to get array with correct dimensions and names
n_exact[] <- R * (w_min/w)^(mu0/hbar) / (hbar * w^n) * n_mult

# ----
#' ### Mizer setup for non-interacting scenario
#'
#' Make sure that the rate of reproduction is R
params@srr <- function(rdi, species_params) {return(species_params$R)}
params@species_params$R <- R
#' We use a step function for the maturity function
params@psi[1, ] <- (w/w_inf)^(1-n)
params@psi[1, w < w_mat] <- 0
params@psi[1, w > w_inf] <- 1
#' We switch off the self-interaction
params@interaction[] <- 0

#' We start the simulation with the exact steady-state solution
sim <- project(params, t_max=50, effort = 0, initial_n = n_exact)
#' We use the biomass plot to judge whether steady-state has been reached
plotBiomass(sim, print_it=FALSE)
#' If all is well, it should stay close to the steady-state solution.
plot(w, sim@n[dim(sim@n)[1],1,]*w, log="xy", type="l", ylab="Biomass density")
lines(w, n_exact[1,]*w, col="blue")

#' We'll use our solution from above as a template
n_init <- sim@n[dim(sim@n)[1], 1, ]
n_init_mat <- matrix(n_init, nrow=1)

# set erepro to meet boundary condition
rdi <- getRDI(params, n_init_mat, params@cc_pp)[1]
gg <- getEGrowth(params, n_init_mat, params@cc_pp)[1]
mumu <- getZ(params, n_init_mat, params@cc_pp, effort = 0)[1]
erepro_final <- erepro * n_init[1] * (gg + params@dw[1]*mumu) / rdi
erepro_final

#' Now let the number of eggs free
params@species_params$erepro <- erepro_final
params@srr <- function(rdi, species_params) {return(rdi)}
#' turn on density dependence
params@chi <- 0.05
params@ddd <- n_init_mat^(params@chi)
#' and run the simulation
sim <- project(params, t_max=100, effort = 0, initial_n = n_init)
plotBiomass(sim, print_it=FALSE)
#' If all is well, it should stay close to the steady-state solution.
plot(w, sim@n[dim(sim@n)[1],1,]*w, log="xy", type="l", ylab="Biomass density")
lines(w, n_init*w, col="blue")

# ----
#' ## Putting several species together
#' 
#' We set up 11 species with egg size between $10^-4$ and $10^-2$, equally
#' spaced in log space
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
    erepro = erepro_final,
    sel_func = "knife_edge", # not used but required
    knife_edge_size = 1000
)

params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
                      kappa = kappa, min_w = min_w, max_w = max_w, no_w = no_w, 
                      min_w_pp = min_w_pp, w_pp_cutoff = min(w_mat), r_pp = r_pp,
                      chi = chi)

w <- params@w
params@srr <- function(rdi, species_params) {return(rdi)}
initial_n <- params@psi
initial_n[] <- 0
for (i in 1:no_sp) {
    initial_n[i, params@species_params$w_min_idx[i]:(params@species_params$w_min_idx[i]+length(n_init)-1)] <-
        n_init * (w_min[1]/w_min[i])^lambda
    params@psi[i, ] <- (w/w_inf[i])^(1-n)
    # We have to be very careful when comparing weights
    params@psi[i, w < (w_mat[i]-1e-10)] <- 0
    params@psi[i, w > (w_inf[i]-1e-10)] <- 1
}

#' Adjust abundance to get community spectrum to line up with background spectrum
v <- sqrt(min(w_mat)*max(w_mat))
v_idx <- length(w[w<v])
initial_n <- initial_n*(kappa*w[v_idx]^(-lambda))/sum(initial_n[,v_idx])

#' Use this rescaled solution as denominator in density dependence
nn <- initial_n
nn[nn==0] <- 1
params@ddd <- nn^(params@chi)

#' set background death so that it combines with the predation death to
#' the perfect power-law death.
#' This could lead to negative death
m2 <- getM2(params, initial_n, params@cc_pp)
for (i in 1:no_sp) {
    params@mu_b[i, ] <- mu0 * w^(n-1) - m2[i, ]
}

#' set the plankton background so that it combines with the community
#' to give the perfect power law
plankton_vec <- (kappa*w^(-lambda))-colSums(initial_n)
plankton_vec[plankton_vec<0] <- 0
plankton_vec[min(which(plankton_vec==0)):length(plankton_vec)] <- 0
params@cc_pp[params@w_full>=w[1]] <- plankton_vec

sim <- project(params, t_max=50 ,effort = 0, initial_n = initial_n)
#' We plot the species biomass over time to see that we are in the steady state.
plotBiomass(sim, print_it = FALSE)

#' We compare the solution at the final time (in black) with the 
#' initial solution (in blue).
t <- dim(sim@n)[1]
plot(w, sim@n[t,1,], log="xy", type="l", ylim=c(10^(4),10^(15)))
lines(w, sim@n[1,1,], col="blue")
for (i in 2:no_sp){
    lines(w, sim@n[t,i,])
    lines(w, sim@n[1,i,], col="blue")
}
lines(w, sim@n_pp[1, 401:1101], col="green")
lines(w, colSums(sim@n[1,,])+sim@n_pp[1, 401:1101], col="red") 

