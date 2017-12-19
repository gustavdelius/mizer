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
r_pp <- 10^18  # Choosing a high value because we want the plankton to stay
# at its power-law steady state
n <- 2/3
p <- n
q <- 0.95
lambda <- 2+q-n
erepro <- 0.1
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
                      min_w_pp = min_w_pp, w_pp_cutoff = w_inf, r_pp = r_pp,
                      chi = chi)

gamma <- params@species_params$gamma[1]
#' Mizer determined gamma with chi=0, so the feeding level will not be exactly
#' as set, but it will be close enough.
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
params@mu_b[1, ] <- mu0 * w^(n-1) * params@ddd
#' The rate at which energy is available for growth and reproduction
#' is $\bar{h} w^n$, with $\bar{h} = \alpha h f_0 - k_s$.
hbar <- alpha * h * f0 - ks
#' For juvenile fish we have
#' $$ 
#' N(w) = \frac{1}{\bar{h}w^n} \left(
#'      R^{-\chi} - \frac{\mu_0}{n \bar{h}^{\chi+1}} (w^{-n\chi}-w_0^{-n\chi})
#'    \right)^{-1/\chi}
#' $$
n_exact <- params@psi  # Just to get array with correct dimensions and names
n_exact[] <- 1/(hbar * w^n) * (
    R^(-chi) - mu0/(n*hbar^(chi+1)) * (w^(-n*chi)-w_min^(-n*chi))
)^(-1/chi)

#' Make sure that the rate of reproduction is R
params@srr <- function(rdi, species_params) {return(species_params$R)}
params@species_params$R <- R
#' We use a step function for the maturity function
params@psi[1, ] <- (w/w_inf)^(1-n)
params@psi[1, w < w_mat] <- 0
params@psi[1, w > w_inf] <- 1
#' We switch off the self-interaction
params@interaction[] <- 0

#' We run the simulation and compare the results (in black) to the above exact
#' solution (in blue)
params@srr <- function(rdi, species_params) {return(species_params$R)}
params@psi[1, ] <- (w/w_inf)^(1-n)
params@psi[1, w < w_mat] <- 0
params@psi[1, w > w_inf] <- 1

#' Now let the number of eggs free
params@srr <- function(rdi, species_params) {return(rdi)}
sim <- project(params, t_max=200, effort = 0, initial_n = n_exact)

#' ## Putting several species together
#' 
#' We'll use our solution from above as a template
n_init <- sim@n[dim(sim@n)[1],1,]
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
    gamma = gamma,
    erepro = erepro,
    sel_func = "knife_edge", # not used but required
    knife_edge_size = 1000
)

#' First without interaction
params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
                      kappa = kappa, min_w = min_w, max_w = max_w, no_w = no_w, 
                      min_w_pp = min_w_pp, w_pp_cutoff = max_w, r_pp = r_pp,
                      chi = chi)

w <- params@w
params@srr <- function(rdi, species_params) {return(rdi)}
params@interaction[] <- 0
initial_n <- params@psi
for (i in 1:no_sp) {
    initial_n[i, params@species_params$w_min_idx[i]:(params@species_params$w_min_idx[i]+length(n_init)-1)] <-
        n_init * (w_min[no_sp]/w_min[i])^lambda
    params@psi[i, ] <- (w/w_inf[i])^(1-n)
    params@psi[i, w < w_mat[i]] <- 0
    params@psi[i, w > w_inf[i]] <- 1
    params@mu_b[i, ] <- mu0 * w^(n-1) * params@ddd[11]
}
sim <- project(params, t_max=5 ,effort = 0, initial_n = initial_n)
#' We plot the species biomass over time to see that we are in the steady state.
plotBiomass(sim, print_it = FALSE)

params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
                      kappa = kappa, min_w = min_w, max_w = max_w, no_w = no_w, 
                      min_w_pp = min_w_pp, w_pp_cutoff = min(w_mat), r_pp = r_pp,
                      chi = chi)

w <- params@w
params@srr <- function(rdi, species_params) {return(rdi)}
initial_n <- params@psi

for (i in 1:no_sp) {
    initial_n[i, params@species_params$w_min_idx[i]:(params@species_params$w_min_idx[i]+length(n_init)-1)] <-
        n_init * (w_min[no_sp]/w_min[i])^lambda
    params@psi[i, ] <- (w/w_inf[i])^(1-n)
    params@psi[i, w < w_mat[i]] <- 0
    params@psi[i, w > w_inf[i]] <- 1
}
m2 <- getM2(params, initial_n, params@cc_pp)
for (i in 1:no_sp) {
    params@mu_b[i, ] <- mu0 * w^(n-1) * params@ddd[11] - m2[i, ]
}

sim <- project(params, t_max=5 ,effort = 0, initial_n = initial_n)
#' We plot the species biomass over time to see that we are in the steady state.
plotBiomass(sim, print_it = FALSE)

#' We again compare the analytic solution (in black) with the numeric solution
#' (in blue) and see that the blue covers the black almost perfectly.
t <- dim(sim@n)[1]
plot(w, sim@n[t,1,], log="xy", type="l", ylim=c(10^(4),10^(15)))
lines(w, sim@n[1,1,], col="blue")
for (i in 2:no_sp){
    lines(w, sim@n[t,i,])
    lines(w, sim@n[1,i,], col="blue")
}
lines(w, sim@n_pp[1, 401:1101], col="green")
lines(w, colSums(sim@n[1,,])+sim@n_pp[1, 401:1101], col="red") 

plot(w, params@mu_b[1, ], log="xy", type="l")
lines(w, mu0 * w^(n-1) * params@ddd, col="blue")
