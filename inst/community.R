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
    knife_edge_size = 1000,
    R = R
)

params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
                      kappa = kappa, min_w = w_min, max_w = w_inf, no_w = no_w, 
                      min_w_pp = min_w_pp, w_pp_cutoff = w_inf, r_pp = r_pp,
                      chi = chi)

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
f0 <- 0.6027546
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
#' We run the simulation and compare the results (in black) to the above exact
#' solution (in blue)
params@srr <- function(rdi, species_params) {return(species_params$R)}
params@psi[1, ] <- (w/w_inf)^(1-n)
params@psi[1, w < w_mat] <- 0
params@psi[1, w > w_inf] <- 1
#params@srr <- function(rdi, species_params) {return(rdi)}
sim <- project(params, t_max=50, effort = 0, initial_n = n_exact)
plot(w, sim@n[dim(sim@n)[1],1,]*w, log="xy", type="l", ylab="Biomass density")
lines(w, n_exact[1,]*w, col="blue")
plot(sim)
getFeedingLevel(sim)
