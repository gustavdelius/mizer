#' ---
#' title: "Comparing analytic and numeric results"
#' output: html_document
#' ---
library(mizer)

#' We want to investigate how well the mizer numerics can reproduce the
#' analytic solution that we can obtain in the case of non-interacting
#' species. 
#' 

#' ## Without density dependence
# ----
#' ### Set parameters 
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

w_min <- 1e-2
w_inf <- 1e3
w_mat <- 10^2.4  # About 251, a quarter of w_inf
min_w_pp <- 1e-7  # Only have to make sure the smallest fish are perfectly fed
# Chose number of gridpoints so that w_mat and w_inf lie on gridpoints
no_w <- log10(w_inf/w_min)*10+1  

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

# ----
#' ### Comparison of solutions
#' 
#' If all is well, it should stay close to the steady-state solution.
plot(w, sim@n[dim(sim@n)[1],1,]*w, log="x", type="l", ylab="Biomass density")
lines(w, n_exact[1,]*w, col="blue")
#' There is a sizeable relative error at small sizes. Also the numerics are
#' off substantially at the maximum size.
#' We plot the relative error but exclude that maximum size.
relative_error <- (n_exact[1,]-sim@n[dim(sim@n)[1],1,])/n_exact[1,]
plot(w[1:(no_w-1)], relative_error[1:(no_w-1)], type="l", log="x",
     xlab="w", ylab="Relative error")

#' Let's us increase the number of steps by a factor of 10
no_w <- log10(w_inf/w_min)*100+1
params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
                      kappa = kappa, min_w = w_min, max_w = w_inf, no_w = no_w, 
                      min_w_pp = min_w_pp, w_pp_cutoff = w_inf, r_pp = r_pp)
# New grid, so we need to recalculate various things
w <- params@w
params@mu_b[1, ] <- mu0 * w^(n-1)
params@srr <- function(rdi, species_params) {return(species_params$R)}
params@species_params$R <- R
params@psi[1, ] <- (w/w_inf)^(1-n)
params@psi[1, w < w_mat] <- 0
params@psi[1, w > w_inf] <- 1
params@interaction[] <- 0
n_mult <- (1 - (w/w_inf)^(1-n))^(pow-1) * (1 - (w_mat/w_inf)^(1-n))^(-pow)
n_mult[w < w_mat] <- 1
n_exact <- params@psi
n_exact[] <- R * (w_min/w)^(mu0/hbar) / (hbar * w^n) * n_mult
sim <- project(params, t_max=50, effort = 0, initial_n = n_exact)
plot(w, sim@n[dim(sim@n)[1],1,]*w, log="x", type="l", ylab="Biomass density")
lines(w, n_exact[1,]*w, col="blue")
#' Now the agreement is much better
relative_error <- (n_exact[1,]-sim@n[dim(sim@n)[1],1,])/n_exact[1,]
plot(w[1:(no_w-1)], relative_error[1:(no_w-1)], type="l", log="x",
     xlab="w", ylab="Relative error")
#' We see that the error has also gone down by a factor of 10. So the numerical
#' method appears to be first-order in stepsize. It will be worthwhile to think
#' about how to improve this.

# ----
#' ## With density dependence
#' 
#' Now we repeat the same calculation with the same parameters, except that
#' we turn on density dependence
chi <- 0.05
ddd <- (kappa * w_mat^(-lambda))^chi
params@chi <- chi
params@ddd <- ddd
params@mu_b[1, ] <- mu0 * w^(n-1) * ddd
#' Now there is no elementary analytic solution for the abundance of mature 
#' fish (the solution involves a hypergeometric function), but for juvenile
#' fish we have
#' $$ 
#' N(w) = \frac{1}{\bar{h}w^n} \left(
#'      R^{-\chi} - \frac{\mu_0}{n \bar{h}^{\chi+1}} (w^{-n\chi}-w_0^{-n\chi})
#'    \right)^{-1/\chi}
#' $$
n_exact[] <- 1/(hbar * w^n) * (
        R^(-chi) - mu0/(n*hbar^(chi+1)) * (w^(-n*chi)-w_min^(-n*chi))
      )^(-1/chi)
#' We run the simulation and compare the results (in black) to the above exact
#' solution (in blue)
sim <- project(params, t_max=5, effort = 0, initial_n = n_exact)
plot(w, sim@n[dim(sim@n)[1],1,]*w, log="xy", type="l", ylab="Biomass density")
lines(w, n_exact[1,]*w, col="blue")
#' Note that we used a logarithmic scale also on the y axis because now the
#' biomass is no longer approximately constant for juveniles but changes
#' over several orders of magnitude.
#' 
#' It is clear from the shape of the mature spectrum that the system has not
#' yet had time to reach steady-state. So we run it for longer.
sim <- project(params, t_max=50, effort = 0, initial_n = n_exact)
plot(w, sim@n[dim(sim@n)[1],1,]*w, log="xy", type="l", ylab="Biomass density")
lines(w, n_exact[1,]*w, col="blue")
#' We plot the relative error in the juvenile region
relative_error <- (n_exact[1,]-sim@n[dim(sim@n)[1],1,])/n_exact[1,]
plot(w[w<w_mat], relative_error[w<w_mat], type="l", log="x",
     xlab="w", ylab="Relative error")

#' The number of eggs produced is close to the number of eggs R used.
rdi <- getRDI(params, matrix(sim@n[dim(sim@n)[1], , ], nrow = 1), params@cc_pp)[1,1]
rdi/R
#' Let us see what happens when we take away the restriction to the number
#' of eggs and instead use all eggs.
params@srr <- function(rdi, species_params) {return(rdi)}
sim <- project(params, t_max=200, effort = 0, initial_n = n_exact)
#' We use the biomass plot to judge whether steady-state has been reached
plotBiomass(sim, print_it=FALSE)
#' This has settled down nicely.
#' 
#' We plot the solution (in black) together with the analytic juvenile solution
#' appropriate for the new rate of egg production
rdi <- getRDI(params, matrix(sim@n[dim(sim@n)[1], , ], nrow = 1), params@cc_pp)[1,1]
n_exact[] <- 1/(hbar * w^n) * (
  rdi^(-chi) - mu0/(n*hbar^(chi+1)) * (w^(-n*chi)-w_min^(-n*chi))
    )^(-1/chi)
plot(w, sim@n[dim(sim@n)[1],1,]*w, log="xy", type="l", ylab="Biomass density")
lines(w, n_exact[1,]*w, col="blue")

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
  params@mu_b[i, ] <- mu0 * w^(n-1) * ddd
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
lines(w, colSums(sim@n[t,,]), col="red")
lines(w, kappa*w^(-lambda), col="orange")
lines(w, 7e10*w^(-lambda), col="green")
#' We overlay on that plot the community spectrum (in red), the background
#' spectrum (in orange). This shows that the community is too high compared
#' to the background. 


