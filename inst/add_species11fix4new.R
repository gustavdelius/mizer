#' ---
#' title: "Gurnard in static background"
#' author: "Richard Southwell"
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
#f0 <- 0.6
#alpha <- 0.4
r_pp <- 10^18  # Choosing a high value because we want the plankton to stay
# at its power-law steady state
n <- 2/3
p <- n
q <- 0.8
lambda <- 2+q-n
#erepro <- 0.1
R <- 1e10  # The rate of reproduction

#beta <- 100
#sigma <- 1.3
#h <- 30
#ks <- 4
kappa <- 7e10

#w_min <- 1e-2
#w_inf <- 1e3
#w_mat <- 10^2.4  # About 251, a quarter of w_inf
#min_w_pp <- 1e-7  # Only have to make sure the smallest fish are perfectly fed
# Chose number of gridpoints so that w_mat and w_inf lie on gridpoints
#no_w <- log10(w_inf/w_min)*100+1  

#######################################

# Gurnard Parameters
beta <- 283
sigma <- 1.8
k_vb <- 0.266
alpha <- 0.6
h <- 19.37727
gamma <- 2.948552e-11
#z0 <- 0.06863713 
z0 <- 0
ks <- 3.875454
erepro <- 0.1
#' used (3.17) to get f0 from gamma. Note this should be corrected when chi>0
f0 <- gamma/(gamma+h*exp(-((lambda-2)^2)*sigma*sigma/2)*(beta^(2-lambda))/(sqrt(2*pi)*kappa*sigma))

# # # # weight specifiers

no_w <- 2000
w_min <- 1e-3
w_inf <- 668
w_mat <- 39
min_w_pp <- exp(-5*sigma)*w_min/beta  # Only have to make sure the smallest fish are perfectly fed

## ## ## ## ## ##

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
  gamma = gamma
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
n_exact[length(n_exact)] <- n_exact[length(n_exact)-1]
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
plot(w, sim@n[dim(sim@n)[1],1,], log="xy", type="l", ylab="Abundance")
lines(w, n_exact[1,], col="blue")


#' we shall turn on prey-switching before letting R become dynamic
#' we do this to keep things stable, but it remains an interesting question 
#' what is the stability when $\chi=0$, which we should indeed study

# ----
#' ## With density dependence
#' 
#' Now we repeat the same calculation with the same parameters, except that
#' we turn on density dependence
chi <- 0.05
#ddd <- (kappa * w_mat^(-lambda))^chi
ddd <- array(1,dim=c(1,length(w)))
ddd[1,] <- (kappa * w_mat^(-lambda))^chi
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
sim <- project(params, t_max=500, effort = 0, initial_n = n_exact)
plot(w, sim@n[dim(sim@n)[1],1,], log="xy", type="l", ylab="Abundance")
lines(w, n_exact[1,], col="blue")

#' The Gurnard parameters I'm using seem to make w_inf size fish build up 
#' (i.e., the abundance curves up at RHS), this is partially because I turned
#'  the background death off. But it makes me wonder if an analytic expression
#'  for the solution in the presence of additional weight-independent death is 
#'  possible. Then we could investigate how much such death is required to stop 
#'  N(W_inf)->inf 
#' It might also be useful for understanding the effects of fishing

#' Before we allow dynamic reproductive influx, we want to make sure the system 
#' naturally produces the R we want. Several paramaters can be used to do this
#' indeed, looking at our $\chi=0$ solution, we see the solution shape just depends 
#' on $n$ and $\mu_0/\hbar$ and $\hbar,$ (and $w_0, w_{mat}, w_{inf}$) and can expand our expressions for these 
#' quantities in terms of the other parameters, and think about the different ways to 
#' change these quantities to vary the RDI integral, to get R correct. However here we 
#' just use the basic way of alterting erepro  
rdi_old_erepro <- getRDI(params, matrix(sim@n[dim(sim@n)[1], , ], nrow = 1), params@cc_pp)[1,1]
params@species_params$erepro <- params@species_params$erepro * R/rdi_old_erepro
#' now let us check that we are now making the correct number of eggs
rdi <- getRDI(params, matrix(sim@n[dim(sim@n)[1], , ], nrow = 1), params@cc_pp)[1,1]

rdi/R
#' now we shall allow the number of eggs to be dynamic

params@srr <- function(rdi, species_params) {return(rdi)}
sim <- project(params, t_max=500, effort = 0, initial_n = n_exact)
plot(sim)
plot(w, sim@n[dim(sim@n)[1],1,], log="xy", type="l", ylab="Abundance")
lines(w, n_exact[1,], col="blue")

#' For these type of solutions that curl up, a large no_w seems to be required
#'  to avoid numerical errors when determining RDI

#' Things to do:
#' figure out how to get gurnard to exist when R is not fixed
#' turn on chi, predict effect
#' put gurnard in dynamic background
#' automate solution
#' Investigate the growth curves of these Gurnard
#' get analytic formulas for f0 and mu0 when chi>0
#' document this with a nice markdown doc, including pde solution details too
#' 

g <- getEGrowth(params, matrix(sim@n[dim(sim@n)[1], , ], nrow = 1), params@cc_pp)[1,]

#' calculate growth curve using ODE solving routine
g_fn <- approxfun(w, g)
myodefun <- function(t, state, parameters){
  return(list(g_fn(state)))
}
age <- (0:500)
library(deSolve)
weight <- ode(y = w[1], times = age, func = myodefun, parms = 1)[,2]
plot(age,weight)
#' calculate growth curve in another way via direct integration
#lines(cumsum((w[2:no_w]-w[1:(no_w-1)])/g[1:(no_w-1)]),w[1:(no_w-1)],xlab="age",ylab="weight",col="red")
f0
params@species_params$erepro

#' Next step
#' 
#' Program in exact soln as ddd slot
#' Repeat community construction
#'

##params@ddd==(kappa*w_mat^(-lambda))^chi