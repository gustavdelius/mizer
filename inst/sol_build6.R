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
r_pp <- 1  # Choosing a high value because we want the plankton to stay
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
#ks <- 3
kappa <- 7e10


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
# gamma was set up in a wierd way in analytic_numerical_comparison, and calulated via Mizerparams
gamma <- 4.999819e-11
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
                      chi = 0)
#gamma <- params@species_params$gamma[1]
w <- params@w


##########


mu0 <- (1-f0) * sqrt(2*pi) * kappa * gamma * sigma *
  (beta^(n-1)) * exp(sigma^2 * (n-1)^2 / 2)
hbar <- alpha * h * f0 - ks
pow <- mu0/hbar/(1-n)
n_mult <- (1 - (w/w_inf[1])^(1-n))^(pow-1) * (1 - (w_mat[1]/w_inf[1])^(1-n))^(-pow)
n_mult[w < w_mat[1]] <- 1
n_mult[w >= w_inf[1]] <- 0
n_exact <- w  # Just to get array with correct dimensions and names

#use R=1 initially # n_exact[] <- R * (w_min/w)^(mu0/hbar) / (hbar * w^n) * n_mult
n_exact <- ((w_min[1]/w)^(mu0/hbar) / (hbar * w^n)) * n_mult
n_exact[w < w_min[1]] <- 0

initial_n <- params@psi
initial_n[,] <- 0
w_inf_idx <- w_inf
for (i in 1:no_sp) {
  #initial_n[i, params@species_params$w_min_idx[i]:(params@species_params$w_min_idx[i]+length(n_init)-1)] <-
  #  n_init * (w_min[no_sp]/w_min[i])^lambda
  ##initial_n[i, params@species_params$w_min_idx[i]:(params@species_params$w_min_idx[i]+length(n_exact)-1)] <-
  ##  n_exact * (w_min[no_sp]/w_min[i])^lambda
  # 500 = w_inf_idx-w_min_idx
  
  w_inf_idx[i] <- length(w[w<=w_inf[i]])
  initial_n[i, params@species_params$w_min_idx[i]:(params@species_params$w_min_idx[i]+(w_inf_idx[1]-params@species_params$w_min_idx[1]))] <-
    n_exact[params@species_params$w_min_idx[1]:(params@species_params$w_min_idx[1]+(w_inf_idx[1]-params@species_params$w_min_idx[1]))] * (w_min[1]/w_min[i])^lambda
  #params@psi[i, ] <- (w/w_inf[i])^(1-n)
  #params@psi[i, w < w_mat[i]] <- 0
  #params@psi[i, w > w_inf[i]] <- 1
  #params@mu_b[i, ] <- mu0 * w^(n-1) * ddd
}

v <- sqrt(min(w_mat)*max(w_mat))
v_idx <- length(w[w<v])
#n_output <- initial_n*(kappa*w[v_idx]^(-lambda))/sum(sqrt(initial_n[,v_idx]*initial_n[,v_idx+19]))
n_output <- initial_n*(kappa*w[v_idx]^(-lambda))/sum(initial_n[,v_idx])

# more generally, just sum over 1:no_background species

plot(w, kappa*w^(-lambda), type="l", log="xy")
for (i in 1:no_sp) {
  lines(w, n_output[i,], col="blue")
}
lines(w, colSums(n_output), col="red")

# set erepro to meet boundary condition
rdi <- getRDI(params, n_output, params@cc_pp)
gg <- getEGrowth(params, n_output, params@cc_pp)
erepro_final <- rdi
for (i in (1:no_sp)){
  erepro_final[i] <- erepro*(gg[i,params@species_params$w_min_idx[i]]*n_output[i,params@species_params$w_min_idx[i]])/
    rdi[i]
}
params@species_params$erepro <- erepro_final

# check boundary condition is now satisfied
rdi <- getRDI(params, n_output, params@cc_pp)
gg <- getEGrowth(params, n_output, params@cc_pp)
erepro_final <- rdi
rdi_wanted <- rdi
for (i in (1:no_sp)){
  rdi_wanted[i] <- gg[i,params@species_params$w_min_idx[i]]*n_output[i,params@species_params$w_min_idx[i]]
}
rdi_wanted/rdi

############

m2 <- getM2(params, n_output, params@cc_pp)


params2 <- params

for (i in 1:no_sp) {
  params2@psi[i, ] <- (w/w_inf[i])^(1-n)
  params2@psi[i, w < w_mat[i]] <- 0
  params2@psi[i, w > w_inf[i]] <- 1
  params2@mu_b[i, ] <- mu0 * w^(n-1) - m2[i,]
}
#params2@cc_pp <- (kappa*w^(-lambda))-colSums(n_output)
#params2@cc_pp[params2@cc_pp<0] <- 0
plankton_vec <- (kappa*w^(-lambda))-colSums(n_output)
plankton_vec[plankton_vec<0] <- 0
params2@cc_pp[sum(params2@w_full<=w[1]):length(params2@cc_pp)] <- plankton_vec

params2@srr <- function(rdi, species_params) {return(rdi)}

sim <- project(params2, t_max=100 ,effort = 0, initial_n = n_output)
plot(sim)

params@w_full[sum(params@w_full<=w[1])]

plot(params2@w_full, params2@w_full*params2@cc_pp, type="l", log="xy")
lines(params2@w_full, params2@w_full*sim@n_pp[dim(sim@n)[1],], col="blue")

##############

nn <- n_output
nn[nn==0] <- 1
nn

params2@chi <- 0.05
params2@ddd <- nn^(params2@chi)
sim <- project(params2, t_max=100 ,effort = 0, initial_n = n_output)
plot(sim)

