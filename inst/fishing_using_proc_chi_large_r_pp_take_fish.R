

library(progress)
library(mizer)
# first run through just sets up the interaction matrix of the right size
f0 <- 0.6
alpha <- 0.4
r_pp <- 1e-1
n <- 2/3
p <- n
q <- 3/4
lambda <- 2+q-n
kappa <- 7e10
erepro <- 0.1 # Will be changed later
beta <- 100
sigma <- 1.3
h <- 30
ks <- 4
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
  erepro = erepro,
  sel_func = "knife_edge", # not used but required
  knife_edge_size = 1000
)

params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
                      kappa = kappa, min_w = min_w, max_w = max_w, no_w = no_w, 
                      min_w_pp = min_w_pp, w_pp_cutoff = max_w, r_pp = r_pp,
                      chi = 0)

############

inter <- params@interaction
inter[,] <- 1

# now we have the interation matrix, we can run the fishing

set_for_fishing <- function(knife_edge_size=10^2,effort=1,inter=inter,r_pp = 1e-1){
  f0 <- 0.6
  alpha <- 0.4
  r_pp <- r_pp
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
  
  
  dist_sp <- 0.2  # separation between species
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
  
  #effort <- 10^(-2)
  # knifeedge=10^2, effort=10^(-2) is too much effort, and causes some species to slowly decline, at least 
  # over 1500 years
  # changing the effort to 10^(-3), maybe its stablizing, but not in 1500 yrs, maybe largest species
  # is going extinct, its hard to tell
 # effort <- 1
  
  
  
  # knife_edge_size <- 10^2
  
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
    knife_edge_size = knife_edge_size
  )
  
  params <- MizerParams(species_params, p=p, n=n, q=q, lambda = lambda, f0 = f0,
                        kappa = kappa, min_w = min_w, max_w = max_w, no_w = no_w, 
                        min_w_pp = min_w_pp, w_pp_cutoff = max_w, r_pp = r_pp,
                        chi = 0,interaction = inter)
  #' gamma is determined by mizerparams. Note that density dependence is currently off
  gamma <- params@species_params$gamma[1]
  w <- params@w
  
  
  # ----
  #' ### Determine analytic solution
  
  
  mu0 <- (1-f0) * sqrt(2*pi) * kappa * gamma * sigma *
    (beta^(n-1)) * exp(sigma^2 * (n-1)^2 / 2)
  hbar <- alpha * h * f0 - ks
  pow <- mu0/hbar/(1-n)
  n_mult <- (1 - (w/w_inf[1])^(1-n))^(pow-1) * (1 - (w_mat[1]/w_inf[1])^(1-n))^(-pow)
  n_mult[w < w_mat[1]] <- 1
  n_mult[w >= w_inf[1]] <- 0
  n_exact <- w  # Just to get array with correct dimensions and names
  n_exact <- ((w_min[1]/w)^(mu0/hbar) / (hbar * w^n)) * n_mult
  n_exact[w < w_min[1]] <- 0
  
  initial_n <- params@psi
  initial_n[,] <- 0
  w_inf_idx <- w_inf
  for (i in 1:no_sp) {
    w_inf_idx[i] <- length(w[w<=w_inf[i]])
    initial_n[i, params@species_params$w_min_idx[i]:
                (params@species_params$w_min_idx[i]+
                   (w_inf_idx[1]-params@species_params$w_min_idx[1]))] <-
      n_exact[params@species_params$w_min_idx[1]:
                (params@species_params$w_min_idx[1]+
                   (w_inf_idx[1]-params@species_params$w_min_idx[1]))] *
      (w_min[1]/w_min[i])^lambda
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
  
  # ----
  #' ### Setup plankton
  
  plankton_vec <- (kappa*w^(-lambda))-colSums(n_output)
  plankton_vec[plankton_vec<0] <- 0
  plankton_vec[min(which(plankton_vec==0)):length(plankton_vec)] <- 0
  params@cc_pp[sum(params@w_full<=w[1]):length(params@cc_pp)] <- plankton_vec
  initial_n_pp <- params@cc_pp
  # The cc_pp factor needs to be higher than the desired steady state in
  # order to compensate for predation mortality
  m2_background <- getM2Background(params, n_output, initial_n_pp)
  params@cc_pp <- (1+m2_background/params@rr_pp) * initial_n_pp
  
  # ----
  #' ### Setup background death and steplike psi
  
  m2 <- getM2(params, n_output, initial_n_pp)
  
  for (i in 1:no_sp) {
    params@psi[i, ] <- (w/w_inf[i])^(1-n)
    params@psi[i, w < (w_mat[i]-1e-10)] <- 0
    params@psi[i, w > (w_inf[i]-1e-10)] <- 1
    params@mu_b[i, ] <- mu0 * w^(n-1) - m2[i,]
  }
  
  # ----
  #' ### Set erepro to meet boundary condition
  
  rdi <- getRDI(params, n_output, initial_n_pp)
  gg <- getEGrowth(params, n_output, initial_n_pp)
  mumu <- getZ(params, n_output, initial_n_pp, effort = 0)
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
  
  # ----
  #' ### Simulate without Rmax or chi
  
  params@srr <- function(rdi, species_params) {return(rdi)}
  params@chi <- 0.0
  
  nn <- n_output
  nn[nn==0] <- 1
  params@chi <- 0.05
  params@ddd <- nn^(params@chi)
  
  
  #li[[1]] <- params
  #li[[2]] <- n_output
  li <- list("params"=params,"n"=n_output,"n_pp"=initial_n_pp,"mu0"=mu0)
  return(li)
}


###

effort <- 1
knife_edge_size <- 10^2
st <- set_for_fishing(knife_edge_size=knife_edge_size,effort=effort,inter=inter,r_pp = 1e-1)
t_max <- 30
t_save <- t_max/100
sim <- project(st$params, t_max=t_max, dt=0.05, t_save=t_save ,effort = effort, 
               initial_n = st$n, initial_n_pp = st$n_pp)
plot(sim)
times <- seq(0,t_max,length=101)
GY <- getYield(sim)

#####

effort <- 10^(0)
#knife_edge_size <- 0.5*10^2
# knife_edge_size <- 600 coexistance with effort=10
knife_edge_size <- 100

t_max <- 300
st <- set_for_fishing(knife_edge_size=knife_edge_size,effort=effort,inter=inter,r_pp = 1e-1)
#t_max <- 5
t_save <- t_max/100
sim2 <- project(st$params, t_max=t_max, dt=0.05, t_save=t_save ,effort = effort, 
               initial_n = st$n, initial_n_pp = st$n_pp)
#plot(sim2)
plotBiomass(sim2)
times <- seq(0,t_max,length=101)
GY2 <- getYield(sim2)


####

GY[dim(GY)[1],]/GY2[dim(GY)[1],]


#plotYield(sim)
i <- 1
plot(times,GY[,1],type="l",col=gray(1-(0.05+i/no_sp)),ylim=c(min(GY),max(GY)))
for (i in (2:no_sp)){
  lines(times,GY[,i],col=gray(1-(i/no_sp)))
  
}
# using a fishing gear with knife_edge_size = 100 gives GY. This gives better yields 
# for species 6,7,8, then GY2 which uses the smaller knife_edge_size = 50 
# the yield is down on species with w_mat close to lower knife edge=50, or between 50 & 100
# but larger knife edge gives better yield on 

i <- 1
plot(times,GY2[,1],type="l",col=gray(1-(0.05+i/no_sp)),ylim=c(min(GY2),max(GY2)))
for (i in (2:no_sp)){
  lines(times,GY2[,i],col=gray(1-(i/no_sp)))
  
}


#i <- 1
#plot(times,GY[,1]/GY2[,1],type="l",col=gray(0),ylim=c(0.01,100))
#for (i in (2:no_sp)){
#  if (GY[1,i]>0){
#  lines(times,GY[,i]/GY2[,i],col=gray(1-(i/no_sp)))
#  }
#}

plot(times,GY[,no_sp]/GY2[,no_sp],col=gray(1-(i/no_sp)),ylim=c(0.01,10))
for (i in ((no_sp-1):7)){
    lines(times,GY[,i]/GY2[,i],col=gray(1-(i/no_sp)))
}



#############################
params <- st$params
w <- st$params@w

speci <- 1
gNS <- getEGrowth(params, sim@n[dim(sim@n)[1],,], sim@n_pp[dim(sim@n_pp)[1],])[speci,]

g_fnNS <- approxfun(w, gNS)
myodefunNS <- function(t, state, parameters){
  return(list(g_fnNS(state)))
}


ageNS <- (0:60)

library(deSolve)
weightNS <- ode(y = params@species_params$w_min[speci], times = ageNS, func = myodefunNS, parms = 1)[,2]
plot(ageNS,weightNS,type="l",ylim=c(0,max(w_inf)))
for (speci in (2:no_sp)){
  gNS <- getEGrowth(params, sim@n[dim(sim@n)[1],,], sim@n_pp[dim(sim@n_pp)[1],])[speci,]
  
  g_fnNS <- approxfun(w, gNS)
  myodefunNS <- function(t, state, parameters){
    return(list(g_fnNS(state)))
  }
  weightNS <- ode(y = params@species_params$w_min[speci], times = ageNS, func = myodefunNS, parms = 1)[,2]
  lines(ageNS,weightNS,type="l")
}

#########################
########################

internull <- inter
internull[,] <- 0


effort <- 1
knife_edge_size <- 10^2
st <- set_for_fishing(knife_edge_size=knife_edge_size,effort=effort,inter=internull,r_pp = 10^12)

plankton_vec <- (kappa*w^(-lambda))
plankton_vec[plankton_vec<0] <- 0
#plankton_vec[min(which(plankton_vec==0)):length(plankton_vec)] <- 0
st$params@cc_pp[sum(st$params@w_full<=w[1]):length(st$params@cc_pp)] <- plankton_vec
initial_n_pp <- st$params@cc_pp
# The cc_pp factor needs to be higher than the desired steady state in
# order to compensate for predation mortality
m2_background <- getM2Background(st$params, st$n, initial_n_pp)
st$params@cc_pp <- (1+m2_background/st$params@rr_pp) * initial_n_pp


m2 <- getM2(st$params, st$n, initial_n_pp)

for (i in 1:no_sp) {
    st$params@mu_b[i, ] <- st$mu0 * w^(n-1) 
}

t_max <- 30
t_save <- t_max/100
sim <- project(st$params, t_max=t_max, dt=0.05, t_save=t_save ,effort = effort, 
               initial_n = st$n, initial_n_pp = st$n_pp)
plot(sim)
times <- seq(0,t_max,length=101)
GYIN <- getYield(sim)

########## plot (yields in interacting model)/(yields in non-interacting model)

plot(times,GY[,no_sp]/GYIN[,no_sp],col=gray(1-(i/no_sp)),ylim=c(0.1,1.2),type="l")
for (i in ((no_sp-1):7)){
  lines(times,GY[,i]/GYIN[,i],col=gray(1-(i/no_sp)))
}
GY[dim(GY)[1],]/GYIN[dim(GYIN)[1],]


# turn on chi to stablize things, and plot relative yields on the two models

#' turned on chi dependence, and compare fishing with knife edge =100g to knife edge=50g. 
#' # Also, I compare the yield predictions under the interacting and non-interacting
#' # models. I hope I setup the non-interacting case properly. I turned on chi to promote 
#' # stability
#' 
#' 
#' now in sim2 we seem to see sustainable fishing. The grid can be setup by this..
#' dist_sp <- 0.4
#minimum_egg <- -4
#maximum_egg <- -2
#no_sp <- (maximum_egg-minimum_egg)/dist_sp + 1
#species <- 1:no_sp
#x_min <- seq(minimum_egg, by = dist_sp, length.out = no_sp)
#w_min <- 10^x_min
#w_inf <- 10^(x_min+5)
#w_mat <- 10^(x_min+4.4)  # This is about a quarter of w_inf
#min_w <- min(w_min)
#max_w <- max(w_inf)
#no_w <- log10(max_w/min_w)*100+1
#min_w_pp <- 10^(minimum_egg-4)
