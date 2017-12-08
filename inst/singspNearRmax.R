
# alpha = 0.6
# sweet spot is alpha in (0.35,0.4)
params <- set_trait_model(no_sp = 1, min_w_inf = 10, max_w_inf = 10+10^(-10),
                          w_pp_cutoff = 10,z0pre = 0,no_w=1000,p=2/3,alpha=0.35)
params@mu_ext <- large_predation(params,WW=params@w[1])
sim <- project(params, t_max=75, effort = 0)
plot(sim)

params@species_params$r_max
getRDI(params,matrix(sim@n[dim(sim@n)[1],,],nrow = 1),sim@n_pp[dim(sim@n_pp)[1],])
getRDD(params,matrix(sim@n[dim(sim@n)[1],,],nrow = 1),sim@n_pp[dim(sim@n_pp)[1],])

plot(params@w,params@mu_ext[1,],log="x")
params@w

# max_w_inf largest gridpoint of w
# min_w_inf = max size of our species
# Delta = max/min
# Rmax prop log(Delta-1) ?
# Rmax prop log(max_w_inf/min_w_inf)
# Abundance scales with Rmax
mu0 <- (1-params@f0)*sqrt(2*pi)*params@kappa*params@species_params$gamma[1]*
    params@species_params$sigma[1]*
    (params@species_params$beta[1]^(1+params@q-params@lambda))*
    exp((params@species_params$sigma[1]^2)*((1+params@q-params@lambda)^2)/2)

plot(params@w,mu0*params@w^(params@n-1),log="xy")
lines(params@w,getM2(params,matrix(sim@n[dim(sim@n)[1],,],nrow = 1),sim@n_pp[dim(sim@n_pp)[1],]),col="red")
lines(params@w,getZ(params,matrix(sim@n[dim(sim@n)[1],,],nrow = 1),sim@n_pp[dim(sim@n_pp)[1],],effort = 0),col="green")

lines(params@w,params@mu_ext[1,])

params@mu_ext <- large_predation(params,WW=params@w[1])


plot(params@w,params@mu_ext[1,],log="x")

getM2

## now max_w_inf/min_w_inf is small, Rmax is small, and so N is small,
# so new species has low abundance, so we set the min size of our eternal predators
# to be the egg size. Plots seem to indicate eternal predation mortality has correct form

params@species_params

hbar <- params@species_params$alpha[1]*params@species_params$h[1]*
    params@f0-params@species_params$ks[1]
plot(params@w,(1-params@psi[1,])*hbar*params@w^params@n,log="xy", type="l", col="blue")
lines(params@w,getEGrowth(params,matrix(sim@n[dim(sim@n)[1],,],nrow = 1),sim@n_pp[dim(sim@n_pp)[1],]))

##
Njuv <- ((params@w[1]/params@w)^(mu0/hbar))/(hbar*params@w^params@n)
plot(params@w,sim@n[dim(sim@n)[1],,],log="xy")
lines(params@w,Njuv)
sol_mult <- getRDD(params,matrix(sim@n[dim(sim@n)[1],,],nrow = 1),sim@n_pp[dim(sim@n_pp)[1],])/
    (Njuv[1]*hbar*params@w[1]^params@n)

plot(params@w,sol_mult*Njuv,log="xy")

plot(params@w,sim@n[dim(sim@n)[1],,],log="xy")
lines(params@w,sol_mult*Njuv)
pow <- mu0/(hbar*(1-params@n))
Nmat <- Njuv*
    ((1-(params@w/params@species_params$w_inf[1])^(1-params@n))^(pow-1))*
    (1-(params@species_params$w_mat[1]/params@species_params$w_inf[1])^(1-params@n))^(-pow)
lines(params@w,sol_mult*Nmat)

Nsol <- sol_mult*((params@w<params@species_params$w_mat[1])*Njuv+
                      (params@w>=params@species_params$w_mat[1])*Nmat)

plot(params@w,Nsol,log="xy",type="l")
lines(params@w,sim@n[dim(sim@n)[1],,],col="blue")

###
plot(params@w,hbar*params@w^params@n,log="xy", type="l", col="blue")
lines(params@w,getEReproAndGrowth(params,
                                  matrix(sim@n[dim(sim@n)[1],,],nrow = 1),
                                  sim@n_pp[dim(sim@n_pp)[1],],
                                  feeding_level = matrix(.6,nrow=1,ncol=length(params@w))))
params@n
plot(params@w,params@intake_max,log="xy")
lines(params@w,params@species_params$h[1]*params@w^params@n, col="Red")
plot(params@w,Nsol,log="xy",type="l")
lines(params@w,sim@n[dim(sim@n)[1],,],col="blue")
