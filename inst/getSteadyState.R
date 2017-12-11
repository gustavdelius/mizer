# got a function to calculate the steady state, and we try and fix levels of all species
# but there is some issue

getSteadyState <- function(params){
    res <- params@psi
    for (i in (1:dim(res)[1])){
        mu0 <- (1-params@f0)*sqrt(2*pi)*params@kappa*params@species_params$gamma[i]*
            params@species_params$sigma[i]*
            (params@species_params$beta[i]^(1+params@q-params@lambda))*
            exp((params@species_params$sigma[i]^2)*((1+params@q-params@lambda)^2)/2)
        
        hbar <- params@species_params$alpha[i]*params@species_params$h[i]*
            params@f0-params@species_params$ks[i]
        
        Njuv <- ((params@species_params$w_min[i]/params@w)^(mu0/hbar))/(hbar*params@w^params@n)
        
        pow <- mu0/(hbar*(1-params@n))
        Nmat <- Njuv*
            ((1-(params@w/params@species_params$w_inf[i])^(1-params@n))^(pow-1))*
            (1-(params@species_params$w_mat[i]/params@species_params$w_inf[i])^(1-params@n))^(-pow)
        res[i,] <- ((params@w<params@species_params$w_mat[i])*Njuv+
                        (params@w>=params@species_params$w_mat[i])*Nmat)
    }
    return(res)
}

####################
myalpha <- 0.6
#myalpha <- 0.35
params <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 10+10^(-0.5),
                          w_pp_cutoff = 10,z0pre = 0,no_w=1000,p=2/3,alpha=myalpha)
params@mu_ext <- large_predation(params,WW=params@w[1])
sim <- project(params, t_max=75, effort = 0)
plot(sim)
##################
sol <- getSteadyState(params)

# supposing egg size of species 1 is w[1]
hbar <- params@species_params$alpha*params@species_params$h*
    params@f0-params@species_params$ks
sol_mult <- getRDD(params,sim@n[dim(sim@n)[1],,],sim@n_pp[dim(sim@n_pp)[1],])[1]/
    (sol[1,1]*hbar[1]*params@w[1]^params@n)

plot(params@w,sol_mult*sol[1,],log="xy",type="l")
lines(params@w,sim@n[dim(sim@n)[1],1,],col="blue")

params2 <- params
if (dim(params@psi)[1]>1){
    for (i in (2:dim(params@psi)[1])){
        params2@species_params$r_max[i] <- (hbar[i]*params@species_params$w_min[i]^params@n)*
            sol_mult*sol[i,params@species_params$w_min_idx[i]]
    }
}

sim2 <- project(params2, t_max=75, effort = 0)

plot(params@w,sol_mult*sol[1,],log="xy",type="l")
lines(params@w,sim2@n[dim(sim2@n)[1],1,],col="blue")

plot(params@w,sol_mult*sol[2,],log="xy",type="l")
lines(params@w,sim2@n[dim(sim2@n)[1],2,],col="blue")

# solutions for different w* agree in first part, because 
# here juvenile solutions are the same because egg size is constant

#(hbar[1]*params@species_params$w_min[1]^params@n)*
#    sol_mult*sol[1,params@species_params$w_min_idx[1]]
#params2@species_params$r_max

#params@species_params$r_max

# looks like old rmax was quite a lot higher than rdd

#getRDI(params,sim@n[dim(sim@n)[1],,],sim@n_pp[dim(sim@n_pp)[1],])
#getRDD(params,sim@n[dim(sim@n)[1],,],sim@n_pp[dim(sim@n_pp)[1],])


getRDI(params2,sim2@n[dim(sim2@n)[1],,],sim2@n_pp[dim(sim2@n_pp)[1],])
getRDD(params2,sim2@n[dim(sim@n)[1],,],sim2@n_pp[dim(sim2@n_pp)[1],])

# this approach works when alpha =0.6, but for smaller alpha, like alpha =0.35, 
# we found useful before, the Rdd is too much lower than the Rmax, so it does not work
# also, we should probably make the egg size proportional to maturity size, so 
# the species have different juvenile solutions
