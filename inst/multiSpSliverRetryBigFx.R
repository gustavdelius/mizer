# got a function to calculate the steady state, and with many species, over a wide size 
# range, but with low abundance, the results are close but not exact for large size range

getSteadyState <- function(params){
    res <- params@psi
    for (i in (1:dim(res)[1])){
        mu0 <- (1-params@f0)*sqrt(2*pi)*params@kappa*params@species_params$gamma[i]*
            params@species_params$sigma[i]*
            (params@species_params$beta[i]^(1+params@q-params@lambda))*
            exp((params@species_params$sigma[i]^2)*((1+params@q-params@lambda)^2)/2)
        
        hbar <- params@species_params$alpha[i]*params@species_params$h[i]*
            params@f0-params@species_params$ks[i]
        LL <- length(params@w)
        Njuv <- rep(0,LL)
        ei <- params@species_params$w_min_idx[i]
        Njuv[ei:LL] <- ((params@species_params$w_min[i]/params@w[ei:LL])^(mu0/hbar))/(hbar*params@w[ei:LL]^params@n)
        
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
#params <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 10+10^(-0.5),
#                          w_pp_cutoff = 10,z0pre = 0,no_w=1000,p=2/3,alpha=myalpha)
params <- set_trait_model(no_sp = 10, min_w_inf = 10, max_w_inf = 10^5,
                          w_pp_cutoff = 10^5,z0pre = 0,p=2/3,alpha=myalpha)
params@mu_ext <- large_predation(params,WW=params@w[1])

# setup egg sizes
ratioEggMat <- params@species_params$w_min[1]/params@species_params$w_mat[1]
params@species_params$w_min <- params@species_params$w_mat*ratioEggMat
no_s <- dim(params@psi)[1]
for (i in (1:no_s)){
    params@species_params$w_min_idx[i] <- 
        length(params@w[params@w<=params@species_params$w_min[i]])
}


##################
sol <- getSteadyState(params)

# supposing egg size of species 1 is w[1]
hbar <- params@species_params$alpha*params@species_params$h*
    params@f0-params@species_params$ks
sol_mult <- rep(10^(-6),no_s)
true_sol <- sol
for (i in 1:no_s){
    sol_mult[i] <- sol_mult[1]*(params@species_params$w_mat[1]/params@species_params$w_mat[i])^(params@lambda)
    true_sol[i,] <- sol_mult[i]*sol[i,]
    params@species_params$r_max[i] <- true_sol[i,params@species_params$w_min_idx[i]]*
        hbar[i]*params@species_params$w_min[i]^params@n
}

true_sol[is.nan(true_sol)] <- 0
#sim <- project(params, t_max=500, effort = 0)
sim <- project(params, t_max=505, effort = 0,initial_n=true_sol)

plot(sim)

plot(params@w,true_sol[1,],log="xy",type="l")
lines(params@w,sim@n[dim(sim@n)[1],1,],col="blue")

plot(params@w,true_sol[8,],log="xy",type="l")
lines(params@w,sim@n[dim(sim@n)[1],8,],col="blue")


plot(params@w,sim@n[dim(sim@n)[1],1,],log="xy",type="l",ylim=c(10^(-20),10^(-5)))
for (i in 2:dim(sim@n)[2]){
    lines(params@w,sim@n[dim(sim@n)[1],i,],col="blue")
}
##
plot(params@w,true_sol[1,],log="xy",type="l",ylim=c(10^(-30),10^(-5)))
for (i in 2:dim(sim@n)[2]){
    lines(params@w,true_sol[i,],col="blue")
}


#Got muliple species working with variable egg size. H[1]=10^(-6) is used,
#and other H values scale in the same way as N does with w*. 
#Also, the ratio of egg size and maturity size is fixed