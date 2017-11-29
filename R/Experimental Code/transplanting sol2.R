library(mizer)
library(plyr)
source("R/MizerParams-classETRN.R")
source("R/wrapper_functionsETRN.R")

paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10,n=2/3,q=0.8,eta=0.25,
                               k0=10^(50),kappa=0.4, alpha=0.17,
                               h=1614.363,f0=0.5,ks=0,z0=0,gamma=660.2633
                              )
plot(paramsConst@w,paramsConst@mu_et,log="x") 
abline(v=max(paramsConst@w)/100)
# have to supply h and gamma, in order for system to be sensitive to varying f0
paramsConst@mu_et[1:7]
#lambda,WW,f0,kappa,qval,n
sim <- project(paramsConst, t_max=150, effort = 0)
plot(sim)
paramsConst@species_params$z0

# add epro
# make new version of set_trait_model where h is computed, and erepro is added
# make gridpoints align with chosen w* and Wmax
# add mortality term for larger species
# set cutoff of plankton at the maturity (max?) weight of the smallest species 

# z0=k=0
# add erepro
# turn off mort

#alphaE = Sqrt[2*Pi]*gamma*sigma*(beta^(lambda - 2))*
#  Exp[((lambda - 2)^2)*(sigma^2)/2]
#},
#Module[{hval = ((alphaE*kappa/fbar) - alphaE*kappa)},
# use to determine h
# load in analytic soln

################################

kappa = 0.4;
#q = 0.8;
#n = 2/3;
beta = 100;
sigma = 1.3;
gamma = 600.4424;
fbar = 0.5;
k = 0;
we = 0.001;
eta = 1/4;
u = 10;
matsize = 19.35659;
alpha = 0.17

lambda <- 2+0.8-2/3
lambda

###########################

set_trait_model2 <- function(no_sp = 10,
                            min_w_inf = 10,
                            max_w_inf = 1e5,
                            no_w = 100,
                            min_w = 0.001,
                            max_w = max_w_inf * 1.1,
                            min_w_pp = 1e-10,
                            no_w_pp = NA,
                            w_pp_cutoff = 1,
                            k0 = 50, # recruitment adjustment parameter
                            n = 2/3,
                            p = 0.75,
                            q = 0.9, 
                            eta = 0.25,
                            r_pp = 4,
                            kappa = 0.005,
                            lambda = 2+q-n,
                            alpha = 0.6,
                            ks = 4,
                            z0pre = 0.6,
                            h = 30,
                            beta = 100,
                            sigma = 1.3,
                            f0 = 0.5,
                            gamma = NA,
                            knife_edge_size = 1000,
                            gear_names = "knife_edge_gear",
                            erepro2 = 0.1,
                            ...){
  if (!is.na(no_w_pp))
    warning("New mizer code does not support the parameter no_w_pp")
  # If not supplied, calculate gamma using equation 2.1 in A&P 2010
  if(is.na(gamma)){
    alpha_e <- sqrt(2*pi) * sigma * beta^(lambda-2) * exp((lambda-2)^2 * sigma^2 / 2) # see A&P 2009
    gamma <- h * f0 / (alpha_e * kappa * (1-f0)) # see A&P 2009 
  }
  w_inf <- 10^seq(from=log10(min_w_inf), to = log10(max_w_inf), length=no_sp)
  w_mat <- w_inf * eta
  
  # Check gears
  if (length(knife_edge_size) > no_sp){
    stop("There cannot be more gears than species in the model")
  }
  if ((length(knife_edge_size) > 1) & (length(knife_edge_size) != no_sp)){
    warning("Number of gears is less than number of species so gear information is being recycled. Is this what you want?")
  }
  if ((length(gear_names) != 1) & (length(gear_names) != no_sp)){
    stop("Length of gear_names argument must equal the number of species.")
  }
  
  # Make the species parameters data.frame
  trait_params_df <- data.frame(
    species = 1:no_sp,
    w_inf = w_inf,
    w_mat = w_mat,
    h = h, # max food intake
    gamma = gamma, # vol. search rate,
    ks = ks,# standard metabolism coefficient,
    beta = beta,
    sigma = sigma,
    z0 = z0pre * w_inf^(n-1), # background mortality
    alpha = alpha,
    #r_max = r_max,
    sel_func = "knife_edge",
    knife_edge_size = knife_edge_size,
    gear = gear_names,
    erepro = erepro2 # not used but included out of necessity
  )
  # Make the MizerParams
  trait_params <- MizerParams(trait_params_df, min_w = min_w, max_w=max_w, no_w = no_w, min_w_pp = min_w_pp, w_pp_cutoff = w_pp_cutoff, n = n, p=p, q=q, r_pp=r_pp, kappa=kappa, lambda = lambda) 
  # Sort out maximum recruitment - see A&P 2009
  # Get max flux at recruitment boundary, R_max
  # R -> | -> g0 N0
  # R is egg flux, in numbers per time
  # Actual flux at recruitment boundary = RDD = NDD * g0 (where g0 is growth rate)
  # So in our BH SRR we need R_max comparable to RDI (to get RDD)
  # R_max = N0_max * g0 (g0 is the average growth rate of smallest size, i.e. at f0 = 0.5)
  # N0 given by Appendix A of A&P 2010 - see Ken's email 12/08/13
  # Taken from Ken's code 12/08/13 - equation in paper is wrong!
  alpha_p <- f0 * h * beta^(2 * n - q - 1) * exp((2 * n * (q - 1) - q^2 + 1) * sigma^2 / 2)
  alpha_rec <- alpha_p / (alpha * h * f0 - ks)
  # Calculating dw using Ken's code - see Ken's email 12/08/13
  tmpA <- w_inf[1]
  tmpB <- (log10(w_inf[length(w_inf)]) - log10(w_inf[1])) / (no_sp - 1) # Difference between logged w_infs, fine
  dw_winf <- tmpB * tmpA *10^(tmpB*((1:no_sp)-1)) # ?
  N0_max <- k0 * w_inf^(n*2-q-3+alpha_rec) * dw_winf  # Why * dw_winf, not / ? Ken confirms * in email
  # No need to include (1 - psi) in growth equation because allocation to reproduction at this size = 0, so 1 - psi = 1
  g0 <- (alpha * f0 * h * trait_params@w[1]^n - ks * trait_params@w[1]^p)
  r_max <- N0_max * g0
  
  trait_params@species_params$r_max <- r_max
  
  return(trait_params)
}

#############################

paramsConst2 <- set_trait_model2(no_sp = 10, min_w_inf = 10,n=2/3,q=0.8,eta=0.25,
                               k0=10^(50),kappa=0.4, alpha=0.17,
                               h=1614.363,f0=0.5,ks=0,z0=0,erepro2 = 0.1
)
sim2 <- project(paramsConst2, t_max=150, effort = 0)
plot(sim2)

# it seems that we should let paramsConst2 determine gamma via (3.16)
# I should make code that uses alphaE and f0 to determine h, and then gamma

#(1-fbar)*exp(-((log(ww)-log(w)-log(beta))^2)/(2*sigmaval*sigmaval))*gamma*kappa*ww^(q-lambda)

tail(wvec)
logdiff <- log(wvec)[2] -log(wvec)[1]
bigenough <- wvec[length(wvec)]*beta*exp(5*sigmaval)
wvecright <- exp(seq(log(wvec[length(wvec)]),log(bigenough),by=logdiff))
Lr <- length(wvecright)-1


mort_eter <- function(w){
  return(sum(((1-fbar)*exp(-((log(wvecright[1:Lr])-log(w)-log(beta))^2)/(2*sigmaval*sigmaval)
             )*gamma*kappaRval*wvecright[1:Lr]^(qval-lambda))*(wvecright[2:(Lr+1)]-wvecright[1:Lr])))
}
plot(wvec,sapply(wvec,mort_eter),log="xy")

mort_eter <- function(w){
  return(sum(((1-fbar)*exp(-(((log(wvecright[1:Lr])-log(w)-log(beta))^2)/(2*sigmaval*sigmaval))
  )*gamma*kappaRval*wvecright[1:Lr]^(qval-lambda))*(wvecright[2:(Lr+1)]-wvecright[1:Lr])))
}
plot(wvec,sapply(wvec,mort_eter),log="x")

wvecL <- c(wvec,wvecright)

plot(wvecL,sapply(wvecL,mort_eter),log="x")

library(pracma)

direct_mort_eter <- function(w){
  X <- log(w)+log(beta)
  A <- qval-lambda-1
  rett <- 
    gamma*kappaRval*(1-fbar)*sigmaval*sqrt(pi/2)*exp((((A*sigmaval)^2)/2)+A*X)*(
      1-erf((A*sigmaval*sigmaval+X-log(wvec[length(wvec)]))/(sigmaval*sqrt(2)))
    )
  return(rett)
}

(A*sigmaval*sigmaval+X-log(wvec[length(wvec)]))/(sigmaval*sqrt(2))

plot(wvec,sapply(wvec,direct_mort_eter),log="xy")

sum(((1-fbar)*exp(-((log(wvecright[1:Lr])-log(w)-log(beta))^2)/(2*sigmaval*sigmaval)
)*gamma*kappaRval*wvecright[1:Lr]^(qval-lambda))*(wvecright[2:(Lr+1)]-wvecright[1:Lr]))
vw <- 1
plot(((1-fbar)*exp(-((log(wvecright[1:Lr])-log(vw)-log(beta))^2)/(2*sigmaval*sigmaval)
)*gamma*kappaRval*wvecright[1:Lr]^(qval-lambda))*(wvecright[2:(Lr+1)]-wvecright[1:Lr]),log="xy")

vw <-10^7
plot(wvecright[1:Lr],(exp(-((log(wvecright[1:Lr])-log(vw)-log(beta))^2)/(2*sigmaval*sigmaval))),log="xy")

plot(wvecright[1:Lr],exp(-((log(10^7)-log(wvecright[1:Lr]))^2)/(2*sigmaval*sigmaval)),log="xy")

plot(wvecright[1:Lr],exp(-((log(10^7)-log(wvecright[1:Lr]))^2)/(1)),log="xy")


plot(wvecright[1:Lr],exp(-(1-log(wvecright[1:Lr]))^2),log="x")
plot(wvec,exp(-(10-wvec)^2),log="xy")
plot(wvec,exp(-((log(10)-log(wvec))^2)),log="x")
