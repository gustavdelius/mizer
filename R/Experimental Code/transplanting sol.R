library(mizer)
library(plyr)
paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10,n=2/3,q=0.8,eta=0.25,
                               k0=10^(50),kappa=0.4, alpha=0.17,
                               h=1614.363,f0=0.5
                              )
sim <- project(paramsConst, t_max=150, effort = 0)


#alphaE = Sqrt[2*Pi]*gamma*sigma*(beta^(lambda - 2))*
#  Exp[((lambda - 2)^2)*(sigma^2)/2]
#},
#Module[{hval = ((alphaE*kappa/fbar) - alphaE*kappa)},
# use to determine h
# load in analytic soln
h      

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

