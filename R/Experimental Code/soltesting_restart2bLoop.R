library(mizer)
library(plyr)
paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10)
simConst <- project(paramsConst, t_max=150, effort = 0)
plot(simConst)

kappaRval <- 0.005
kappaRstarval <- 0.004

itme <- function(Z){
  kappaRval <- Z[1]
  kappaRstarval <- Z[2]
  


myspno <- 3
nval <- 2/3
qval <- 0.9
lambdastar <- 2+qval-nval
chi <- 0.05
beta <- paramsConst@species_params$beta[1]
sigmaval <- paramsConst@species_params$sigma[1] 

# determine phi
wvec <- paramsConst@w
wvecfull <- paramsConst@w_full
determine_phi <- function(w){
  s <- exp(-(log(beta*wvecfull/w))^2/(2*sigmaval*sigmaval))
  integrand <- s*wvecfull*kappaRstarval*wvecfull^(-lambdastar)
  LL <- length(wvecfull)-1
  riemann <- sum(((wvecfull[2:(LL+1)]-wvecfull[1:LL])*integrand[1:LL]))
  return(riemann)
}
#determine_phi(1)

phivals <- sapply(wvec,determine_phi)
gamma <- paramsConst@species_params$gamma[1]
h <- paramsConst@species_params$h[1]
Ee <- gamma*phivals*wvec^qval
feeding_lvl <- Ee/(Ee+h*wvec^nval)
plot(wvec,feeding_lvl,log="x")
# hopefully dip at LHS is caused because finite w resolution cannot 
# account for beta
fbar <- max(feeding_lvl)
#fbar <- feeding_lvl[length(feeding_lvl)-2]

psi <- paramsConst@psi[myspno,]
alpha <- paramsConst@species_params$alpha[myspno]

growth_rate <- (alpha*fbar*h*wvec^nval)*(1-psi)

################
lambda <- (lambdastar-chi)/(1+chi)
alpha_p <- (1-fbar)*sqrt(2*pi)*kappaRval*gamma*sigmaval*(beta^(1+qval-lambda))*
  exp(((sigmaval^2)*(1+qval-lambda)^2)/2)
alpha_p
mu_p <- alpha_p*wvec^(nval-1)

# check this vs numeric integral

determine_mort <- function(w){
  s <- exp(-(log(beta*w/wvecfull))^2/(2*sigmaval*sigmaval))
  integrand <- (1-fbar)*s*kappaRval*(wvecfull^(-lambda))*gamma*wvecfull^qval
  #integrand <- (1-feeding_lvl)*s*kappaRval*(wvec^(-lambda))*gamma*wvec^qval
  
  LL <- length(wvecfull)-1
  riemann <- sum(((wvecfull[2:(LL+1)]-wvecfull[1:LL])*integrand[1:LL]))
  return(riemann)
}

mortvals <- sapply(wvec,determine_mort)

A<-1
mu_pp <- alpha_p*wvec^(A+1+qval-((lambdastar+A)/(chi+1)))


plot(wvec,mortvals,log="xy")
lines(wvec,mu_pp)
# work out mortality integral again, 
#and check whether an adjustment should be added
maxsize <- paramsConst@species_params$w_inf[myspno]
maxsizeindex <- length(maxsize[wvec<maxsize])
realalpha_p <- max((mortvals/(wvec^(A+1+qval-((lambdastar+A)/(chi+1)))))[1:maxsizeindex])
XX <- (mortvals/(wvec^(A+1+qval-((lambdastar+A)/(chi+1)))))[1:maxsizeindex]
plot(XX)
alpha_p
abline(v=maxsize)
realmu_pp <- realalpha_p *wvec^(A+1+qval-((lambdastar+A)/(chi+1)))

plot(wvec,mortvals,log="xy")
lines(wvec,realmu_pp)
abline(v=maxsize)

integrand <- realmu_pp/(growth_rate^(chi+1))

determine_int <- function(w){
  L <- length(wvec[wvec<w])
  within <- sum(((wvec[2:(L+1)]-wvec[1:L])*integrand[1:L]))
  return(within)
}
intvals <- sapply(wvec[1:maxsizeindex],determine_int)

result_n <- function(C){
  npart <- ((chi*(intvals+C))^(-1/chi))/growth_rate[1:maxsizeindex]
  nstart <- rep(0,length(wvec))
  nstart[1:maxsizeindex] <- npart
  return(nstart)
}

plot(wvec,result_n(11),log="xy")

#plot(wvec[1:71],result_n(1)[1:71],log="xy")

repro_eff <- 0.1
egg_size <- wvec[1]

hbar <- alpha*fbar*h

is_zero <- function(C){
  result_nnn <- result_n(C)
  LHS <- result_nnn[1]*growth_rate[1]
  RHS <- sum(((result_nnn*psi*hbar*wvec^nval)[1:(length(wvec)-1)])*(wvec[2:length(wvec)]-wvec[1:(length(wvec)-1)]))*repro_eff/(2*egg_size)
  return(RHS-LHS)
}
library(pracma)
Cval <- newtonRaphson(is_zero, 1)$root
#plot(sapply(seq(1,7,0.1),is_zero))
#is_zero(5)
#is_zero(5.6699)
#optim(1,is_zero,lower = 0.5,upper=8)
plot(wvec[1:71],result_n(1)[1:71],log="xy")
plot(wvec,result_n(Cval),log="xy")
ngood <- result_n(Cval)
# make function that gives a value of ngood for any input weight
FF <- approxfun(x=wvec, y = result_n(Cval),       method = "linear",
                yleft=0, yright=0, rule = 1, f = 0, ties = mean)
plot(wvec,sapply(wvec,FF),log="xy")


Lambda <- (lambdastar+A)/(chi+1)

matsize <- paramsConst@species_params$w_mat[myspno]
Npred <- function(w){
  ws <- wvecfull
  integrand <- ((matsize/ws)^Lambda)*FF(w*(matsize/ws))
  L <- length(ws)-1
  return(sum((ws[2:(L+1)]-ws[1:L])*integrand[1:(L)]))
}

community <- sapply(wvec,Npred)
lines(wvec,wvec^(-lambdastar))
plot(wvec,sapply(wvec,Npred),log="xy")
plot(sapply(wvec,Npred)/(wvec^(-lambda)))
kappaouth <- sapply(wvec,Npred)/(wvec^(-lambda))[1]
kappaout <- kappaouth[1]

Nprey <- function(w){
  ws <- wvecfull
  integrand <- (((matsize/ws)^Lambda)*FF(w*(matsize/ws)))^(chi+1)
  L <- length(ws)-1
  return(sum((ws[2:(L+1)]-ws[1:L])*integrand[1:(L)]))
}

## redo Nprey
y <- wvec
chii <- chi
integrand <- matsize*y^(Lambda*(chii+1))*ngood^(chii+1)*y^(-2)
LL <- length(wvec)-1
kappastaroutput <- sum(integrand[1:LL]*(wvec[2:(LL+1)]-wvec[1:LL]))

## redo Npred
y <- wvec
chii <- 0
integrand <- matsize*y^(Lambda*(chii+1))*ngood^(chii+1)*y^(-2)
LL <- length(wvec)-1
kappaoutput <- sum(integrand[1:LL]*(wvec[2:(LL+1)]-wvec[1:LL]))


plot(sapply(wvec,Nprey)/(wvec^(-lambdastar)))
kappastarouth <- sapply(wvec,Nprey)/(wvec^(-lambdastar))[1]
kappastarout <- kappastarouth[1]

return(c(kappaoutput,kappastaroutput))
}
itme(itme(itme(c(0.005,0.004))))
itme(itme(c(0.005,0.004)))
