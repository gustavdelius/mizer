library(mizer)
library(plyr)
paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10)
simConst <- project(paramsConst, t_max=150, effort = 0)
plot(simConst)

kappaRval <- 0.005
kappaRstarval <- 0.004
qval <- 0.9

kappaRval <- 5.273115e-11
kappaRstarval <- 2.697401e-18
qval <- 1.5

#kappaRval <- 1.912070e-11 
#kappaRstarval <- 7.758407e-19

go(5.273115e-11, 2.697401e-18,1.5)


myspno <- 3
nval <- 2/3
lambdastar <- 2+qval-nval
chi <- 0.5
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


##h <- paramsConst@species_params$h[1]
Ee <- gamma*phivals*wvec^qval
##feeding_lvl <- Ee/(Ee+h*wvec^nval)

##plot(wvec,feeding_lvl,log="x")
# hopefully dip at LHS is caused because finite w resolution cannot 
# account for beta
##fbar <- max(feeding_lvl)
#fbar <- feeding_lvl[length(feeding_lvl)-2]

##((Ee/fbar)-Ee)/(wvec^nval)

##################
fbar <- 0.5
h <- (((Ee/fbar)-Ee)/(wvec^nval))[1]
##################

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

plot(seq(-5,5,1),sapply(seq(-5,5,1),is_zero))

plot(wvec,result_n(Cval),log="xy")

plot(wvec,result_n(10^100+Cval),log="xy")
is_zero(1)
is_zero(10^100)

rootSolve 
cvalues <- 10^seq(1,100,1)
plot(sapply(cvalues,is_zero))

BCdeviation <- sapply(cvalues,is_zero)

Lc <- length(BCdeviation)-1
cvalues[BCdeviation[2:(Lc+1)]-BCdeviation[1:(Lc)]>0]

midval <-function(v){
  L <- length(v)
  return(v[floor((L-1)/2)+1])
}

midval(cvalues)

newtonRaphson(is_zero, midval(cvalues[BCdeviation[2:(Lc+1)]-BCdeviation[1:(Lc)]>0]))$root


hyper_is_zero<- function(logc)
{return(is_zero(exp(logc)))}
newtonRaphson(is_zero, 1)$root
newtonRaphson(hyper_is_zero, 1)$root

is_zero(1.00964e+11)

is_zero(1.00964e+12)

is_zero(1.00964e+10)

plot(1:200,sapply(1:200,hyper_is_zero))

goodcvals <- cvalues[BCdeviation[2:(Lc+1)]-BCdeviation[1:(Lc)]>0]
Lg <- length(goodcvals)
plot(goodcvals[1:(Lg-1)],sapply(goodcvals[1:(Lg-1)],is_zero),log="x")


newtonRaphson(is_zero, midval(goodcvals[1:(Lg-1)]))$root

is_zero(147550537211)
10^15
is_zero()

log(10^15)
plot(exp(20:40),sapply(exp(20:40),is_zero),log="x")
abline(h=0)

exp(20:40)[sapply(exp(20:40),is_zero)>0][1]
is_zero(2.146436e+14)
is_zero(exp(20:40)[sapply(exp(20:40),is_zero)>0][1])
jj <- sapply(exp(20:40),is_zero)>0
length(jj)

toolow <- exp(20:40)[sapply(exp(20:40),is_zero)<=0]
toolow[length(toolow)]



is_zero(toolow[length(toolow)])
is_zero(exp(20:40)[sapply(exp(20:40),is_zero)>0][1])

########
err <- 0
poscand <- exp((-10):100)
negcand <- -exp((-10):100)

cgroup <- poscand

results <- sapply(cgroup,is_zero)

cgroup[which(results^2==min(results^2))]

is_zero(cgroup[which(results^2==min(results^2))])

cgroup <- poscand
results <- sapply(cgroup,is_zero)
if (length(which(results<0))==length(results)){
  err <- 1
}
if (length(which(results>0))==length(results)){
  err <- 1
}

if (err ==0){
  firstpos <- cgroup[which(results>0)[1]]
  negones <- cgroup[which(results<=0)]
  lastneg <- negones[length(negones)]
  bestguess <- firstpos
  if (is_zero(lastneg)^2<is_zero(bestguess)^2){
    bestguess <- lastneg
  }
}

cgroup[which(results^2==min(results^2))]

is_zero(cgroup[which(results^2==min(results^2))])


bestguess
is_zero(bestguess)

is_zero(2.688117e+48)
2.688117e+43

library(NLRoot)
z<- BFfzero(is_zero, 10^10, 10^30)
z
ff <- function(x){0.1+x^2-x}
z<- BFfzero(ff, -2, 2)

bisection <- function(f, a, b, n = 1000, tol = 1e-7) {
  # If the signs of the function at the evaluated points, a and b, stop the function and return message.
  if (!(f(a) < 0) && (f(b) > 0)) {
    stop('signs of f(a) and f(b) differ')
  } else if ((f(a) > 0) && (f(b) < 0)) {
    stop('signs of f(a) and f(b) differ')
  }
  
  for (i in 1:n) {
    c <- (a + b) / 2 # Calculate midpoint
    
    # If the function equals 0 at the midpoint or the midpoint is below the desired tolerance, stop the 
    # function and return the root.
    if ((f(c) == 0) || ((b - a) / 2) < tol) {
      return(c)
    }
    
    # If another iteration is required, 
    # check the signs of the function at the points c and a and reassign
    # a or b accordingly as the midpoint to be used in the next iteration.
    ifelse(sign(f(c)) == sign(f(a)), 
           a <- c,
           b <- c)
  }
  # If the max number of iterations is reached and no root has been found, 
  # return message and end function.
  print('Too many iterations')
}
bisection(ff,-1,1)

BFfzero(ff,-1,1)
bisection(ff,-1,1)

is_zero(bestguess)

hyper_is_zero(log(bestguess))

BFfzero(hyper_is_zero,log(bestguess/10),log(bestguess*10))

32.43766

BFfzerobetter <- function (f, a, b, num = 10, eps = 1e-05) 
{
  h = abs(b - a)/num
  i = 0
  j = 0
  a1 = b1 = 0
  while (i <= num) {
    a1 = a + i * h
    b1 = a1 + h
    if (f(a1) == 0) {
      print(a1)
      print(f(a1))
    }
    else if (f(b1) == 0) {
      print(b1)
      print(f(b1))
    }
    else if (f(a1) * f(b1) < 0) {
      repeat {
        if (abs(b1 - a1) < eps) 
          break
        x <- (a1 + b1)/2
        if (f(a1) * f(x) < 0) 
          b1 <- x
        else a1 <- x
      }
      print(j + 1)
      j = j + 1
      print((a1 + b1)/2)
      print(f((a1 + b1)/2))
    }
    i = i + 1
  }
  if (j == 0) 
    print("finding root is fail")
  else {print("finding root is successful")
  return((a1 + b1)/2)}
}

ww<-BFfzerobetter(hyper_is_zero,log(bestguess/10),log(bestguess*10))
exp(ww)

www <- BFfzerobetter(is_zero,10^10,10^17)
is_zero(www)
is_zero(www*2)
is_zero(www/2)

is_zero(exp(ww))
bestguess

BFfzerobetter(is_zero,10^10,10^11)

ty <- BFfzerobetter(is_zero,-10^10,10^20)

ty <- BFfzerobetter(is_zero,10^14,10^17)
