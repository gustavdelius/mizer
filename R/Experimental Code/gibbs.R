
library(mizer)
library(FME)
library(ggplot2)
library(reshape2)
library(deSolve)
library(ggplot2)
library(grid)
library(methods)
library(plyr)
library(reshape2)
library("plot3D")
library(rgl)
library("plot3Drgl")
library(optimx)
library(MCMCpack)
set.seed(123)

rinvgamma(1, shape = 3.2, scale= 6)
plot(dinvgamma(1:5, shape = 3.2, scale= 6))
#3.2 is alpha=shape, other is beta=scale

####### posterior of parameter

#### make data choose parameters used


N <- 10
thesd <- 1
thepar <- 0.5
y <- thepar*(1:N)+runif(N, 0, thesd )
y
plot(y)

myshape <- .01
myscale <- .5

inipar <- runif(1,0,1)
inisd <- rinvgamma(1, shape = myshape, scale= myscale)

#### compute sum of squared SS and posterior

SS <- function(apar){
  return(sum((apar*(1:N)-y)^2))
} 

unnormalized_posterior <- function(apar, sd =1){
  return(dunif(apar,0,1)*exp(-SS(apar)/(2*sd^2)))
}

plot((1:100)*0.01,sapply((1:100)*0.01,unnormalized_posterior))



####### metropolis within gibbs
par <- inipar
tsd <- inisd

proposal_sd <- 0.01

T <- 1000
par_vec <- (1:T)
sd_vec <- (1:T)



for (t in (1:T)){
  prop <- rnorm(1,par,proposal_sd)
  u <- runif(1,0,1)
  if (u < unnormalized_posterior(prop,sd=tsd)/unnormalized_posterior(par, sd=tsd)){
    par <- prop
  }
  postshape <- myshape + N/2
  postscale = myscale + SS(par)/2
  tsd <- rinvgamma(1, shape = postshape, scale= postscale)
  par_vec[t] <- par
  sd_vec[t] <- tsd
}

hist(par_vec, prob=TRUE)
res <- 1000
parr <- (1:res)/res
dx <- 0.001
pts <- seq(0,1,dx)
lines(pts,sapply(pts, function(x) unnormalized_posterior(x,sd=tsd))/(dx*sum(sapply(pts, function(x) unnormalized_posterior(x,sd=tsd)))))




#plot(pts,sapply(pts, function(x) unnormalized_posterior(x,sd=tsd))/(dx*sum(sapply(pts, function(x) unnormalized_posterior(x,sd=tsd)))))

