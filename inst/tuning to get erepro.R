
library("BB")


get_erepro <- function(ks){return(set_scaling_model(rfac = Inf, ks=ks)@species_params$erepro)}

erepro_target <- rep(0.7,11)
ks_start <- rep(4.2,11)
erepro_deviation <- function(ks){return(sum((set_scaling_model(rfac = Inf, ks=ks)@species_params$erepro-erepro_target)^2))}
#op <- optim(par = ks_start, fn = erepro_deviation)
            

#model <- nls(erepro_target ~ get_erepro(ks),start = list(ks = ks_start))

model <- nls(erepro_target ~ get_erepro(ks),start = list(ks = ks_start), control = list(maxiter = 5))




#optim(par, fn, gr = NULL, ...,
#      method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
#                 "Brent")



set.seed(20160227)
x<-seq(0,50,1)
y<-((runif(1,10,20)*x)/(runif(1,0,10)+x))+rnorm(51,0,1)
#for simple models nls find good starting values for the parameters even if it throw a warning
m<-nls(y~a*x/(b+x))



###################
datsB <- data.frame(X=TimesWeKnowG-TimesWeKnowG[1], Y= lengthsB)
param1 <- params
W_egg <- param1@w[1]
L_egg <- (W_egg/a)^(1/b)
gett0 <- function(Linf,k){
  return(log(1-L_egg/Linf)/k)
}
vbTyp<-function(X,Linf,k){(Linf*(1-exp(-k*(X-gett0(Linf,k)))))}
obs_k <- MU[[1]]
obs_Linf <- MU[[2]]
fitTypB<-nls(Y~vbTyp(X,Linf,k),data=datsB,start=list(Linf=obs_Linf,k=obs_k))
vbfitB <- coef(fitTypB)
loglikeB <- dmvnorm(c(vbfitB[["k"]],vbfitB[["Linf"]]),MU,SIGMA,log=T)
loglikeB
#return(loglikeB)}
#################



