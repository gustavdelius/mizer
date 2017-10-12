library(deSolve)


sim <- runnit(opout)
L_egg <- (params@w[1]/params_data$a)^(1/params_data$b)
W_egg <- param1@w[1]
L_egg <- (W_egg/a)^(1/b)
gett0 <- function(Linf,k){
  return(log(1-L_egg/Linf)/k)
}
L_inf_from_params <- (params_data$w_inf/params_data$a)^(1/params_data$b)
k_from_params <- params_data$k_vb
gett0(L_inf_from_params,k_from_params)

############################# get simulated time vs length curves

ggB <- getEGrowth(params, sim@n[dim(sim@n)[1],,], sim@n_pp[dim(sim@n_pp)[1],])
sols <- as.list(1:12)
mytimes <- seq(0.5,50,0.1)
weightsB <- matrix(0,nrow=length(mytimes),ncol=dim(sim@n)[2])
lengthsB <- weightsB
L_inf_best_fit <- rep(.7,12)
k_best_fit <- rep(70,12)
for (i in (1:12)){
  # scary try, because least squares will go wrong sometimes
  try({
  gini <- approxfun(params@w, ggB[i,])
  myodefun <- function(t, state, parameters){
    return(list(gini(state)))
  }
  weightsB[,i] <- ode(y = params@w[1], times = mytimes, func = myodefun, parms = 1)[,2]
  lengthsB[,i] <- (weightsB[,i]/params_data$a[i])^(1/params_data$b[i])
  fitTypB<-nls(Y~vbTyp(X,Linf,k),data=data.frame(X=mytimes, Y=lengthsB[,i]),start=list(Linf=L_inf_from_params[i],k=k_from_params[i]))
  L_inf_best_fit[i] <- coef(fitTypB)[["Linf"]]
  k_best_fit[i] <- coef(fitTypB)[["k"]]
  }, silent = T)
}
k_best_fit
plot(mytimes,lengthsB[,1])

