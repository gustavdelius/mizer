library(mizer)
library(plyr)
paramsConst <- set_trait_model(no_sp = 10, min_w_inf = 10)
sim <- project(paramsConst, t_max=150, effort = 0)
n <- sim@n[dim(sim@n)[1],,]
n_pp <- sim@n_pp[dim(sim@n_pp)[1],]
getEGrowth(paramsConst,n,n_pp)
getZ(paramsConst, n, n_pp, 0)

wvec <- paramsConst@w
myspno <- 3
ns <- n[myspno,]
gs <- getEGrowth(paramsConst,n,n_pp)[myspno,]
ms <- getZ(paramsConst, n, n_pp, 0)[myspno,]

outvec <- function(n){
  gs <- getEGrowth(paramsConst,n,n_pp)[myspno,]
  ms <- getZ(paramsConst, n, n_pp, 0)[myspno,]

  L <- length(wvec)
  dif <- rep(0,L)
  for (i in (2:L)){
    dif[i] <- (gs[i]*ns[i]-gs[i-1]*ns[i-1])/(wvec[i]-wvec[i-1])
  }
  
  ratio <- wvec[2]/wvec[1]
  w0 <- wvec[1]/ratio
  
  dif[1] <- (gs[1]*ns[1]-getRDI(paramsConst, n, n_pp)[myspno])/(wvec[1]-w0)
  
  return(-dif-ms*ns)
    
}






