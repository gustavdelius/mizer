x <- cbind(x1 = 3, x2 = c(4:1, 2:5))

x
colSums(x)



setGeneric('getM2dead', function(object, n, n_pp, n_d, pred_rate, ...)
    standardGeneric('getM2dead'))

setMethod('getM2dead', signature(object='MizerParams', n = 'missing', 
                             n_pp='missing', n_d='missing', pred_rate = 'array'),
          function(object, pred_rate){
              if ((!all(dim(pred_rate) == c(nrow(object@species_params),length(object@w),length(object@w_full)))) | (length(dim(pred_rate))!=3)){
                  stop("pred_rate argument must have 3 dimensions: no. species (",nrow(object@species_params),") x no. size bins (",length(object@w),") x no. size bins in community + background (",length(object@w_full),")")
              }
              idx_sp <- (length(object@w_full) - length(object@w) + 1):length(object@w_full)
              interr <- matrix(1, nrow = dim(object@interaction)[1], ncol = dim(object@interaction)[1]) 
              m2 <- t(interr) %*% colSums(aperm(pred_rate, c(2,1,3)),dims=1)[,idx_sp]
              return(m2)
          }
)


setMethod('getM2', signature(object='MizerParams', n = 'matrix', 
                             n_pp='numeric', n_d='numeric', pred_rate = 'missing'),
          function(object, n, n_pp, n_d){
              noSpecies <- dim(object@interaction)[1]
              muVals <- matrix(0, nrow = noSpecies, ncol = length(object@w))
              w <- object@w
              x <- log(w)
              x <- x - x[1]
              dx <- x[2]-x[1]
              feeding_level <- getFeedingLevel(object, n=n, n_pp=n_pp, n_d=n_d)
              no_P <- length(object@smatM[1,])
              muIntermediate <- matrix(0, nrow = noSpecies, ncol = length(object@w))
              for (j in 1:noSpecies){
                  f <- (1-feeding_level[j,])*object@search_vol[j,]*n[j,]*w
                  f <- c(f[1:length(x)], rep(0, no_P-length(x)))
                  mortalityIntegral <- dx*Re(fft((object@fsmatM[j,])*fft(f), inverse=TRUE)/no_P)
                  muIntermediate[j, ] <- c(mortalityIntegral[(no_P-1):no_P], mortalityIntegral[1:(length(x)-1-1)])
              }
              interr <- matrix(1, nrow = dim(object@interaction)[1], ncol = dim(object@interaction)[1])
              for (i in 1:noSpecies){
                  for (j in 1:noSpecies){
                      muVals[i, ] <- muVals[i, ]+interr[j,i]*muIntermediate[j, ]
                  }
                  
              }	    
              rownames(muVals) <- rownames(object@interaction)
              return(muVals)
          }
)
