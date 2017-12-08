#' Extra Mortality From Large Predators
#'
#' Computes the predation mortality that would be induced by extra predators
#' beyond the maximum size of a fish included in the simulation, under the assumption
#' that these new large predators have an abundance as described by the 
#' idealized power law,
#'
#' @param res the mizer params object.
#' @export
large_predation <- function(res,WW=max(params@species_params$w_inf)){
    betaval <- res@species_params$beta[1]
    lambdaval <- res@lambda
    sigmaval <- res@species_params$sigma[1]
    wvec <- res@w
    gammaval <- res@species_params$gamma[1]
   # WW <- max(params@species_params$w_inf)
    fbar <- res@f0
    kappaval <- res@kappa
    qval <- res@q
    vec <- (-betaval^(1 - lambdaval + qval))*exp((1/2)*(
        (1 - lambdaval + qval)^2)*sigmaval^2)*(-1 + fbar)*
        gammaval*kappaval*sqrt(pi/2)*sigmaval*wvec^(1 - lambdaval + qval)*
        (sqrt(1/sigmaval^2)*sigmaval +
             2*pnorm(((1 - lambdaval + qval)*sigmaval^2 + log(betaval) +
                          log(wvec) -log(WW))/(sqrt(2)*sigmaval),mean=0,sd=sqrt(0.5))-1)
    output <- res@mu_ext
    for (i in 1:dim(res@mu_ext)[1])
    {
        output[i,] <- vec
    }
    return(output)
}


