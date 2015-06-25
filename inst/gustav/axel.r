# Here I will try to reproduce some of the calculation in the paper
# Rossberg, Axel G., Jennifer E. Houle, and Kieran Hyder. 
# “Stock–recruitment Relations Controlled by Feeding Interactions Alone.” 
# Canadian Journal of Fisheries and Aquatic Sciences 70, no. 10 
# (July 19, 2013): 1447–55. doi:10.1139/cjfas-2012-0531.

# To set up the parameters I modify the set_trait_model() function.
# To give every species its own resource spectrum I have to hack the core
# mizer code. I want to do that hack somewhat general, so I will allow an
# arbitrary number of resource spectra. The number of resource spectra
# will be specified in the new no_pp slot in the MizerParams object. There
# will also be an extra interaction_pp matrix to the MizerParams object that
# gives how species i interacts with resource j.
# I have to give an extra index to the n_pp array in the MizerSim object and 
# then modify all places that update or use n_pp.
# This is all work in progress

#' Sets up parameters for a trait-based model without stock-recruitment relation
#'
#' This functions creates a \code{MizerParams} object similar to that set up
#' by \code{set_trait_model} but with a no artificial stock-recruitment
#' relationship and with a randomly generated species interaction matrix
#' TODO: and with several independent resource spectra.
#'
#' @param no_sp The number of species in the model. The default value is 40. 
#'              About half of these on average is expected to go extinct
#' @param min_w_inf The asymptotic size of the smallest species in the community.
#' @param max_w_inf The asymptotic size of the largest species in the community.
#' @param no_w The number of size bins in the community spectrum.
#' @param min_w The smallest size of the community spectrum.
#' @param max_w The largest size of the community spectrum. 
#'              Default value is the largest w_inf in the community x 1.1.
#' @param min_w_pp The smallest size of the background spectrum.
#' @param no_w_pp The number of the extra size bins in the background spectrum 
#'                (i.e. the difference between the number of sizes bins in the 
#'                community spectrum and the full spectrum).
#' @param w_pp_cutoff The cut off size of the background spectrum. Default value is 1.
#' @param n Scaling of the intake. Default value is 2/3. 
#' @param p Scaling of the standard metabolism. Default value is 0.75. 
#' @param q Exponent of the search volume. Default value is 0.9. 
#' @param eta Factor to calculate \code{w_mat} from asymptotic size.
#' @param r_pp Growth rate of the primary productivity. Default value is 4. 
#' @param kappa Carrying capacity of the resource spectrum. Default value is 0.005. 
#' @param lambda Exponent of the resource spectrum. Default value is (2+q-n). 
#' @param alpha The assimilation efficiency of the community. The default value is 0.6
#' @param ks Standard metabolism coefficient. Default value is 4.
#' @param z0pre The coefficient of the background mortality of the community. 
#'              z0 = z0pre * w_inf ^ (n-1). The default value is 0.6.
#' @param h Maximum food intake rate. Default value is 30.
#' @param beta Preferred predator prey mass ratio. Default value is 100.
#' @param sigma Width of prey size preference. Default value is 1.3.
#' @param f0 Expected average feeding level. Used to set \code{gamma}, the factor for the search volume. The default value is 0.5.
#' @param gamma Volumetric search rate. Estimated using \code{h}, \code{f0} and \code{kappa} if not supplied.
#' @param knife_edge_size The minimum size at which the gear or gears select species. Must be of length 1 or no_sp.
#' @param gear_names The names of the fishing gears. A character vector, the same length as the number of species. Default is 1 - no_sp.
#' @param ... Other arguments to pass to the \code{MizerParams} constructor.
#' @export
#' @return An object of type \code{MizerParams}
#' @seealso \linkS4class{MizerParams}
#' @references K. H. Andersen and M. Pedersen, 2010, Damped trophic cascades driven by fishing in model marine ecosystems. Proceedings of the Royal Society V, Biological Sciences, 1682, 795-802.
#' @examples
#' \dontrun{
#' }
set_axel_model <- function(no_sp = 40,
                            min_w_inf = 10,
                            max_w_inf = 1e5,
                            no_w = 100,
                            min_w = 0.001,
                            max_w = max_w_inf * 1.1,
                            min_w_pp = 1e-10,
                            no_w_pp = round(no_w)*0.3,
                            w_pp_cutoff = 1,
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
                            ...){
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
        erepro = 1 # not used but included out of necessity
    )
    
    # Make species interaction matrix
    # as specified in eq.(2) in Rossberg et.al (2013)
    sigmat <- 2*sqrt(2*log(no_sp))
    interaction <- matrix(rlnorm(no_sp^2, -sigmat^2/2, sigmat),nrow=no_sp, ncol=no_sp)
    
    # Make the MizerParams
    trait_params <- MizerParams(trait_params_df, interaction, min_w = min_w, max_w=max_w, no_w = no_w, min_w_pp = min_w_pp, w_pp_cutoff = w_pp_cutoff, n = n, p=p, q=q, r_pp=r_pp, kappa=kappa, lambda = lambda) 
    
    # Set trivial stock-recruitment relationship
    trait_params@srr <- function(rdi, species_params){
        return(rdi)
    }
    
    return(trait_params)
}

# Try it out
params <- set_axel_model(no_sp = 10, min_w_inf = 10, max_w_inf = 1e5)
# Ignore the warning that the interaction matrix should have entries between 0 and 1
summary(params)
sim <- project(params, t_max=75, effort = 0)
plot(sim)