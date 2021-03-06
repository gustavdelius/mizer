% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MizerParams-class.R
\name{setPlankton}
\alias{setPlankton}
\title{Set up plankton}
\usage{
setPlankton(
  params,
  rr_pp = NULL,
  cc_pp = NULL,
  r_pp = params@plankton_params[["r_pp"]],
  kappa = params@plankton_params[["kappa"]],
  lambda = params@plankton_params[["lambda"]],
  n = params@plankton_params[["n"]],
  w_pp_cutoff = params@plankton_params[["w_pp_cutoff"]],
  plankton_dynamics = NULL
)
}
\arguments{
\item{params}{A MizerParams object}

\item{rr_pp}{Optional. Vector of plankton intrinsic growth rates}

\item{cc_pp}{Optional. Vector of plankton intrinsic carrying capacity}

\item{r_pp}{Coefficient of the intrinsic growth rate}

\item{kappa}{Coefficient of the intrinsic carrying capacity}

\item{lambda}{Scaling exponent of the intrinsic carrying capacity}

\item{n}{Scaling exponent of the intrinsic growth rate}

\item{w_pp_cutoff}{The upper cut off size of the plankton spectrum. 
Default is 10 g.}

\item{plankton_dynamics}{Function that determines plankton dynamics by
calculating the plankton spectrum at the next time step from the current
state.}
}
\value{
A MizerParams object
}
\description{
Sets the intrinsic plankton growth rate and the intrinsic plankton carrying
capacity as well as the name of the function used to simulate the plankton
dynamics
}
\section{Setting plankton dynamics}{

By default, mizer uses a semichemostat model to describe the plankton
dynamics in each size class independently. This semichemostat dynamics is implemented
by the function \code{\link{plankton_semichemostat}}. You can change the
plankton dynamics by writing your own function, modelled on
\code{\link{plankton_semichemostat}}, and then passing the name of your
function in the \code{plankton_dynamics} argument.

The \code{rr_pp} argument is a vector specifying the intrinsic plankton
growth rate for each size class. If it is not supplied, then the intrinsic growth
rate \eqn{r_p(w)} at size \eqn{w}
is set to \deqn{r_p(w) = r_p\, w^{n-1}.}{r_p(w) = r_p w^{n-1}}
The values of \eqn{r_p} and \eqn{n} are taken from the \code{r_pp}
and \code{n} arguments.

The \code{cc_pp} argument is a vector specifying the intrinsic plankton
carrying capacity for each size class. If it is not supplied, then the intrinsic carrying
capacity \eqn{c_p(w)} at size \eqn{w}
is set to \deqn{c_p(w) = \kappa\, w^{-\lambda}}{c_p(w) = \kappa w^{-\lambda}}
for all \eqn{w} less than \code{w_pp_cutoff} and zero for larger sizes.
The values of \eqn{\kappa} and \eqn{\lambda} are taken from the \code{kappa}
and \code{lambda} arguments.
}

\seealso{
Other functions for setting parameters: 
\code{\link{setBMort}()},
\code{\link{setFishing}()},
\code{\link{setInitial}()},
\code{\link{setIntakeMax}()},
\code{\link{setInteraction}()},
\code{\link{setMetab}()},
\code{\link{setParams}()},
\code{\link{setPredKernel}()},
\code{\link{setReproduction}()},
\code{\link{setSearchVolume}()}
}
\concept{functions for setting parameters}
