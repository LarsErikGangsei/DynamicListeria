% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Simulation_func.R
\name{Simulation_func}
\alias{Simulation_func}
\title{Functions for simulating pH and listeria concentrations}
\usage{
Simulation_func(Enviorment, Simpar, Par_list = NULL, Interaction = TRUE)
}
\arguments{
\item{Enviorment:}{Data frame of size n x 5 with elements: \cr
- t: Times (in days) for observations.\cr
- Temperature: Temperatures \cr
- Acetate: Added acetate in ?\cr
- Lactate: Added lactate in ?\cr
- aW: Wateractivity. \cr}

\item{Simpar:}{List with the following parameters\cr
- mu_opt: \eqn{\mu_{opt}}, see Gangsei et.a. for details.\cr
- sigma_mu:  \eqn{\sigma_{\mu}^2},see Gangsei et.a. for details.\cr 
- y0: \eqn{y_0} Listeria start consentration.
- beta2: \eqn{\beta^2}, see Gangsei et.a. for details.\cr
- sigma_b: \eqn{\sigma_{\beta}^2}, see Gangsei et.a. for details.\cr
- carb: \eqn{x_c}, carbohydrate level, see Gangsei et.a. for details.\cr
- pH0: \eqn{pH_0}, start pH (pH at time t = 0), , see Gangsei et.a. for details.\cr
- pHmin: \eqn{pH_{min}}, the minimum pH value (minimum cardinal pH value)
- TempPar: Cardinal parametersfor temperature for "pH -bacterias".}

\item{Par_list:}{Named list which is passed to "Cardinal_func". If set to NULL
the values in data "CardinalPar" is used.}

\item{Interaction:}{Passed directly to "Cardinal_func". 
If TRUE, interactions are calculated based on CM elements.}
}
\value{
A data frame of size n x 5 containing elements\cr
- t: Times in days, directly from input "Enviorment"
- beta2: Simulated values for \eqn{\beta^2}
- pH: Simulated pH values\cr
- delta_mu: Random errors for instant listeria growth rate
- int_mu: \eqn{\int_0^t h(t) dt}, the overall growth potential without 
random error.\cr
- y_t: \eqn{y_t}, i.e simulated listeria concentrations at times defined by t.
}
\description{
Functions for simulating pH and listeria concentrations
}
\examples{
# Load data

}
\author{
Lars Erik Gangsei
}
