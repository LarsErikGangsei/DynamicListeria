% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_Listeria.R
\name{predict_Listeria}
\alias{predict_Listeria}
\title{Functions for simulating pH and listeria concentrations}
\usage{
predict_Listeria(
  InputListeria,
  pHpar,
  Par_list = NULL,
  Interaction = TRUE,
  log10Return = TRUE
)
}
\arguments{
\item{t:}{time (in days), length n}

\item{Temperature:}{in Celsius}

\item{Acetate:}{Acetatelevel}

\item{Lactate:}{lactatelevel}

\item{aW:}{Wateractivity}

\item{y0:}{Listeria levels at time 0}

\item{mu_opt:}{value for \eqn{\mu_{opt}}}

\item{pHpar:}{List with the following elements
- pH0: \eqn{pH_0}, start pH (pH at time t = 0), see Gangsei et.a. for details.\cr
- pHmin: \eqn{pH_{min}}, the minimum pH value (minimum cardinal pH value)
- beta2: value for \eqn{\beta^2}
- TempparpH: Cardinal parameters for temperature for "pH -bacterias".}

\item{Par_list:}{Named list which is passed to "Cardinal_func". If set to NULL
the values in data "CardinalPar" is used.}

\item{Interaction:}{Passed directly to "Cardinal_func".}

\item{log10Return:}{TRUE (default) or FALSE if results are to be returned at log10 scale.
If TRUE, interactions are calculated based on CM elements.}
}
\value{
Predicted Listeria densities
}
\description{
Functions for simulating pH and listeria concentrations
}
\author{
Lars Erik Gangsei
}
