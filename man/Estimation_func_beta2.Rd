% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estimation_func_beta2.R
\name{Estimation_func_beta2}
\alias{Estimation_func_beta2}
\title{Functions for simulating pH and listeria concentrations}
\usage{
Estimation_func_beta2(t, pH, Temp, carb, pH0, pHmin, Temppar = NULL)
}
\arguments{
\item{t:}{Times (in days) for observations, length m.}

\item{Temp:}{Temperatures in Celsius, length m (or scalar)}

\item{pH:}{observed pH, length m}

\item{carb:}{Carbohydrate level,length m, or scalar}

\item{Temppar:}{Named list which is passed to "Cardinal_func" for temperature.}

\item{pH0:}{Start pH, length m or scalar}

\item{pHmin:}{Minimum pH, scalar}
}
\value{
A list with 5 elements \cr
- beta2: Estimated value for \eqn{\beta^2}\cr
}
\description{
Functions for simulating pH and listeria concentrations
}
\author{
Lars Erik Gangsei
}
