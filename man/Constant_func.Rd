% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Constant_func.R
\name{Constant_func}
\alias{Constant_func}
\title{Function for calculating constant environmental factors over time.}
\usage{
Constant_func(t, Constantpar)
}
\arguments{
\item{t:}{Vector of length n, representing times for when temperature is evaluated.}

\item{Constantpar:}{list with element beta0, the constant value to be returned}
}
\value{
A vector of length n with all elements equal to beta0.
}
\description{
Function for calculating constant environmental factors over time.
}
\references{
Used in Gangsei, ....
}
\author{
Lars Erik Gangsei
}
