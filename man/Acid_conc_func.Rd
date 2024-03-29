% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Acid_conc_func.R
\name{Acid_conc_func}
\alias{Acid_conc_func}
\title{Function which returns the acid concentration in millimolar}
\usage{
Acid_conc_func(XX_list, Par_list, pH = NULL)
}
\arguments{
\item{XX_list:}{Named data.frame with elements "pH", "DeltaL", 
"wHA" and "wSalt", each element containing one vector with the value 
to be evaluated.}

\item{Par_list:}{Named list with elements "pKa", "MHA", "MSalt"}
}
\value{
A scalar i.e. the acid concentration in millimolar
}
\description{
Function which returns the acid concentration in millimolar
}
\examples{
# Load data

}
\author{
Lars Erik Gangsei
}
