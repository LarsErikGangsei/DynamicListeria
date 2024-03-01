#' Function for calculating pH-value

#' @param t: Vector of length n, representing times for when pH is evaluated.
#' @param pHpar: A list with the following scalar elements.\cr
#' - beta2 >= 0: pH change rate.\cr
#' - carb: carbohydrate level\cr
#' - temp: temperature\cr
#' - pH0: Start pH (pH at time 0).\cr
#' - pHmin: Cardinal value for minimum pH-value.\cr
#' - TempPar: A list with cardinal parameters for temperature. Passed to CM_func. 
#' If NULL Parameters from data CardinalPar is used.

#' @return A vector of length n with pH observations. 
#'
#' @author  Lars Erik Gangsei
#' 
#' @references Used in Gangsei, .... 
#' 
#' @examples
#' # Load data
#' 
#' @export
pH_func <- function(t,pHpar)
{
  nn <- length(t)
  # Simulate pH
  if(is.null(pHpar$TempPar)){pHpar$TempPar <- CardinalPar$Temperature}
  if(pHpar$beta2==0){pH <- rep(pHpar$pH0,nn)}else{
  pH <- pHpar$pH0 -(pHpar$pH0-pHpar$pHmin)*(1 - 
              exp(-pHpar$beta2*log(pHpar$carb+1)*CM_func(pHpar$temp,pHpar$TempPar)*t))}
  return(pH)
}
