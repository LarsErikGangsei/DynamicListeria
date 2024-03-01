#' Functions for predicting pH 

#' @param t: Times (in days) for observations, length m.
#' @param Temp: Temperatures in Celsius, length m (or scalar)
#' @param carb: Carbohydrate level,length m, or scalar
#' @param beta2: value for \eqn{\beta^2}
#' @param pH0: Start pH, length m or scalar
#' @param pH_min: Minimum pH, scalar
#' @param Temppar: Named list which is passed to "Cardinal_func" for temperature.
#' 
#' @return predicted length m
#' 
#' @author  Lars Erik Gangsei
#' 
#' @examples
#' 
#' @export
#' 

predict_pH <- function(t,Temp,carb,pH0,pH_min,beta2,Temppar = NULL)
{
  nn <- length(t)
  if(length(Temp)==1){Temp <- rep(Temp,nn)}
  if(is.null(Temppar)){Temppar <- CardinalPar$Temperature}
  g_vec <- log(carb+1)*CM_func(x =Temp,Par = Temppar)
  return(pH0-(pH0-pH_min)*(1-exp(-beta2*t*g_vec)))
}