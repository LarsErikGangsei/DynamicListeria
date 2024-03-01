#' Functions for Estimating pH parameter beta2 

#' @param t: Times (in days) for observations, length m.
#' @param pH: observed pH, length m
#' @param Temp: Temperatures in Celsius, length m (or scalar)
#' @param carb: Carbohydrate level,length m, or scalar
#' @param pH0: Start pH, length m or scalar
#' @param pHmin: Minimum pH, scalar
#' @param Temppar: Named list which is passed to "Cardinal_func" for temperature.


#' @return A list with 5 elements \cr
#' - beta2: Estimated value for \eqn{\beta^2}\cr
#'
#' @author  Lars Erik Gangsei
#' 
#' @examples
#' 
#' @export
#' 

Estimation_func_beta2 <- function(t,pH,Temp,carb,pH0,pHmin,Temppar = NULL)
{
  nn <- length(t)
  if(length(Temp)==1){Temp <- rep(Temp,nn)}
  if(is.null(Temppar)){Temppar <- CardinalPar$Temperature}
  g_vec <- log(carb+1)*CM_func(x =Temp,Par = Temppar)
  h_vec <- ((log(pH0-pHmin)-log(pH-pHmin))/(g_vec*t))
  h_vec <- h_vec[((abs(h_vec)<Inf)&(!is.na(h_vec))&!is.nan(h_vec))]
  beta2 <- max(c(0,mean(h_vec)))
}