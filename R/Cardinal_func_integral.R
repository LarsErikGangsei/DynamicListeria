#' Function calculating gamma value over time

#' @param Environment: List with the following elements (scalars)\cr
#' - t: Time\cr
#' - Temp: : one numeric value, to be passed to Constant_func() for temperature.\cr
#' - aW: one numeric value, to be passed to Constant_func() for aW.\cr
#' - Lactate: one numeric value, to be passed to Constant_func() for lactate.\cr
#' - Acetate: one numeric value, to be passed to Constant_func() for acetate.\cr
#' - deltaL: one numeric value, to be passed to Constant_func() for deltaL.\cr
#' - carb: one numeric value for the carbohydrate level
#' @param pHpar: List with necessary elements to be sent to pH_func()
#' @param Par_list: List with elements (known parameters), same names as XX_list.\cr 
#' -If CM function is to be used names are: xmin,xmax,xopt and n (all scalars)\cr
#' -If SR function is to be used names are: alpha, mic\cr
#'  If Par_list is set to NULL the values in data "CardinalPar" is used.
#' @param Interaction: - If TRUE, interactions are calculated based on CM elements
#' @param log: If the cardinal values are to be returned on log scale

#' @return A vector of length n with elements >=0 and <=1, reflecting cardinal 
#' values at times defined by input argument t.
#'
#' @author  Lars Erik Gangsei
#' 
#' @examples
#' Load data
#' 
#' @export
Cardinal_func_integral <- function(Environment,pHpar,Par_list = NULL,
                                   Interaction = TRUE,log=FALSE)
{
  if(Environment$t==0){res <- 0}else{
  if(is.null(Par_list)){Par_list <- CardinalPar}
  res <- integrate(f= Cardinal_func_time,lower = 0,upper = Environment$t,
                   Par_list = Par_list,pHpar=pHpar,Temppar=Environment$Temperature,
                   carb = Environment$carb,pH0 =Environment$pH0, 
                   aWpar=Environment$aW,lacpar = Environment$Lactate,
                   acepar = Environment$Acetate,deltaLpar = Environment$DeltaL,
                   Interaction = Interaction,log = log)
  if(res$message=='OK'){res <- res$value}else{res <- NA}}
  return(res)
}
