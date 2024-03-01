#' Functions for predicting Listeria concentrations 

#' @param InputListeria: data.frame (n x 10) with observations for Listeria containing at least:
#' the following elements\cr
#' - t: time (in days)\cr
#' - Temperature: in Celsius\cr
#' - Acetate: Acetate level\cr
#' - Lactate: lactate level\cr
#' - aW: Wateractivity\cr
#' - carb: \eqn{x_c}, carbohydrate level, see Gangsei et.a. for details.\cr
#' - DeltaL: Free water in the matrix.\cr
#' - y0: Listerialevels at time 0
#' - pH0: Start pH.\cr
#' - mu_opt: value for \eqn{\mu_{opt}}
#' @param pHpar: List with the following elements
#' - pH0: \eqn{pH_0}, start pH (pH at time t = 0), see Gangsei et.a. for details.\cr
#' - pHmin: \eqn{pH_{min}}, the minimum pH value (minimum cardinal pH value)
#' - beta2: value for \eqn{\beta^2}
#' - TempparpH: Cardinal parameters for temperature for "pH -bacterias".
#' @param Par_list: Named list which is passed to "Cardinal_func". If set to NULL
#' the values in data "CardinalPar" is used.
#' @param Interaction: Passed directly to "Cardinal_func". 
#' @param log10Return: TRUE (default) or FALSE if results are to be returned at log10 scale.
#' If TRUE, interactions are calculated based on CM elements.
#' 
#' @return Predicted Listeria densities
#'
#' @author  Lars Erik Gangsei
#' 
#' @examples
#' 
#' @export
#' 

predict_Listeria <- function(InputListeria,
                             pHpar,Par_list = NULL,
                             Interaction=TRUE,log10Return = TRUE)
{
  # Basis for estimation of mu_opt based on dynamic enviorment
  nn <- dim(InputListeria)[1]
 
  if(is.null(Par_list)){Par_list <- CardinalPar}
  
  int_vec <-  sapply(split(dplyr::select(InputListeria,
                      t,carb,Temperature,aW,Lactate,Acetate,DeltaL,pH0),1:nn),
                     Cardinal_func_integral,
                     pHpar=pHpar,Par_list= Par_list,
                     Interaction = Interaction,log=FALSE)
  
log10Res <- log10(InputListeria$y0) + InputListeria$mu_opt*int_vec

if(log10Return==TRUE){return(log10Res)}else(return(10^log10Res))
}