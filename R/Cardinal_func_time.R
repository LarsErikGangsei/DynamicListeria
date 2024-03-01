#' Function calculating gamma value over time

#' @param t: vector of lengt n with times (in days) for when to calculate cardinal value
#' @param Par_list: List with elements (known parameters), same names as XX_list.\cr 
#' -If CM function is to be used names are: xmin,xmax,xopt and n (all scalars)\cr
#' -If SR function is to be used names are: alpha, mic\cr
#'  If Par_list is set to NULL the values in data "CardinalPar" is used.
#' @param pHpar: List with necessary elements to be sent to pH_func()
#' @param Temppar: List with element beta0 to be passed to Constant_func()
#' @param aWpar: List with element beta0 to be passed to Constant_func()
#' @param lacpar: List with element beta0 to be passed to Constant_func()
#' @param acepar: List with element beta0 to be passed to Constant_func()
#' @param deltaLpar: List with element beta0 to be passed to Constant_func()
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
Cardinal_func_time <- function(t,Par_list = NULL,pHpar,carb=NULL,
                               pH0=NULL,Temppar,aWpar,lacpar,
                               acepar,deltaLpar=1,Interaction = TRUE,log=FALSE)
{
  if(is.null(Par_list)){Par_list <- CardinalPar}
  nn <- length(t)
  if(!is.null(carb)){pHpar$carb <- carb}
  if(!is.null(Temppar)){pHpar$temp <- Temppar}
  if(!is.null(pH0)){pHpar$pH0 <- pH0}
  XX_df <- data.frame(Temperature = Constant_func(t,Temppar),
                     pH = pH_func(t,pHpar),
                     aW = Constant_func(t,aWpar),
                     Acetate = Constant_func(t,acepar),
                     Lactate = Constant_func(t,lacpar),
                     DeltaL = Constant_func(t,deltaLpar))
 
  
  XX_list <- split(XX_df,1:nn)
  
  res <- sapply(XX_list,Cardinal_func,Par_list = Par_list,
                Interaction = Interaction,log = log,
                Acid_Param = Par_list_Acid)
  return(res)
}
