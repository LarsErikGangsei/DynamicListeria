#' Function for calculating constant environmental factors over time. 
#' 
#' @param t: Vector of length n, representing times for when temperature is evaluated.
#' @param Constantpar: list with element beta0, the constant value to be returned 
#' 
#' @return A vector of length n with all elements equal to beta0. 
#'
#' @author  Lars Erik Gangsei
#' 
#' @references Used in Gangsei, .... 
#' 
#' @examples
#' 
#' @export
Constant_func <- function(t,Constantpar)
{
  nn <- length(t)
  if(is.element(class(Constantpar),c('numeric','integer'))){return(rep(Constantpar,nn))}else{
  return(rep(Constantpar$beta0,nn))}
}
