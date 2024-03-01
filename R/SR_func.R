#' Function calculating SR value

#' @param x: scalar. The environmental factor (typically temperature, pH or water activity)
#' to be evaluated.
#' @param Par: List with the following elements "mic" and "alpha" (all scalars) or list with 
#' elements "xmin", "xopt"

#' @return A scalar >=0 and <=1. 
#'
#' @author  Lars Erik Gangsei
#' 
#' @references Augustin, J., Zuliani, V., Cornu, M., and Guillier, L. (2005). 
#' Growth rate and growth probability ofListeria monocytogenes in dairy, meat and seafood 
#' products in suboptimal conditions. Journal of Applied Microbiology, 99(5):1019-1042.
#' doi: https://doi.org/10.1111/j.1365-2672.2005.02710.x 
#' 
#' @examples
#' # Load data
#' 
#' @export
SR_func <- function(x,Par)
{
  #print(x)
  if(class(x)=='data.frame'){x <- x[,1]}
  x <- as.numeric(x)
  #print(Par)
  if(is.element('mic',names(Par)))
  {
  res <- 1 - (x/Par$mic)^Par$alpha
  res <- res*(x<Par$mic)}else{
    res <- ((x-Par$xmin)/(Par$xopt-Par$xmin))*(x>Par$xmin)*(x<=Par$xopt)
  }
  return(res)
}
