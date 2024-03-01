#' Function calculating CM value

#' @param x: vector of length n. The environmental factor (typically temperature, pH or water activity)
#' to be evaluated.
#' @param Par: List with the following elements "xmin", "xmax" ,"xopt" and "n" (all scalars)

#' @return a vector of length n with all elements in interval  >=0 and <=1. 
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
CM_func <- function(x,Par)
{
  nominator <- ((x-Par$xmin)^Par$n)*(x-Par$xmax)
  denominator <- ((Par$xopt-Par$xmin)^(Par$n-1))*((Par$xopt-Par$xmin)*(x-Par$xopt)-
                  (Par$xopt-Par$xmax)*((Par$n-1)*Par$xopt+Par$xmin-Par$n*x))
  res <- nominator/denominator
  res <- res*((x<Par$xmax)&(x>Par$xmin))
  return(res)
}
