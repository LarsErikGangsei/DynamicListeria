#' Function calculating interaction between temperature, pH and wateractivity in cardinal modelling

#' @param XX_list: Named list (normally "Temp", "pH", "aW","Lactat" and "Acetat"), each element containing one scalar 
#' with the value to be evaluated.
#' @param Par_list: Named list with same names as XX_list (normally "Temp", "pH", "aW", '"Lactat" and "Acetat") 
#' each list  contains elements (at least) "xmin", "xmax" and "xopt" (all scalars)

#' @return A scalar >=0 and <=1. 
#'
#' @author  Lars Erik Gangsei
#' @references Augustin, J., Zuliani, V., Cornu, M., and Guillier, L. (2005). 
#' Growth rate and growth probability ofListeria monocytogenes in dairy, meat and seafood 
#' products in suboptimal conditions. Journal of Applied Microbiology, 99(5):1019-1042.
#' doi: https://doi.org/10.1111/j.1365-2672.2005.02710.x   
#' 
#' @examples
#' Load data
#' 
#' @export

Interaction_func <- function(XX_list,Par_list)
{
  
  Par_list <- Par_list[names(XX_list)]
  
  phi_vec <- c(mapply(function(x,y){((y$xopt-x)/(y$xopt-y$xmin))^3},
                      XX_list[c('Temperature','pH','aW')],
                      Par_list[c('Temperature','pH','aW')]),
               unlist((1-SR_func(XX_list['Lactate'],Par_list[[c('Lactate')]])*
                  SR_func(XX_list['Acetate'],Par_list[[c('Acetate')]]))^2))
  
  
  psi <- sum(phi_vec*(1-phi_vec)/(2*prod(1-phi_vec)))
  if(is.nan(psi)){psi<-1}
  
  res <- as.numeric(psi<=0.5) + as.numeric((psi>0.5)&(psi<1))*2*(1-psi)
  return(as.vector(res))
}