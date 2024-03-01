#' Function which returns the acid concentration in millimolar

#' @param XX_list: Named data.frame with elements "pH", "DeltaL", 
#' "wHA" and "wSalt", each element containing one vector with the value 
#' to be evaluated.
#' @param Par_list: Named list with elements "pKa", "MHA", "MSalt"
#' @param pH: The pH
#' 
#' @return A scalar i.e. the acid concentration in millimolar 
#'
#' @author  Lars Erik Gangsei
#' 
#' @examples
#' # Load data
#' 
#' @export
#' 

Acid_conc_func <- function(XX_list,Par_list,pH=NULL)
{
  if(is.null(pH)){pH <- XX_list$pH}
  if(is.element('wHA',names(XX_list)))
  {if(!is.null(XX_list$wHA)){cp <- (XX_list$wHA*10^4)/Par_list$MHA}else{
    cp <- XX_list$wSalt/Par_list$MSalt}
    }else{cp <- XX_list$wSalt/Par_list$MSalt}
  
  return(cp/((10^(pH-Par_list$pKa)+1)*XX_list$DeltaL))
}