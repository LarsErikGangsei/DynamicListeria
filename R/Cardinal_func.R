#' Function calculating gamma value

#' @param XX_list: List with named elements with values to be evaluated
#' @param Par_list: List with elements (known parameters), same names as XX_list.\cr 
#' -If CM function is to be used names are: xmin,xmax,xopt and n (all scalars)\cr
#' -If SR function is to be used names are: alpha, mic\cr
#'  If Par_list is set to NULL the values in data "CardinalPar" is used.
#' @param Interaction: - If TRUE, interactions are calculated based on CM elements \cr
#' @param log: If the cardinal value is to be returned on log scale
#' @param Acid_Param: List with parameters for acids
#' @param Acid_Weights: Logical, if TRUE (default) acids are to be recalculated to 
#' consentrations

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
Cardinal_func <- function(XX_list,Par_list = NULL,Interaction = TRUE,log=FALSE,
                          Acid_Param = NULL,Acid_Weights = TRUE)
{
  if(is.null(Par_list)){Par_list <- CardinalPar}
  CM_names <- names(Par_list)[sapply(Par_list,function(x) 
    sum(is.element(names(x),c('xmin','xopt','xmax')))>0)]
  SR_names <- names(Par_list)[sapply(Par_list,function(x) 
    sum(is.element(names(x),c('alpha','mic')))>0)]
  
  if((Acid_Weights==TRUE)&&(sum(unlist(XX_list[SR_names]))>0))
  {
    if(is.null(Acid_Param)){Acid_Param <- Par_list_Acid}
    for(ss in SR_names)
    {
      XX_list[[ss]] <- Acid_conc_func(XX_list = list(pH = XX_list$pH,
                                                   wSalt = XX_list[[ss]],
                                                   wHA = NULL,
                                                   DeltaL = XX_list$DeltaL),
                                    Par_list= Acid_Param[[ss]])
    }
  }
  
  CM_val <- mapply(FUN=CM_func,XX_list[CM_names],Par_list[CM_names])
  SR_val <- mapply(FUN=SR_func,XX_list[SR_names],Par_list[SR_names])
  
  if(Interaction==TRUE)
  {Interact_val <- Interaction_func(XX_list,Par_list)
  }else{Interact_val <- 1}
  if(log==FALSE)
  {res <- prod(c(CM_val,SR_val,Interact_val))}else{
    res <- sum(log(c(CM_val,SR_val,Interact_val)))}
  return(res)
}
