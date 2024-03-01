#' Function for estimating Shelf Life of product

#' @param Temperature: in Celsius
#' @param Acetate: Acetatelevel
#' @param Lactate: lactatelevel
#' @param aW: Wateractivity
#' @param y0: Listeria levels at time 0
#' @param mu_opt: value for \eqn{\mu_{opt}}
#' @param pHpar: List with the following elements
#' - carb: carbohydrate level\cr
#' - pH0: \eqn{pH_0}, start pH (pH at time t = 0), , see Gangsei et.a. for details.\cr
#' - pHmin: \eqn{pH_{min}}, the minimum pH value (minimum cardinal pH value)
#' - beta2: value for \eqn{\beta^2}
#' - TempparpH: Cardinal parameters for temperature for "pH -bacterias".
#' @param Par_list: Named list which is passed to "Cardinal_func". If set to NULL
#' the values in data "CardinalPar" is used.
#' @param Interaction: Passed directly to "Cardinal_func". 
#' @param log10Return: TRUE (default) or FALSE if results are to be returned at log10 scale.
#' If TRUE, interactions are calculated based on CM elements.
#' 
#' @return Predicted Shelf life
#'
#' @author  Lars Erik Gangsei
#' 
#' @examples
#' 
#' @export
#' 

Shelflife_predict <- function(Temperature,Acetate = 0,Lactate = 0,aW = 0.99,
                             y0 = 1,y_max = 100,DeltaL=1,mu_opt,pHpar,
                             Par_list = NULL,Interaction=TRUE)
{
  min_func <- function(x,Temperature,Acetate,Lactate,aW,y0,y_max,mu_opt,pHpar,
                       Par_List,Interaction)
  { n <- length(x)
  Input <- data.frame(t = x,
                      carb = rep(InputListeria$carb[1],n),
                      Temperature= rep(pHpar$carb[1],n),
                      aW= rep(aW,n),
                      Lactate= rep(Lactate,n),
                      Acetate= rep(Acetate,n),
                      DeltaL= rep(DeltaL,n),
                      pH0= rep(pHpar$pH0,n),
                      y0 = rep(y0,n),
                      mu_opt=rep(mu_opt,n))
    return((predict_Listeria(Input,
                            pHpar = pHpar,
                            Par_list = Par_list,Interaction=Interaction,
                            log10Return = TRUE)-(log10(y_max)-log10(y0)))^2)}
  
  res <- nlm(min_func,p = 1,Temperature = Temperature,Acetate = Acetate,
             Lactate = Lactate,aW = aW,y0 = y0,y_max = y_max,mu_opt=mu_opt,
             pHpar = pHpar,Par_List = Par_List,Interaction=Interaction)
  
  if(res$code==1){return(res$estimate)}else{return(Inf)}
}