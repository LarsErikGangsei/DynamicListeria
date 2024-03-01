#' Functions for simulating pH and listeria concentrations 

#' @param Enviorment: Data frame of size n x 5 with elements: \cr
#' - t: Times (in days) for observations.\cr
#' - Temperature: Temperatures \cr
#' - Acetate: Added acetate in ?\cr
#' - Lactate: Added lactate in ?\cr
#' - aW: Wateractivity. \cr
#' @param Simpar: List with the following parameters\cr
#' - mu_opt: \eqn{\mu_{opt}}, see Gangsei et.a. for details.\cr
#' - sigma_mu:  \eqn{\sigma_{\mu}^2},see Gangsei et.a. for details.\cr 
#' - y0: \eqn{y_0} Listeria start consentration.
#' - beta2: \eqn{\beta^2}, see Gangsei et.a. for details.\cr
#' - sigma_b: \eqn{\sigma_{\beta}^2}, see Gangsei et.a. for details.\cr
#' - carb: \eqn{x_c}, carbohydrate level, see Gangsei et.a. for details.\cr
#' - pH0: \eqn{pH_0}, start pH (pH at time t = 0), , see Gangsei et.a. for details.\cr
#' - pHmin: \eqn{pH_{min}}, the minimum pH value (minimum cardinal pH value)
#' - TempPar: Cardinal parametersfor temperature for "pH -bacterias".
#' @param Par_list: Named list which is passed to "Cardinal_func". If set to NULL
#' the values in data "CardinalPar" is used.
#' @param Interaction: Passed directly to "Cardinal_func". 
#' If TRUE, interactions are calculated based on CM elements.
#' 
#' @return A data frame of size n x 5 containing elements\cr
#' - t: Times in days, directly from input "Enviorment"
#' - beta2: Simulated values for \eqn{\beta^2}
#' - pH: Simulated pH values\cr
#' - delta_mu: Random errors for instant listeria growth rate
#' - int_mu: \eqn{\int_0^t h(t) dt}, the overall growth potential without 
#' random error.\cr
#' - y_t: \eqn{y_t}, i.e simulated listeria concentrations at times defined by t.
#' 
#'
#' @author  Lars Erik Gangsei
#' 
#' @examples
#' # Load data
#' 
#' @export
#' 

Simulation_func <- function(Enviorment,Simpar,Par_list=NULL,
                            Interaction = TRUE)
{
  if(is.null(Par_list)){Par_list <- CardinalPar}
  if(is.null(Simpar$TempPar)){Simpar$TempPar <- CardinalPar$Temperature}
  
  ## Simulate pH
  nn <- dim(Enviorment)[1]
  
  pHpar <- lapply(split(data.frame(beta2 = pmax(0,rep(Simpar$beta2,nn)+
                          rnorm(nn,0,sqrt(Simpar$sigma_b))),
                      carb = rep(Simpar$carb,nn),
                      temp = Enviorment$Temperature,
                      pH0 = rep(Simpar$pH0,nn),
                      pHmin = rep(Simpar$pHmin,nn),
                      TempPar = rep(NA,nn)),1:nn),
                  function(x){x <- as.list(x);x$TempPar = Simpar$TempPar;return(x)})
  
  
  
  
  ## Get integrals from the simulated beta2
  res <- data.frame(t = Enviorment$t,
                    Beta2 = sapply(pHpar,function(x) x$beta2),
                    pH = mapply(pH_func,split(Enviorment$t,1:nn),pHpar),
                    delta_mu = rnorm(nn,0,sqrt(Simpar$sigma_mu)),
                    int_mu = rep(NA,nn),
                    yt = rep(NA,nn))
  
  Par_list_Exp <- vector('list',nn) 
  Par_list_Exp <- lapply(Par_list_Exp,function(x) x <- Par_list)
  
  res$int_mu <- mapply(Cardinal_func_integral,split(Enviorment,1:nn),
                      pHpar,Par_list_Exp,MoreArgs = list(Interaction = Interaction,
                      log=FALSE))
  
  res$yt <- Simpar$y0*10^(Simpar$mu_opt*res$int_mu+res$delta_mu*res$t)
  return(res)

  
}