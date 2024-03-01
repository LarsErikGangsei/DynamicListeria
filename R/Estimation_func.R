#' Functions for estimating pH and listeria concentrations 
#' @param InputListeria: data.frame (n x 11) with observations for Listeria containing at least:
#' the following elements\cr
#' - t: time (in days)\cr
#' - Temperature: in Celsius\cr
#' - Acetate: Acetate level\cr
#' - Lactate: lactate level\cr
#' - aW: Wateractivity\cr
#' - carb: \eqn{x_c}, carbohydrate level, see Gangsei et.a. for details.\cr
#' - DeltaL: Free water in the matrix.\cr
#' - yt: Listeria density at times t (cfu's)\cr
#' - y0: Listerialevels at time 0
#' - pH0: Start pH.\cr
#' @param beta2: Parameter giving rate of pH change.
#' @param pHmin: Minimum pH.
#' @param Par_list: Named list which is passed to "Cardinal_func". If set to NULL
#' the values in data "CardinalPar" is used.
#' @param Interaction: Passed directly to "Cardinal_func". 
#' If TRUE, interactions are calculated based on CM elements.
#' 
#' @return A list with 5 elements \cr
#' - beta2: Estimated value for \eqn{\beta^2}\cr
#' - mu_opt_w_D: Estimated value for \eqn{mu_{Opt}} based on weighted estimate and dynamic environment\cr
#' - mu_opt_f_D: Estimated value for \eqn{mu_{Opt}} based on unweighted estimate and dynamic environment\cr
#' - mu_opt_w_ND: Estimated value for \eqn{mu_{Opt}} based on weighted estimate and non-dynamic environment\cr
#' - mu_opt_f_ND: Estimated value for \eqn{mu_{Opt}} based on unweighted estimate and non-dynamic environment\cr
#' 
#'
#' @author  Lars Erik Gangsei
#' 
#' @import dplyr
#' 
#' @examples
#' 
#' @export
#' 

Estimation_func <- function(InputListeria,beta2 = NULL,pHmin = 4.3,
                            Par_List = NULL,Interaction=TRUE)
{
   # Basis for estimation of mu_opt based on dynamic enviorment
  nn <- dim(InputListeria)[1]

int_vec <-  sapply(split(dplyr::select(InputListeria,
                      t,carb,Temperature,aW,Lactate,Acetate,DeltaL,pH0),1:nn),
                      Cardinal_func_integral,
                      pHpar=list(beta2 = beta2,pHmin = pHmin),
                      Interaction = Interaction,log=FALSE)
  

ww <- (int_vec/InputListeria$t)^2
ww[abs(ww)==Inf] <- NA
ww <- ww/sum(ww,na.rm=TRUE)

# Basis for estimation of mu_opt based on non dynamic enviorment.
XX_list <- split(dplyr::rename(InputListeria[c('Temperature','Acetate',
                      'Lactate','aW','pH0','DeltaL')],pH = pH0),f = 1:nn)

int_vec_ND <-  sapply(XX_list,Cardinal_func,Par_list = Par_List,
                      Interaction = Interaction,log=FALSE)*InputListeria$t


ww_ND <- (int_vec_ND/InputListeria$t)^2
ww_ND[abs(ww_ND)==Inf] <- NA
ww_ND <- ww_ND/sum(ww_ND,na.rm=TRUE)


res <- list(beta2 = beta2,
            mu_opt_w_D = sum((ww*(log(InputListeria$yt)-log(InputListeria$y0))/
                             (log(10)*int_vec))[!is.na(ww)],na.rm=TRUE),
            mu_opt_f_D = mean(((log(InputListeria$yt)-log(InputListeria$y0))/
                              (log(10)*int_vec))[!is.na(ww)],na.rm = TRUE),
            mu_opt_w_ND = sum((ww_ND*(log(InputListeria$yt)-log(InputListeria$y0))/
                               (log(10)*int_vec_ND))[!is.na(ww_ND)],na.rm=TRUE),
            mu_opt_f_ND = mean(((log(InputListeria$yt)-log(InputListeria$y0))/
                                (log(10)*int_vec_ND))[!is.na(ww_ND)],na.rm = TRUE))


return(res)
}