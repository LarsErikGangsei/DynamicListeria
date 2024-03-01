################################################################################
##                                                                            ##
## R-script reproducing results presented in manuscript:                      ##
##      "Extension of cardinal growth models for Listeria monocytogenes to    ##
##       handle dynamically changing environments during shelf-life".         ##
##                                                                            ##
## Author (R-script):  - Lars Erik Gangsei (Animalia, NMBU),                  ##
##                       email: lars.erik.gangsei@animalia.no                 ## 
##                       phone: +47 950 61 231                                ##
##                                                                            ##
## Adresses: Norwegian University of Life Sciences (NMBU),                    ##
##           Faculty of Chemistry, Biotechnology and Food Science,            ##
##           P.O. Box 5003,                                                   ##
##           NO-1432 Ås, Norway                                               ##
##                                                                            ##
##           Animalia AS,                                                     ##
##           P.O. Box 396 - Økern,                                            ##
##           NO-0513 Oslo, Norway                                             ##
##                                                                            ##
##                                                                            ##
##                                                                            ##
## Manuscript resubmitted to Food Control in 29th of February 2024            ##

## 1) Clean workspace, load necessary packages ----------------------------   ##

#  Clean workspace              
rm(list = ls())

# Load packages
packages <- c('xtable','devtools','tidyverse','plot3D')

if (length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages())))
}

for(pp in packages){library(pp,character.only=TRUE)}


# Assuming current directory is to "DynamicListeria"  
devtools::load_all()

## 2) Figures for illustrating effect of pH on Listeria growth -------------- ##
# Estimate beta2
InputpH <- dplyr::select(Listeria_pH_data,
                         Day,pH,Temp,Carbohydrates,Packing)%>%
           dplyr::rename(t = Day,carb = Carbohydrates)%>%
           dplyr::mutate(pH0 = mean(pH[t==0]),pHmin = 4.3)

beta2hat <- Estimation_func_beta2(t=InputpH$t,pH = InputpH$pH,Temp = InputpH$Temp,
                      carb = InputpH$carb,pH0= InputpH$pH0,
                      pHmin = InputpH$pHmin,Temppar = CardinalPar$Temperature)


# Estimate mu_opt
InputListeria <- dplyr::select(Listeria_pH_data,Day,Temp,cfu_Listeria,
                               Perc_Liq_phase,aW,Carbohydrates,pH,Packing)%>%
                 dplyr::rename(yt = cfu_Listeria,t = Day,DeltaL = Perc_Liq_phase,
                    Temperature = Temp,carb = Carbohydrates)%>%
                 dplyr::mutate(y0 = mean(yt[t==0]),Acetate = 0,Lactate = 0,
                           pH0 = mean(pH[t==0]),pHmin = 4.3)%>%
                 dplyr::select(-pH)

Parameter_Est <- Estimation_func(InputListeria = InputListeria,beta2 = beta2hat,
                                 pHmin = 4.3,Par_List = NULL,Interaction=TRUE)

# Figure showing observed and predicted Listeria growth
png(width = 4*1280, height = 4*720, units = 'px',res = 300, 
    filename = 'ListeriaGrowth_real_and_predicted.png')

par(mar = c(8, 8, 8, 2) + 0.1)
plot(InputListeria$t*(1+(77/21-1)*(InputListeria$Temperature==12)),
     log10(InputListeria$yt/InputListeria$y0),
     col = gsub('6','blue',gsub('12','red',as.character(InputListeria$Temperature))),
     axes = FALSE,pch = as.numeric(as.factor(InputListeria$Packing))+16,
     cex = 1.5,xlab ='',ylab ='')#,
     #main = 'Listeria growt vs. time')
box()

mtext('Time (days) for temperature 12 celsius',side = 3,line = 5,cex = 2)
mtext('Time (days) for temperature 6 celsius',side = 1,line = 5,cex = 2)
mtext('Relative Listeria density, log scale',side = 2,line = 5,cex = 2)

ShelfLife <- rep(NA,4)
names(ShelfLife) <- c('mu_opt_w_D6','mu_opt_w_ND6','mu_opt_w_D12','mu_opt_w_ND12')

# Function for plotting curves with predicted values
Curve_func <- function(x,tt,mm){
  n <- length(x)
  Input <- data.frame(t = (1/(1+(77/21-1)*(tt==12)))*x,
                      carb = rep(InputListeria$carb[1],n),
                      Temperature= rep(InputListeria$Temperature[1],n),
                      aW= rep(InputListeria$aW[1],n),
                      Lactate= rep(InputListeria$Lactate[1],n),
                      Acetate= rep(InputListeria$Acetate[1],n),
                      DeltaL= rep(InputListeria$DeltaL[1],n),
                      pH0= rep(InputListeria$pH0[1],n),
                      y0 = rep(1,n),
                      mu_opt=rep(Parameter_Est[[mm]],n))
  predict_Listeria(Input,
                   pHpar = list(beta2 = ifelse(mm =='mu_opt_w_ND',0,Parameter_Est$beta2),
                                pH0 = InputpH$pH0[1],pHmin = InputpH$pHmin[1],
                                carb = InputpH$carb[1],
                                temp = tt,Temppar = NULL),
                   Par_list = NULL,Interaction=TRUE,log10Return = TRUE)
}

for(tt in c(6,12)){
  for(mm in c('mu_opt_w_D','mu_opt_w_ND'))
  {
    curve(Curve_func(x,tt = tt,mm=mm),
      add = TRUE,col = ifelse(tt==6,'blue','red'),lwd = 2,
      lty = (2+(mm == 'mu_opt_w_ND')))

    ShelfLife[paste(mm,tt,sep='')] <- Shelflife_predict(Temperature = tt,
                    Acetate = InputListeria$Acetate[1],
                    Lactate = InputListeria$Lactate[1],aW = InputListeria$aW[1],
                    y0 = 1,mu_opt=Parameter_Est[[mm]],
                    DeltaL = InputListeria$DeltaL[1],
                    pHpar = list(beta2 = ifelse(mm =='mu_opt_w_ND',0,Parameter_Est$beta2),
                    pH0 = InputpH$pH0[1],pHmin = InputpH$pHmin[1],
                    carb = InputpH$carb[1],temp = tt,TempparpH = NULL),
                    Par_list = NULL,Interaction=TRUE)
  
  }
}
abline(h = 2,col = 'black',lty = 1,lwd = 2)
arrows(x0 = ShelfLife*c(1,1,77/21,77/21),y0 = 2,y1=c(-0.2,-0.2,6,6),
       col = rep(c('blue','red'),each = 2),lty = c(2,3,2,3),lwd =2)

axis(1,at = c(0,ShelfLife[1:2],25,50,75),cex.axis = 2,
     labels = round(c(0,ShelfLife[1:2],25,50,75),1),las= 2)
axis(2,at = seq(0,6,by = 1),cex.axis = 2,)
axis(3,at = c(0,ShelfLife[3:4],10,15,20)*77/21,cex.axis = 2,
     labels = round(c(0,ShelfLife[3:4],10,15,20),1),las= 2)

legend('bottomright',col = c('blue','red','gray','gray','gray','gray','gray'),
       legend = c('6 degrees celsius', '12 degrees celsius',
                  'Fluctating environment','Constant environment','Air','CO2','Nitrogen'),
       pch = c(15,15,NA,NA,17:19),lty = c(rep(NA,2),2,3,rep(NA,3)),cex = 2,
       lwd = 2)

dev.off()

# Figure showing observed pH and predicted pH
png(width = 4*1280, height = 4*720, units = 'px',res = 300, 
    filename = file.path(GitPath,'Figures/pH_real_and_predicted.png'))

par(mar = c(8, 8, 8, 2) + 0.1)
plot(InputpH$t*(1+(77/21-1)*(InputpH$Temp==12)),
     InputpH$pH,col = gsub('6','blue',gsub('12','red',as.character(InputpH$Temp))),
     axes = FALSE,pch = as.numeric(as.factor(InputpH$Packing))+16,
     cex = 1.5,xlab ='',ylab ='',ylim = c(4,8),xlim = c(0,80))

box()

mtext('Time (days) for temperature 12 celsius',side = 3,line = 5,cex = 2)
mtext('Time (days) for temperature 6 celsius',side = 1,line = 5,cex = 2)
mtext('pH',side = 2,line = 5,cex = 2)

curve(predict_pH(x,Temp=6,carb = InputpH$carb[1],pH0=InputpH$pH0[1],
                 pH_min=InputpH$pHmin[1],beta2=Parameter_Est$beta2,
                 Temppar = NULL), col ='blue',lty = 2,lwd = 2,add = TRUE)
curve(predict_pH((21/77)*x,Temp=12,carb = InputpH$carb[1],pH0=InputpH$pH0[1],
                 pH_min=InputpH$pHmin[1],beta2=Parameter_Est$beta2,
                 Temppar = NULL), col ='red',lty = 2,lwd = 2,add = TRUE)

axis(1,at = c(0,ShelfLife[1:2],25,50,75),cex.axis = 2,
     labels = round(c(0,ShelfLife[1:2],25,50,75),1),las= 2)
axis(2,at = seq(4,8,by = 0.5),cex.axis = 2)
axis(3,at = c(0,ShelfLife[3:4],10,15,20)*77/21,cex.axis = 2,
     labels = round(c(0,ShelfLife[3:4],10,15,20),1),las= 2)

#abline(h = CardinalPar$pH$xmin,lwd = 2)

legend('topright',col = c('blue','red','gray','gray','gray'),
       legend = c('6 degrees celsius', '12 degrees celsius',
                  'Air','CO2','Nitrogen'),pch = c(15,15,17:19),cex = 2)
dev.off()

## 3) Simulation study ------------------------------------------------------ ##
Simpar <- list(mu_opt = 6.31,sigma_mu = NA,y0 = 10^2,beta2 = NA,
               sigma_b = NA,carb = 11,pH0 = 7.51,pHmin = 4.62,
               TempPar = NULL)

Enviorment <- data.frame(t= c(0,5,20,50,100,0,3,15,30,75,0,1,3,8,20,0,1,2,5,10),
                         Temperature= rep(c(4,6,12,20),each = 5),
                         Acetate= rep(0,20),
                         Lactate= rep(0,20),
                         aW = rep(0.99,20),
                         yt = rep(NA,20),
                         y0 = rep(Simpar$y0,20),
                         pH0 = rep(Simpar$pH0,20),
                         pHmin = rep(Simpar$pHmin,20),
                         carb = rep(Simpar$carb,20),
                         DeltaL = rep(1,20))

VarComps <- data.frame(sigma_mu = c(0,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5)/10,
                       sigma_beta = c(0,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5),
                       beta2 = c(0,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1,0.5)*10)

rr <- 50
mm <- dim(VarComps)[1]
nn <- mm^3*rr


InputpH <- list(t = Enviorment$t, pH = rep(NA,dim(Enviorment)[1]),
                Temp = Enviorment$Temperature,
                carb = Enviorment$carb,
                aW = Enviorment$aW,
                pH0 = rep(Simpar$pH0,dim(Enviorment)[1]),
                pHmin = rep(Simpar$pHmin,dim(Enviorment)[1]))

Par_List_pH <- list(TempPar = NULL,beta2 = NULL)

set.seed(15)
Seed_vec <- sample(2*nn,nn,replace = FALSE)


SimRes <- data.frame(sigma_mu = rep(VarComps$sigma_mu,each = rr*mm^2),
                     sigma_beta = rep(rep(VarComps$sigma_beta,each = rr*mm),mm),
                     beta2 = rep(rep(VarComps$beta2,each = rr),mm^2),
                     beta_hat = rep(NA,nn),
                     mu_opt_w_D = rep(NA,nn),
                     mu_opt_f_D = rep(NA,nn),
                     mu_opt_w_ND = rep(NA,nn),
                     mu_opt_f_ND = rep(NA,nn))

seed_nr <- 1
for(ii in 1:nn)
  {
  Simpar$sigma_mu <- SimRes$sigma_mu[ii]
  Simpar$sigma_b <- SimRes$sigma_beta[ii]
  Simpar$beta2 <- SimRes$beta2[ii]
  set.seed(Seed_vec[seed_nr])
  
  Simdata <- Simulation_func(Enviorment,Simpar,Par_list=NULL,Interaction = TRUE)
  InputpH$pH <- Simdata$pH
  Enviorment$yt <- Simdata$yt
  
  beta2hat <- Estimation_func_beta2(t=Simdata$t,pH = Simdata$pH,Temp = Enviorment$Temperature,
                                    carb = Enviorment$carb,pH0= Enviorment$pH0,
                                    pHmin = Enviorment$pHmin,
                                    Temppar = CardinalPar$Temperature)
  
  Par_est <- Estimation_func(InputListeria = Enviorment,pHmin = Enviorment$pHmin[1],
                             beta2 = beta2hat,Par_List = NULL,
                             Interaction=TRUE)
  SimRes$beta_hat[ii] <- Par_est$beta2
  SimRes$mu_opt_w_D[ii] <- Par_est$mu_opt_w_D
  SimRes$mu_opt_f_D[ii] <- Par_est$mu_opt_f_D
  SimRes$mu_opt_w_ND[ii] <- Par_est$mu_opt_w_ND
  SimRes$mu_opt_f_ND[ii] <- Par_est$mu_opt_f_ND
  seed_nr <- seed_nr+1
  print(ii)
}



## 4) Evaluate Simulation Study --------------------------------------------- ##
Methods <- data.frame(PrintName = c('Dynamical environment, weighted estimate',
                                     'Dynamical environment, unweighted estimate',
                                     'Static environment, weighted estimate',
                                     'Static environment, unweighted estimate'),
                      DataName = c('mu_opt_w_D','mu_opt_f_D','mu_opt_w_ND',
                                   'mu_opt_f_ND'))

# Regression models
SimMods <- vector('list',4)
names(SimMods) <- Methods$DataName

for(ii in 1:dim(Methods)[1])
{
  dn <- Methods$DataName[ii]
  #print(dn)
  Innput_var <- dplyr::group_by(SimRes,sigma_mu,sigma_beta,beta2)%>%
    dplyr::summarise(Var_est = sqrt(var(.data[[dn]][abs(.data[[dn]])<Inf],na.rm=TRUE)))
  
    SimMods[[ii]] <- list(Var_mod = lm('Var_est~sigma_mu*sigma_beta*beta2',
                                     data = Innput_var),
                        Mean_mod = NULL,
                        Var_Summary = quantile(Innput_var$Var_est,probs = c(0.025,0.5,0.975)))
    
   SimMods[[ii]]$Mean_mod <- lm(paste(dn,'~beta2*sigma_mu*sigma_beta',sep = ''),
                               data = SimRes,
                               weights = 1/pmax(0.001,
                              predict(SimMods[[ii]]$Var_mod,newdata = SimRes)))
}

## Make table for variances in mu_opt with different methods.
Tab_Var <- round(sapply(SimMods,function(x) x[['Var_Summary']]),3)
colnames(Tab_Var) <- Methods$PrintName

xtable(Tab_Var)

## Regression tables

Reg_tab <- rbind(sapply(lapply(SimMods,function(x) x[['Var_mod']]),
       function(x) c(sum(anova(x)[,'Sum Sq'],na.rm = TRUE),
                     summary(x)$r.squared,
                     summary(x)$fstatistic[1])),
       sapply(lapply(SimMods,function(x) x[['Mean_mod']]),
              function(x) c(sum(anova(x)[,'Sum Sq'],na.rm = TRUE),
                            summary(x)$r.squared,
                            summary(x)$fstatistic[1])))

print(xtable(data.frame(c('','$\\widehat{\\mu_{opt}}$','','','$\\mu_{opt}$',''),
              rep(c('Sum of Squares','$R^2$','F-value'),2),
              Reg_tab)),include.rownames = FALSE,include.colnames = FALSE)


## 5) Make 3D figures ------------------------------------------------------- ##

sigma_mu_smooth <- VarComps$sigma_mu[1]
beta2_smooth <- VarComps$beta2[1]
for(ii in 2:mm)
{
  sigma_mu_smooth <- c(sigma_mu_smooth,seq(VarComps$sigma_mu[ii-1],
                                      VarComps$sigma_mu[ii],length.out = 11)[-1])
  beta2_smooth <- c(beta2_smooth,seq(VarComps$beta2[ii-1],
                                         VarComps$beta2[ii],length.out = 11)[-1])
}

df_predict <- expand.grid(sigma_mu_smooth,beta2_smooth)%>%
  dplyr::rename(sigma_mu = Var1,beta2=Var2)%>%
  dplyr::mutate(sigma_beta = 0.0001)
 
for(ii in 1:dim(Methods)[1])
{
  dn <- Methods$DataName[ii]
 # png(width = 4*680, height = 4*680, units = 'px',res = 300, 
 # filename = file.path(GitPath,'Figures',paste('Fig3D_',dn,'.png',sep='')))
  
  PmatMean <- matrix(predict(SimMods[[dn]]$Mean_mod,newdata = df_predict),
                        ncol = mm^2,nrow = mm^2)
  
  PmatSe <- matrix(qnorm(0.975)*sqrt(pmax(0,predict(SimMods[[dn]]$Var_mod,
            newdata = df_predict))),ncol = mm^2,nrow = mm^2)
  
  persp3D(x = (1:mm^2)/mm, y = (1:mm^2)/mm, 
          z = PmatMean, zlim = c(0,10),xlab = '',ylab = '', 
          zlab = '',col = 'blue',theta = 120,phi = 5,border = 'black',
          box = FALSE,alpha = 1,cex.main = c(3,3,0.1,0.1)[ii],
          main = c('Weighted estimate','Unweighted estimate','','')[ii])
  
  mtext(c('Dynamic E.  model','','Constant E. model','')[ii], 
        side = 2, line = 1, las = 0,cex = 3)
  
  # X-axis
  arrows3D(0,0,0,mm,0,0,colkey = FALSE,add=TRUE,col = 'black',
           lwd = 1.5)
  text3D(4,0,-1, expression(sigma[mu]^2),add=TRUE,cex.txt = 2.5)
  text3D(1:mm,rep(0,mm),rep(0.75,mm),
         as.character(VarComps$sigma_mu),cex.txt = 0.8,add=TRUE)
  
  # Y -axis
  arrows3D(mm,0,0,mm,mm,0,
           colkey = FALSE,add=TRUE,col = 'black',lwd = 1.5)
  text3D(mm,5,-1, expression(beta^2),add=TRUE,cex.txt = 2.5)
  text3D(rep(mm,mm),1:mm,rep(0.4,mm),as.character(VarComps$beta2),
         cex.txt = 0.8,add=TRUE)
 
  # Z - axis
  arrows3D(mm,0,0,mm,0,10,colkey = FALSE,add=TRUE,col = 'black',lwd = 1.5)
  text3D(rep(mm,10),rep(-0.75,10),1:10,as.character(1:10),cex.txt = 1.2,add=TRUE)
  text3D(mm,0.5,9, expression(mu[opt]), add = TRUE,cex.txt = 2.5)
  
  
  # Add surfaces with approximate confidence levels
  persp3D(x = (1:mm^2)/mm, y = (1:mm^2)/mm,z = PmatMean-PmatSe,
          add=TRUE,border = NA,facets = TRUE,alpha = 0.5,axes = FALSE,box = FALSE,
          colkey = FALSE,col = 'red')
  
  persp3D(x = (1:mm^2)/mm, y = (1:mm^2)/mm,z = PmatMean+PmatSe,
          add=TRUE,border = NA,facets = TRUE,alpha = 0.5,axes = FALSE,box = FALSE,
          colkey = FALSE,col = 'red')
  
  #dev.off()
  
}


