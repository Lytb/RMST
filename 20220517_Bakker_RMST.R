################################################################################
###################################Load packages################################
################################################################################
    rm(list = ls())
    pacman::p_load(survRM2,survival,flexsurv,RISCA)

################################################################################
####################################Load data###################################
################################################################################
#Load Colon Dataset available in survival package.
    data_surv <- colon[colon$etype==2 & c(colon$rx== "Obs"|colon$rx== "Lev"),c("rx","age","time","status")]



#Function.Model Fit#############################################################
 
f_fittedmodels<-   function(distr,data_colon, treatment){
                   surv_model <- function(distr){tryCatch({flexsurvreg(formula = Surv(time,status) ~ 1, data= data_colon[data_colon$rx ==treatment,], dist = paste(distr))}, error= function(e){"NULL"})
                        }
                        
                   ModelOutput <- lapply(distr,surv_model)
                   return(ModelOutput)
      }
     
    
distributions <- c("exp","weibull", "gompertz","gamma","lnorm", "llogis")
Models_st_func <- f_fittedmodels(distributions,data_surv, "Obs")

#Function.Extract Data for Plot#################################################

f_pred <- function(distr, models_fitted){
          prediction <- function(models_fitted){summary(models_fitted,t=seq(0, 5000, by=1))[[1]][,"est"] }                 
          pred <- lapply(models_fitted,prediction)
          return(pred)
    }
    
pred_func <- f_pred(distributions,Models_st_func)

#plotCIs for Plot
f_pred_lci <- function(distr, models_fitted){
  
          prediction_lci <- function(models_fitted)summary(models_fitted,t=seq(0, 5000, by=1))[[1]][,"lcl"]
          pred_lci <- lapply(models_fitted,prediction_lci)
          return(pred_lci)
}

pred_func_lci <- f_pred_lci(distributions,Models_st_func)


f_pred_uci <- function(distr, models_fitted){
  
          prediction_uci <- function(models_fitted)summary(models_fitted,t=seq(0, 5000, by=1))[[1]][,"ucl"]
          pred_uci <- lapply(models_fitted,prediction_uci)
          return(pred_uci)
}

pred_func_uci <- f_pred_uci(distributions,Models_st_func)

#Function.Plot##################################################################

f_plot <- function(distr,data_colon, treatment, title,pred){

          labels <- c("KM curve",paste(distr))
          r_col <- rainbow(length(distr))
          p1 <- plot(survfit(Surv(time,status)~1, data = data_colon[data_colon$rx==treatment,]), lty = c(1), xlab = "Time in Days", ylab = "Overall Survival", main= title, xlim=c(0.01, 5000), conf.int = F)
                    for (i in 1:c(length(distr))){
                    lines(x=seq(0, 5000, by=1),y=pred[[i]], col = r_col[i], lty = c(1))
                    }
          legend("topright", legend = labels, col = c("black",r_col[1:length(distr)]), lty = rep(1,length(distr)+1), cex = 0.5, bty="n")
      
          return(p1)

      }

plot_func <- f_plot(distributions, data_surv, "Obs","Overall Survival for Observation Arm",pred_func)

#Plot with CIs
f_plot_withCI <- function(distr,data_colon, treatment, title,pred,pred_lci,pred_uci){
  
            labels <- c("KM curve",paste(distr))
            r_col <- rainbow(length(distr))
            p1 <- plot(survfit(Surv(time,status)~1, data = data_colon[data_colon$rx==treatment,]), lty = c(1),xlab = "Time in Days", ylab = "Overall Survival", main= title, xlim=c(0.01, 5000), conf.int = T)
            for (i in 1:c(length(distr))){
              lines(x=seq(0, 5000, by=1),y=pred[[i]], col = r_col[i], lty = c(1))
            }
            for (i in 1:c(length(distr))){
              lines(x=seq(0, 5000, by=1),y=pred_lci[[i]], col = r_col[i], lty = c(2))
            }
            for (i in 1:c(length(distr))){
              lines(x=seq(0, 5000, by=1),y=pred_uci[[i]], col = r_col[i], lty = c(2))
            }
            legend("topright", legend = labels, col = c("black",r_col[1:length(distr)]), lty = rep(1,length(distr)+1), cex = 0.5, bty="n")
            
            return(p1)
  
}

plot_func_withCI <- f_plot_withCI(distributions, data_surv, "Obs","Overall Survival for Observation Arm",pred_func, pred_func_lci,pred_func_uci)

#Function.RMST##################################################################
f_rmst <- function(fittedm,rmst_time){
          set.seed(1)
          f_est_rmst <- function(fittedm){tryCatch({round(unlist(print(summary(fittedm,type = "rmst", t=c(rmst_time))[[1]])),3)}, error=function(e) {"NULL"})}
              
          rmst_all <- lapply(fittedm, f_est_rmst)
          return(rmst_all)
    }
    
rmst_all_func <- f_rmst(Models_st_func,2190)
attributes(rmst_all_func) <- list(names= paste(distributions),row.names=c("time","est","lci","uci"),class="data.frame")
formattable(rmst_all_func)


#RMST KM survRM2################################################################
#https://cran.r-project.org/web/packages/survRM2/vignettes/survRM2-vignette3-2.html

levels(data_surv$rx) <- c("0","1","3")
data_surv$rx <- as.numeric(as.character(data_surv$rx))

set.seed(1)
obj <- rmst2(time=data_surv$time,status = data_surv$status, arm=data_surv$rx,tau=2190)
SurvRM2_rmst <- c(round(obj[[4]]$rmst[1],4),round(obj[[4]]$rmst[3],4),round(obj[[4]]$rmst[4],4))

#RMST KM RISCA##################################################################
#Identical to SurvRM2 but doesn't report Confidence intervals

km_colon <- summary(survfit(Surv(time,status) ~ 1, data= data_surv[data_surv$rx== 0,]))
risca_rmst <- round(RISCA::rmst(km_colon$time,km_colon$surv,max.time = 2190),4)

#summary(survfit(Surv(time,status) ~ 1, data= data_surv[data_surv$rx== 0,]),type = "rmst", t=c(2190))
  

################################ADDITIONAL INFO#################################
#Restricted vs mean survival time###############################################
#https://stackoverflow.com/questions/43173044/how-to-compute-the-mean-survival-time


#Other packages#################################################################
#RMST using Pseudovalues########################################################
#https://stat.ethz.ch/R-manual/R-devel/library/survival/html/pseudo.html
  #Non Parametric

#RMST:Restricted Mean Survival Time Methods Jonathan Wessen
#https://rdrr.io/github/scientific-computing-solutions/RMST/f/inst/doc/vignette.pdf
    
#Extrapolating values in plots
#https://stackoverflow.com/questions/52001230/r-flexsurv-estimating-a-survival-function
    