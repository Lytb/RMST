################################################################################
###################################Load packages################################
################################################################################
    rm(list = ls())
    pacman::p_load(ggplot2, survRM2,survival,flexsurv,RISCA,survminer)

################################################################################
####################################Load data###################################
################################################################################
#Load Colon Dataset available in survival package.
    data_surv <- colon[colon$etype==2,c("rx","age","time","status")]

################################################################################
############################RMST WITHOUT FUNCTIONS##############################
#Fit models#####################################################################
#save the distributions that will be fitted as a vector            
    distributions <- c("exp","weibull", "gompertz","gamma","lnorm", "llogis")
    
    Models_st <- list(length(distributions))
        for(i in 1:length(distributions)){
            fitted_distrsp<- tryCatch({flexsurvreg(formula = Surv(time,status) ~ 1, data= subset(data_surv,rx=="Obs"), dist = paste(distributions[i]))},error=function(e){"NULL"})#use tryCatch to address the issue that sometimes models don't converge.
            Models_st[[i]]<-fitted_distrsp
    }
    

#ExtractData for Plot###########################################################
#https://stackoverflow.com/questions/52001230/r-flexsurv-estimating-a-survival-function
    pred <- setNames(data.frame(matrix(ncol = length(distributions)+1, nrow = length(seq(0, 3000, by=1)))), c("time", paste(distributions)))
    pred[,1] <-seq(0, 3000, by=1)#select the time for which to plot survival
        for (i in 1:length(distributions)) {
            pred[,i+1] <- summary(Models_st[[i]],t=seq(0, 3000, by=1))[[1]][,"est"]#Generate the new values for the novel timepoints
          }  
          
      

#Plots##########################################################################
    labels <- c("KM curve",paste(distributions))#labels for the plot
    r_col <- rainbow(length(distributions))#select the colours for the curves
    plot(survfit(Surv(time,status)~1, data = data_surv[data_surv$rx== "Obs",]), lty = c(1),xlab = "Time in Days", ylab = "Overall Survival", main= "Overall Survival for Observation Arm", xlim=c(0.01, 3000), conf.int = FALSE)#Plot the KM estimate
        for (i in 2:c(length(distributions)+1)){#plot the extrapolated survival curves for each distribution
            lines(x=pred[,1],y=pred[,i], col = r_col[i], lty = c(1))
    }
    legend("topright", legend = labels, col = c("black",r_col[1:length(distributions)]), lty = rep(1,length(distributions)+1), cex = 0.5, bty="n")#add the legend

#rmst parametric################################################################
    set.seed(1) # to ensure comparability between what is done here and when using a function to estimate these results
    rmst_all <- data.frame(matrix(nrow = 1, ncol = 5)) #create data.frame to store results and name columns
    colnames(rmst_all) <- c("model","time","est","lci","uci") 
    
        for (i in 1:length(distributions)){
            rmst_ <- c(c(paste(distributions[i])),round(unlist(print(summary(Models_st[[i]],type = "rmst", t=2190)[[1]])),4)) #create object that includes name of the model fitted and the rmst (incl lower and upperci)
            rmst_all[i,] <- rmst_
    }

#RMST KM survRM2################################################################
#https://cran.r-project.org/web/packages/survRM2/vignettes/survRM2-vignette3-2.html
#data_surv_sub <- data_surv[data_surv$rx=="Obs",]
    data_surv_sub <- data_surv
    levels(data_surv_sub$rx) <- c(0,1,2)
    set.seed(1)
    obj <- rmst2(time=data_surv_sub$time,status = data_surv_sub$status, arm=data_surv_sub$rx,tau=2190)
    rmst_all[length(distributions)+1,] <- c("SurvRM2",2190,round(obj[[4]]$rmst[1],4),round(obj[[4]]$rmst[3],4),round(obj[[4]]$rmst[4],4))

#RMST KM RISCA##################################################################
#Identical to SurvRM2 but doesn't report Confidence intervals
    km_colon <- summary(survfit(Surv(time,status) ~ 1, data= data_surv[data_surv$rx== "Obs",]))
    rmst_all[c(length(distributions)+2),] <- c("Risca rmst",2190,round(RISCA::rmst(km_colon$time,km_colon$surv,max.time = 2190),4),"","")


####################################Functions###################################
#Function.Model Fit#############################################################
    f_fittedmodels<-   function(distr,data_colon, treatment){
  
    ModelOutput <- list(length(distr))
        for(i in 1:length(distr)){
            fitted_distr<- tryCatch({flexsurvreg(formula = Surv(time,status) ~ 1, data= data_colon[data_colon$rx ==treatment,], dist = paste(distr[i]))}, error= function(e){"NULL"})#fit for each of the distributions
            ModelOutput[[i]] <- fitted_distr #save the distributions
        }
        return(ModelOutput)
  
  }



#Function.Extract Data for Plot#################################################
    f_pred <- function(distr, models_fitted){
    pred[,1] <-seq(0, 3000, by=1)
  
        for (i in 1:length(distr)) {
            pred[,i+1] <- summary(models_fitted[[i]],t=seq(0, 3000, by=1))[[1]][,"est"]#Generate the new values for the novel timepoints
        }  
        return(pred)
  }



#Function.Plot##################################################################
    f_plot <- function(distr,data_colon, treatment, title,pred){

        labels <- c("KM curve",paste(distr))
        r_col <- rainbow(length(distr))
        p1 <- plot(survfit(Surv(time,status)~1, data = data_colon[data_colon$rx==treatment,]), lty = c(1),xlab = "Time in Days", ylab = "Overall Survival", main= title, xlim=c(0.01, 3000), conf.int = FALSE)
            for (i in 2:c(length(distr)+1)){
            lines(x=pred[,1],y=pred[,i], col = r_col[i], lty = c(1))
            }
        legend("topright", legend = labels, col = c("black",r_col[1:length(distr)]), lty = rep(1,length(distr)+1), cex = 0.5, bty="n")
      
        return(p1)

    }



#Function.RMST##################################################################
    f_rmst <- function(distr,fittedm,rmst_time){
    set.seed(1)
    rmst_all <- data.frame(matrix(nrow = length(distr), ncol = 5))
    colnames(rmst_all) <- c("model","time","est","lci","uci") #Add column name
        for (i in 1:length(fittedm)){
            rmst_all[i,] <- c(paste(distributions[i]),tryCatch({round(unlist(print(summary(fittedm[[i]],type = "rmst", t=c(rmst_time))[[1]])),3)
    }, error=function(e) {"NULL"}))
        }
        return(rmst_all)
    }
    
##############################RESULTS USING FUNCTIONS###########################
#Standard parametric distributions in flexsurv##################################
            
    distributions <- c("exp","weibull", "gompertz","gamma","lnorm", "llogis")
    Models_st_func <- f_fittedmodels(distributions,data_surv, "Obs")
    pred_func <- f_pred(distributions,Models_st_func)
    plot_func <- f_plot(distributions, data_surv, "Obs","Overall Survival for Observation Arm",pred_func)
    rmst_all_func <- f_rmst(distributions,Models_st_func,2190)


################################ADDITIONAL INFO#################################
#Restricted vs mean survival time###############################################
#https://stackoverflow.com/questions/43173044/how-to-compute-the-mean-survival-time


#Other packages#################################################################
#RMST using Pseudovalues########################################################
#https://stat.ethz.ch/R-manual/R-devel/library/survival/html/pseudo.html
  #Non Parametric

#RMST:Restricted Mean Survival Time Methods Jonathan Wessen
#https://rdrr.io/github/scientific-computing-solutions/RMST/f/inst/doc/vignette.pdf