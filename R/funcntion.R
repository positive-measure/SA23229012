
#' @import dplyr
#' @import splines
#' @import aftgee
#' @import survival
#' @import ggplot2
#' @import twosamples
#' @import coda
#' @import DAAG
#' @import microbenchmark
#' @import bootstrap
#' @import boot
#' @importFrom Rcpp evalCpp
#' @importFrom VGAM rgumbel
#' @importFrom quantreg rq
#' @importFrom stats rnorm rbinom rgamma runif coef quantile
#' @useDynLib SA23229012
NULL


#' @title Generate heteroscedastic survival data
#' @description Generate heteroscedastic survival data
#' @param N the number of samples
#' @param p dimension of samples
#' @param type The type of error distribution selected,"normal","extreme"
#' @param heteroskedasticity Whether the generated data has variance heterogeneity, defaults to TRUE
#' @return dataframe
#' @examples 
#' \dontrun{
#'         set.seed(100)
#'         N<-500
#'         p<-2
#'         data_1 <- f_gener_data(N,p,"extreme",heteroskedasticity=TRUE)
#' }
#' @export
f_gener_data<-function(N,p,type,heteroskedasticity=TRUE)
{
  f <- function(w,k) 
  {
    return((exp(k*w-exp(1)^w))/gamma(k))
  }
  
  if(p==2)
  {
    
    #data_1<-data.frame(rbinom(N,1,0.5),rbinom(N,1,0.5))
    data_1<-data.frame(rnorm(N,0,1),rbinom(N,1,0.5))
    list_col<-paste("X",1:p,sep="")
    
    ff <- function(x,a)
    {
      return(x-c(a))
    }
    
    
    colnames(data_1) <- list_col
    
    each_mean <-t(as.matrix(apply(data_1,2,mean)))
    mean_m <- matrix(rep(t(each_mean),N),ncol=ncol(each_mean),byrow=TRUE)
    data_1 <- data_1-mean_m
    #data_1 <- scale(data_1)
    
    data_1$Y <-exp(data_1["X1"]-data_1["X2"])
    colnames(data_1) <- c(list_col,"Y")
    status<-matrix(1,N,1)
    status<-data.frame(status)
  }
  
  if(type=="normal")
  {
    sample_y1<-function(x)
    {
      if(heteroskedasticity)
      {
        return(rnorm(1,mean=0,sd=x))
      }
      else{
        return(rnorm(1,mean=0,sd=1))
      }
    }
    
    Sigma_l<-abs(log(as.matrix(data_1["Y"])))+0.5
    
    eps <- apply(Sigma_l,1,sample_y1)
  }
  
  if(type=="extreme")
  {
    sample_y1<-function(x)
    {
      
      if(heteroskedasticity)
      {
      return(rgumbel(1, location = 0, scale =x)+digamma(1)* (x))
      }
      else
      {
        return(rgumbel(1, location = 0, scale =1)+digamma(1)* (1))
        
      }
    }
    
    
    Sigma_l<-((abs(log(as.matrix(data_1["Y"])))))^2+0.5
    #Sigma_l<-((abs(log(as.matrix(data_1["Y"])))))
    
    Sigma_l[Sigma_l<0.1]<-0.1
    eps <- apply(Sigma_l,1,sample_y1)
    
  }
  
  
  time<-data_1["Y"]*exp(eps)
  time[time>100]=100
  data_1["Y"]<-time
  data_1["status"]<-status
  status_n <- status
  censor_time <-runif(N,0,10)
  for (kk in 1:N) 
  {
    
    if(censor_time[kk]<time[kk,1])
    {
      time[kk,1]<-censor_time[kk]
      status_n[kk,1]<-0
    }
    
    
    
  }
  data_1["Y"]<-time
  data_1["status"]<-status_n
  
  
  return(data_1)
  
}


f_get_basisnumber <- function(data,beta)
{
  data_x1 <- data%>%select(-c("Y","status"))
  data_x2<- as.matrix(cbind(data.frame(rep(1,dim(data)[1])),data_x1))
  mu <- data_x2 %*% beta
  
  s1<-get_basisfunc(mu)
  beta_dim <- dim(s1)[2]
  return(beta_dim)
}

get_basisfunc <- function(mu0)
{
  seq1 <-quantile(mu0,c(0.5))
  s1<-as.matrix(bs(mu0,knots=seq1,degree=2,intercept=TRUE)[,])
  return(s1)
}


#' @title Obtain non-parametric estimates of tau quantiles
#' @description This function will return the coefficients under the spline basis function, and thus the user can obtain the tau quantile non-parametric estimator of the survival data by himself. And when using this function, you need to know those quantiles that are smaller than tau
#' @param tau The quantiles you want to estimate
#' @param data survival data
#' @param beta the estimator that has been obtained from AFT model
#' @param tau_pre Previous quantiles, vector
#' @param beta_pre The previous quantile estimator,Matrix,one column corresponds to a quantile
#' @return vector
#' @examples 
#' \dontrun{
#'      set.seed(100)
#'      N<-500
#'      p<-2
#'      data_1 <- f_gener_simple(N,p,"extreme")
#'      library(aftgee)
#'      library(survival)
#'      fit1 <- aftgee::aftgee(Surv(Y, status) ~ ., data = data_1)
#'      beta_0 <-matrix(fit1$coefficients[,1])
#'      tau_0<-0.04
#'      tau_list <- matrix(sort(seq(0,1-tau_0,tau_0), decreasing = TRUE),nrow = 1)
#'      basis_number<-f_get_basisnumber(data_main,beta_0)
#'      tau_pre_out<-matrix(c(0),nrow=1)
#'      beta_pre_out<-matrix(c(rep(-10,basis_number)),nrow=basis_number)
#'      #initial
#'      for (j in 1:(length(tau_list)-1)) 
#'      {
#'      beta_up<-f_update_basiscoef(tau=tau_list[1,length(tau_list)-j],
#'      data=data_1,beta=beta_0,
#'      tau_pre=tau_pre_out,beta_pre=beta_pre_out)
#'      beta_pre_out <- cbind(beta_up,beta_pre_out)
#'      tau_pre_out <- matrix(tau_list[1,(length(tau_list)-j):length(tau_list)],nrow=1)
#'      }
#' }
#' @export
f_update_basiscoef <- function(tau,data,beta,tau_pre,beta_pre)
{
  H <-function(x)
  {
    return(-log(1-x))
  }
  dim(data)[2]
  status1<-as.matrix(data["status"])
 
  
  data_x1 <- data%>%select(-c("Y","status"))
  data_x2<- as.matrix(cbind(data.frame(rep(1,dim(data)[1])),data_x1))
  mu <- data_x2 %*% beta
  lnY <- as.matrix(log(data["Y"]))

  s1<-get_basisfunc(mu)
  beta_dim <- dim(s1)[2]
  delta <- as.matrix(data["status"])
  time <- as.matrix(data["Y"])
  mu_matrix<-as.matrix(mu[,rep(1,each=length(tau_pre))])
  time_matrix<-as.matrix(time[,rep(1,each=length(tau_pre))])
  
  ts1 <- time_matrix-exp(mu_matrix+s1 %*% beta_pre)
  ts1[ts1<0]<-0
  ts1[ts1>0]<-1
  ts2<-colSums(ts1)
  tau_contain0<-as.matrix(tau_pre)
  ifelse(length(tau_pre)==1,{tau_containn<-matrix(c(tau),nrow=1)},
         {tau_containn<-matrix(c(tau,tau_pre[1,1:(length(tau_pre)-1)]),nrow=1)})
  delta_H<-apply(tau_containn, 2, H)-apply(tau_contain0, 2, H)
  delta_H<- matrix(delta_H,nrow = 1)
  
  
  
  censor_con<-2*ts1%*%t(delta_H)
  
  data_forquant <-s1
  y_forquant <-lnY-mu
  
  for (i in 1:dim(data_x1)[1])
  {
    if(status1[i,1]==0)
    {
      data_forquant[i,]<-censor_con[i,1]*s1[i,]
      #不同的数据集，取值不同,一般取个上界
      y_forquant[i,1]<-2
    }
  }
  
  
  data_qreg<-data.frame(cbind(data_forquant,y_forquant))
  
  colnames(data_qreg) <- c(paste("X",1:beta_dim,sep=""),"Y")
  
  model <- rq(Y ~ X1+X2+X3+X4-1., data = data_qreg, tau =tau )
  
  coefficients <- coef(model)
  

  coefficients_matrix <- matrix(coefficients, ncol = 1)
  
  return(coefficients_matrix)
}


