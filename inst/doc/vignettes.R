## ----setup,include=FALSE------------------------------------------------------
library(SA23229012)

## -----------------------------------------------------------------------------
set.seed(100)
N<-500
p<-2
data_1 <- f_gener_data(N,p,"normal")
head(data_1)
data_2 <- f_gener_data(N,p,"normal",heteroskedasticity=FALSE)

## ----fig.width=10, fig.height=6, out.width="80%", fig.align='center',warning=FALSE,message=FALSE,collapse=TRUE----
library(aftgee)
library(dplyr)

fit1 <- aftgee(Surv(Y, status) ~ ., data = data_1)
beta_0 <-matrix(fit1$coefficients[,1])
print(computeCindex(as.numeric(beta_0), as.matrix(data_1)))

# draw pictures,demonstrate heterogeneity in variance 
data_main<-data_1
data_x1 <- data_main%>%select(-c("Y","status"))
data_x2_plot<- as.matrix(cbind(data.frame(rep(1,dim(data_main)[1])),data_x1))
mu <- data_x2_plot %*% beta_0
lnY_plot <- as.matrix(log(data_main["Y"]))
status<-data_main["status"]
plot(mu,lnY_plot,ann = F, xlab = expression("Mu (μ)"), ylab = expression("Log(Y)"))
points(mu[status == 0], lnY_plot[status == 0], col = "blue", pch = 16)  
points(mu[status == 1], lnY_plot[status == 1], col = "red", pch = 16)  
abline(a = 0, b = 1, col = "red")
mtext("The line represents the predicted mean, red=uncensoed, blue=censored", side = 1, line = 2, col = "black", cex = 1.2)

## ----include=FALSE------------------------------------------------------------
library(splines)
f_get_basisnumber <- function(data,beta)
{
  data_x1 <- data%>%select(-c("Y","status"))
  data_x2<- as.matrix(cbind(data.frame(rep(1,dim(data)[1])),data_x1))
  mu <- data_x2 %*% beta
  seq1 <-quantile(mu,c(0.5))
  s1<-as.matrix(bs(mu,knots=seq1,degree=2,intercept=TRUE)[,])
  beta_dim <- dim(s1)[2]
  return(beta_dim)
}
basis_number<-f_get_basisnumber(data_1,beta_0)

## -----------------------------------------------------------------------------
library(aftgee)
library(dplyr)

tau_0<-0.05
tau_list <- matrix(sort(seq(0,1-tau_0,tau_0), decreasing = TRUE),nrow = 1)
tau_pre_out<-matrix(c(0),nrow=1)
#basis_number means number of basis function
beta_pre_out<-matrix(c(rep(-10,basis_number)),nrow=basis_number)
#initialize
for (j in 1:(length(tau_list)-1)) 
    {
    beta_up<-f_update_basiscoef(tau=tau_list[1,length(tau_list)-j],
    data=data_1,beta=beta_0,
    tau_pre=tau_pre_out,beta_pre=beta_pre_out)
    
    beta_pre_out <- cbind(beta_up,beta_pre_out)
    
    tau_pre_out <- matrix(tau_list[1,(length(tau_list)-j):length(tau_list)],nrow=1)
    }

## ----include=FALSE------------------------------------------------------------
print(beta_pre_out)

## ----fig.width=10, fig.height=6, out.width="80%", fig.align='center',warning=FALSE,message=FALSE----
library(ggplot2)
data_x1 <- data_1%>%select(-c("Y","status"))
data_x2_plot<- as.matrix(cbind(data.frame(rep(1,dim(data_1)[1])),data_x1))
mu <- data_x2_plot %*% beta_0
seq1<- quantile(mu,c(0.5))
s1<-as.matrix(bs(mu,knots=seq1,degree=2,intercept=TRUE)[,])
mu_matrix <- cbind(mu,mu)
data_plot<- cbind(data.frame(mu),data.frame(mu_matrix+s1 %*% beta_pre_out[,c(5,15)]))
data_plot$lnY <- log(data_1$Y)
data_plot$status <- data_1$status

scatter_plot <- ggplot(data_plot, aes(x = mu, y = lnY, color =as.character(status))) +
    geom_point() +
    labs( x = expression(beta * x) , y = "ln T",color = "Status") +
    scale_color_manual(values = c("0" = "blue", "1" = "red"))  # 
  curve_plot <- scatter_plot +
    geom_smooth(aes(y = X1), color = "red") +
    geom_smooth(aes(y = X2), color = "green") +
    labs(title = "Scatter Plot with Quartiles Curves,Upper quartile and lower quartile",x = expression(beta * x), y = "Y")

print(curve_plot)

