## ----echo=TRUE,eval=TRUE------------------------------------------------------
my_sample<-function(data,prob)
{
  n<-dim(data)[1]
  
  if(sum(prob)!=1)
      {
    return("输入抽样概率错误！")
      }
    
  if(length(prob)!=n)
      {
    return("输入抽样概率错误！")
      }
    sample_prob<-prob
 
    
  cum_prob<-cumsum(sample_prob)
  
  data_s <-data
  

  for (i in 1:n) 
  {
    rd <-runif(1,0,1)
    index<-findInterval(rd,cum_prob)
    data_s[i,] <- data[i,]
  }
  
  return(data_s)
}


## ----echo=TRUE,eval=TRUE,collapse=TRUE----------------------------------------

data1<- data.frame(c(1,2,3),c(2,3,4),c(3,4,5)) 
colnames(data1)<-c("X1","X2","X3")
data_sample1 <-my_sample(data1,c(1/3,1/3,1/3))
sam_prob <-c(seq(1:dim(data1)[1]))/(3*(3+1)/2) #构造一个抽样概率
data_sample2 <-my_sample(data1,sam_prob)
data_sample3 <-my_sample(data1,c(1,2)) #函数会抛出错误

print(data_sample1)
print(data_sample3)

## -----------------------------------------------------------------------------
# 定义生成标准Laplace分布随机样本的函数
mysample2 <- function(n) {
  # 生成n个均匀分布在[0, 1]之间的随机数
  u <- runif(n)
  
  # 使用逆变换方法生成标准Laplace分布的随机样本
  laplace_sample <- -sign(u - 0.5) * log(1 - 2 * abs(u - 0.5))
  
  return(laplace_sample)
}

# 生成一个大小为1000的标准Laplace分布随机样本
sample_size <- 1000
laplace_sample <- mysample2(sample_size)

head(laplace_sample)

## ----echo=TRUE,eval=TRUE------------------------------------------------------

# 定义生成Beta分布随机样本的函数
generate_beta_sample <- function(n, a, b) {
  # 定义接受-拒绝采样函数
  accept_reject <- function() {
    while (TRUE) {
      # 生成两个均匀分布随机数 U1 和 U2
      u1 <- runif(1)
      u2 <- runif(1)
      
      # 计算接受概率
    f_u1 <- dbeta(u1, a, b)  # Beta(a, b)的概率密度函数
      
    #这是密度函数的最大值
      g_u1 <- 12
   
      if (u2 <= f_u1 / g_u1) {
        return(u1)
      }
    }
  }
  
  # 生成随机样本
  beta_sample <- numeric(n)
  for (i in 1:n) {
    beta_sample[i] <- accept_reject()
  }
  
  return(beta_sample)
}

# 生成Beta(3, 2)分布的随机样本
sample_size <- 1000
a <- 3
b <- 2
beta_sample <- generate_beta_sample(sample_size, a, b)


hist(beta_sample, probability = TRUE, breaks = 30, ylim = c(0, 2.5), main = "Histogram of Beta(3, 2) Sample")

# 叠加理论Beta(3, 2)密度曲线
x <- seq(0, 1, length.out = 1000)
y <- dbeta(x, a, b)
lines(x, y, col = "red", lwd = 2)

# 添加图例
legend("topright", legend = "Beta(3, 2) Density", col = "red", lwd = 2)

## -----------------------------------------------------------------------------
# 定义生成随机变量样本的函数
generate_random_sample <- function(n) {
  random_sample <- numeric(n)
  
  for (i in 1:n) {
    # 生成Uniform(-1, 1)分布的三个随机数
    u1 <- runif(1, min = -1, max = 1)
    u2 <- runif(1, min = -1, max = 1)
    u3 <- runif(1, min = -1, max = 1)
    
    # 判断 
    if (abs(u3) >= abs(u2) && abs(u3) >= abs(u1)) {
      random_sample[i] <- u2
    } else {
      random_sample[i] <- u3
    }
  }
  
  return(random_sample)
}

# 生成模拟随机样本
sample_size <- 10000
random_sample <- generate_random_sample(sample_size)

# 绘制直方图密度估计
hist(random_sample, probability = TRUE,ylim = c(0, 1), breaks = 15)

# 添加理论rescaled Epanechnikov核密度曲线
x <- seq(-1, 1, length.out = 1000)
epanechnikov_density <- (3/4) * (1 - x^2)
lines(x, epanechnikov_density, col = "red", lwd = 2)

# 添加图例
legend("topright", legend = "Rescaled Epanechnikov Kernel Density", col = "red", lwd = 2)

## ----collapse=TRUE------------------------------------------------------------
# 生成从卡方分布抽样的随机数
n=20
S<-0
M<-10000
alpha=0.05
for (i in 1:M) 
  {
  
  chisq_sample <- rchisq(n, df = 2)
  t_alpha_over_1 <-mean(chisq_sample)+ qt(alpha / 2, df = n-1)*var(chisq_sample)/sqrt(n)
  t_alpha_over_2 <-mean(chisq_sample)+ qt(1 - alpha / 2, df = n-1)*var(chisq_sample)/sqrt(n)
  
  if(t_alpha_over_1<2 & t_alpha_over_2>2)
  {
    S<-S+1
  }
  
}
print(S/M)

## ----collapse=TRUE------------------------------------------------------------
set.seed(123)  

# 模拟参数
n_simulations <- 10000  
sample_size <- 30  
alpha <- 0.05 


simulate_t_test <- function(distribution, mu0) {
  rejections <- 0
  
  for (i in 1:n_simulations) {
   
    if (distribution == "chi-square") {
      data <- rchisq(sample_size, df = 1)
    } else if (distribution == "uniform") {
      data <- runif(sample_size, min = 0, max = 2)
    } else if (distribution == "exponential") {
      data <- rexp(sample_size, rate = 1)
    }
    

    t_result <- t.test(data, mu = mu0)
    

    p_value <- t_result$p.value
    if (p_value < alpha) {
      rejections <- rejections + 1
    }
  }
  

  return(rejections / n_simulations)
}

# 三种不同分布的模拟试验
distribution_names <- c("chi-square", "uniform", "exponential")
mu0_values <- c(1, 1, 1)  # 对应三种分布的均值

for (i in 1:length(distribution_names)) {
  distribution_name <- distribution_names[i]
  mu0 <- mu0_values[i]
  
  error_rate <- simulate_t_test(distribution_name, mu0)
  
  cat(paste("Distribution:", distribution_name, "\n"))
  cat(paste("True Mean (mu0):", mu0, "\n"))
  cat(paste("Simulated Empirical Type I Error Rate:", error_rate, "\n\n"))
}


## ----collapse=TRUE------------------------------------------------------------

alpha<-0.1
m<-1000

#初始化
count_Bf<-numeric(3)
count_BH<-numeric(3)

for (i in 1:m) {
  

p1<-runif(950,0,1)
p2<-rbeta(50,0.1,1)
p_vec<-c(p1,p2)

p_Bf<-p.adjust(p_vec,method ='bonferroni')
p_BH<-p.adjust(p_vec,method ='BH')

#计算FWER
if(min(p_Bf[1:950])<alpha)
{
  count_Bf[1]<-count_Bf[1]+1
}

if(min(p_BH[1:950])<alpha)
{
  count_BH[1]<-count_BH[1]+1
}


#计算FDR，先初始化错误率
error_rate_Bf<-0
error_rate_BH<-0


  error_rate_Bf<-sum(p_Bf[1:950]<alpha)/sum(p_Bf[1:1000]<alpha)
  error_rate_BH<-sum(p_BH[1:950]<alpha)/sum(p_BH[1:1000]<alpha)

  count_Bf[2]<-error_rate_Bf+count_Bf[2]
  count_BH[2]<-error_rate_BH+count_BH[2]
  
  
  
  #计算TPR

  
  count_Bf[3]<-sum(p_Bf[951:1000]<alpha)/50+count_Bf[3]
  count_BH[3]<-sum(p_BH[951:1000]<alpha)/50+count_BH[3]
  
  }

print("两种方法的FWER,FER,TPR依次为:")
print(count_Bf/m)
print(count_BH/m)


## ----collapse=TRUE------------------------------------------------------------

lambda <- 2
num_bootstrap_samples <- 1000
M<-1000


for (id in c(5,10,20)) {
  
sample_size <- id

# 生成原始指数分布的数据


# 定义函数执行Bootstrap抽样并计算统计量
bootstrap_mean <- function(data, num_samples) 
{
  bootstrap_means <- numeric(num_samples)
  for (s in 1:num_samples) {
    resampled_data <- sample(data, replace = TRUE)
    bootstrap_means[s] <- 1/mean(resampled_data)
  }
  return(mean(bootstrap_means))
}

bootstrap_var <- function(data, num_samples) {
  bootstrap_means <- numeric(num_samples)
  for (s in 1:num_samples) {
    resampled_data <- sample(data, replace = TRUE)
bootstrap_means[s] <- 1/mean(resampled_data)
  }
  bootstrap_means<-bootstrap_means-lambda
  return(var(bootstrap_means))
}


# 执行Bootstrap抽样
bootstrap_bias<-numeric(M)
bootstrap_vars<-numeric(M)
for (i in 1:M) {
  
  original_data <- rexp(sample_size, rate = lambda)
  
bootstrap_bias[i] <- bootstrap_mean(original_data, num_bootstrap_samples)-lambda
bootstrap_vars[i] <- bootstrap_var(original_data, num_bootstrap_samples)

}
# 打印结果
cat("样本大小:",sample_size,"\n")
cat("Bootstrap bias:",mean(bootstrap_bias) , "bias真值", 2/(sample_size-1), "\n")
cat("Bootstrap方差:", mean(bootstrap_vars), "真实方差",lambda*sample_size/(sample_size*sqrt(sample_size-2)), "\n")

}


## ----collapse=TRUE------------------------------------------------------------
library("bootstrap")
law<-law82[c("LSAT","GPA")]
M<-1000
C<-numeric(M)
for(i in 1:M)
{
  sample_data<-sample(law,replace=TRUE)
 
  C[i]<- cor(sample_data[,1],sample_data[,2])
}

print(mean(C))
cat("置信区间为:(",mean(C)-1.96*sd(C),",",mean(C)+1.96*sd(C),")")



## ----collapse=TRUE------------------------------------------------------------
# 原始数据

data <- c(3,5,7,18,43,85,91,98,100,130,230,487)

# 设置Bootstrap参数
n <- length(data)
num_bootstrap_samples <- 1000
alpha <- 0.05  # 95% confidence interval

# 初始化变量来存储Bootstrap样本均值
bootstrap_sample_means <- numeric(num_bootstrap_samples)

# 执行Bootstrap抽样
for (i in 1:num_bootstrap_samples) {
  resampled_data <- sample(data, replace = TRUE)
  bootstrap_sample_means[i] <- mean(resampled_data)
}

# Standard Normal Method
se <- sd(bootstrap_sample_means)  # Standard error of sample means
z_alpha_over_2 <- qnorm(1 - alpha / 2)  # Z quantile

standard_normal_lower <- mean(bootstrap_sample_means) - z_alpha_over_2 * se
standard_normal_upper <- mean(bootstrap_sample_means) + z_alpha_over_2 * se



# Basic Bootstrap Method
bias <- mean(bootstrap_sample_means) - mean(data)
bias_corrected_mean <- mean(bootstrap_sample_means) - bias
basic_bootstrap_se <- sd(bootstrap_sample_means - bias)  # Bias-corrected SE

basic_lower <- bias_corrected_mean - qnorm(1 - alpha / 2) * basic_bootstrap_se
basic_upper <- bias_corrected_mean + qnorm(1 - alpha / 2) * basic_bootstrap_se

# Percentile Method
percentile_lower <- quantile(bootstrap_sample_means, alpha / 2)
percentile_upper <- quantile(bootstrap_sample_means, 1 - alpha / 2)



# BCa (Bias-Corrected and Accelerated) Method

library(boot)
theta.b <- numeric(num_bootstrap_samples)
  theta.hat <- mean(data)
  #bootstrap
  for (b in 1:num_bootstrap_samples) {
    i <- sample(1:n, size = n, replace = TRUE)
    theta.b[b] <- mean(data[i])
  }
  
  stat <- function(dat, index) {
    mean(dat[index]) }
  
  boot.BCa <-
    function(x, th0, th, stat, conf = .95) {
   
      x <- as.matrix(x)
      n <- nrow(x) #observations in rows
      N <- 1:n
      alpha <- (1 + c(-conf, conf))/2
      zalpha <- qnorm(alpha)
      # the bias correction factor
      z0 <- qnorm(sum(th < th0) / length(th))
      # the acceleration factor (jackknife est.)
      th.jack <- numeric(n)
      for (i in 1:n) {
        J <- N[1:(n-1)]
        th.jack[i] <- stat(x[-i, ], J)
      }
      L <- mean(th.jack) - th.jack
      a <- sum(L^3)/(6 * sum(L^2)^1.5)
      # BCa conf. limits
      adj.alpha <- pnorm(z0 + (z0+zalpha)/(1-a*(z0+zalpha)))
      limits <- quantile(th, adj.alpha, type=6)
      return(list("est"=th0, "BCa"=limits))
    }
  
  op<-boot.BCa(data, th0 = theta.hat, th = theta.b, stat = stat)
  
  bca_lower<-op$BCa[1]
  bca_upper<-op$BCa[2]
  
  cat("Standard Normal Method CI: [", standard_normal_lower, ", ", standard_normal_upper, "]\n")
  
  
cat("Basic Bootstrap Method CI: [", basic_lower, ", ", basic_upper, "]\n")

cat("Percentile Method CI: [", percentile_lower, ", ", percentile_upper, "]\n")

cat("BCa Method CI: [", bca_lower, ", ", bca_upper, "]\n")


## ----collapse=TRUE------------------------------------------------------------
 library(bootstrap)

data(scor)
data<-scor

cov_matrix <- cov(data)
eigen_result <- eigen(cov_matrix)


eigenvalues <- eigen_result$values


theta <- max(eigenvalues) / sum(eigenvalues)

# 创建一个函数
jackknife_estimate <- function(data) {
  n <- length(data)
  theta_estimates <- numeric(n)
  
  for (i in 1:n) {
 
    pseudo_data <- data[-i]

    pseudo_theta <- min(pseudo_data) / sum(pseudo_data)
    
    theta_estimates[i] <- pseudo_theta
  }
  
 
  bias <- (n - 1) * mean(theta_estimates) - theta
  
  
  se <- sqrt(((n - 1) / n) * sum((theta_estimates - mean(theta_estimates))^2))
  
  return(list(bias = bias, se = se))
}


jackknife_result <- jackknife_estimate(eigenvalues)

# 输出Jackknife估计的偏差和标准误差
cat("Jackknife Estimate of Bias:", jackknife_result$bias, "\n")
cat("Jackknife Estimate of Standard Error:", jackknife_result$se, "\n")


## ----collapse=TRUE------------------------------------------------------------

library(DAAG)

calculate_prediction_error1 <- function(model, data) {
  predicted_values <- predict(model, newdata = data)
  actual_values <- data$magnetic
  error <- mean((predicted_values - actual_values)^2)  
  return(error)
}

calculate_prediction_error2 <- function(model, data) {
  predicted_values <- predict(model, newdata = data)
  actual_values <- data$magnetic
  error <- mean((exp(predicted_values) - actual_values)^2)  # Mean Squared Error
  return(error)
}

data(ironslag)
data<-ironslag
# Define the proposed models
model_linear <- lm(magnetic ~ chemical, data = data)
model_quadratic <- lm(magnetic ~ chemical + I(chemical^2), data = data)
model_exponential <- lm(log(magnetic) ~ chemical, data = data)
model_log_log <- lm(log(magnetic) ~ log(chemical), data = data)

# Set the number 
num_iterations <- nrow(data) * (nrow(data) - 1) / 2


errors_linear <- numeric(num_iterations)
errors_quadratic <- numeric(num_iterations)
errors_exponential <- numeric(num_iterations)
errors_log_log <- numeric(num_iterations)


set.seed(123)  
iteration <- 1
for (i in 1:(nrow(data) - 1)) {
  for (j in (i + 1):nrow(data)) {
    validation_data <- data[c(i, j), ]
    training_data <- data[-c(i, j), ]
    
    # Calculate prediction errors for each model
    errors_linear[iteration] <- calculate_prediction_error1(model_linear, validation_data)
    errors_quadratic[iteration] <- calculate_prediction_error1(model_quadratic, validation_data)
    errors_exponential[iteration] <- calculate_prediction_error2(model_exponential, validation_data)
    errors_log_log[iteration] <- calculate_prediction_error2(model_log_log, validation_data)
    
    iteration <- iteration + 1
  }
}

# Calculate mean prediction errors for each model
mean_error_linear <- mean(errors_linear)
mean_error_quadratic <- mean(errors_quadratic)
mean_error_exponential <- mean(errors_exponential)
mean_error_log_log <- mean(errors_log_log)


cat("Mean Prediction Errors (Leave-Two-Out Cross-Validation):\n")
cat("Linear Model:", mean_error_linear, "\n")
cat("Quadratic Model:", mean_error_quadratic, "\n")
cat("Exponential Model:", mean_error_exponential, "\n")
cat("Log-Log Model:", mean_error_log_log, "\n")




## ----collapse=TRUE------------------------------------------------------------
library(twosamples)
attach(chickwts)
x <- sort(as.vector(weight[feed == "soybean"]))
y <- sort(as.vector(weight[feed == "linseed"]))
detach(chickwts)


cv_statistic_observed <- cvm_test(x,y)[1]


pooled_data <- c( x, y)

num_permutations <- 1000

cv_statistics_permuted <- numeric(num_permutations)


set.seed(123)  
for (i in 1:num_permutations) {
  permuted_data <- sample(pooled_data, replace = FALSE)
  
  # Split the permuted data into two samples
  permuted_sample1 <- permuted_data[1:length(x)]
  permuted_sample2 <- permuted_data[(length(x) + 1):length(permuted_data)]
  
  # Calculate the Cramér-von Mises statistic for the permuted samples
  cv_statistics_permuted[i] <- cvm_test(permuted_sample1, permuted_sample2)[1]
}


p_value <- sum(cv_statistics_permuted >= cv_statistic_observed) / num_permutations

# Display the results
cat("原样本的 Cramér-von Mises 统计量:", cv_statistic_observed, "\n")
cat("Permutation Test p-value:", p_value, "\n")

## ----collapse=TRUE------------------------------------------------------------
maxout <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  return(max(c(outx, outy)))
}


n1 <- 20
n2 <- 50
mu1 <- mu2 <- 0
sigma1 <-1
sigma2 <-1
m <- 1000

p_list <- numeric(100)
for (j in 1:100)
{
  

x <- rnorm(n1, mu1, sigma1)
y <- rnorm(n2, mu2, sigma2)
cv_statistic_observed<-maxout(x,y)



pooled_data <- c(x, y)

num_permutations <- 1000

cv_statistics_permuted <- numeric(num_permutations)


for (i in 1:num_permutations) {
  permuted_data <- sample(pooled_data, replace = FALSE)
  
  # Split the permuted data into two samples
  permuted_sample1 <- permuted_data[1:length(x)]
  permuted_sample2 <- permuted_data[(length(x) + 1):length(permuted_data)]
  
  # Calculate the Cramér-von Mises statistic for the permuted samples
  cv_statistics_permuted[i] <- maxout(permuted_sample1,permuted_sample2)
}

p_list[j]=sum(cv_statistics_permuted > cv_statistic_observed) / num_permutations
}
mean(p_list)

cat("原样本的 Cramér-von Mises 统计量:", cv_statistic_observed, "\n")
cat("Permutation Test p-value:", p_value, "\n")

## ----collapse=TRUE------------------------------------------------------------

laplace_pdf <- function(x) {
  return(0.5 * exp(-abs(x)))
}


metropolis_sampler <- function(iterations, proposal_sd) {

  chain <- numeric(iterations)
  acceptance_count <- 0
  current_value <- rnorm(1)  
  

  for (i in 1:iterations) {
 
    proposal_value <- current_value + rnorm(1, mean = 0, sd = proposal_sd)
    

    acceptance_ratio <- min(1, laplace_pdf(proposal_value) / laplace_pdf(current_value))
 
    if (runif(1) < acceptance_ratio) {
      current_value <- proposal_value
      acceptance_count <- acceptance_count + 1
    }
    
    chain[i] <- current_value
  }
  
  # 计算接受率
  acceptance_rate <- acceptance_count / iterations
  
  return(list(chain = chain, acceptance_rate = acceptance_rate))
}


iterations <- 10000
proposal_sd_values <- c(0.1, 1, 10)  



par(mfrow = c(3, 1), mar = c(4, 4, 1, 1))


for (proposal_sd in proposal_sd_values) {
 
  result <- metropolis_sampler(iterations, proposal_sd)
  

  plot(result$chain, type = "l", col = "blue", xlab = "Iterations", ylab = "Chain Values",
       main = "", ylim = c(-4, 4))
  

  cat("Proposal SD:", proposal_sd, "\tAcceptance Rate:", result$acceptance_rate, "\n")
}


## -----------------------------------------------------------------------------

gibbs_sampler <- function(n, rho = 0.9) {
 
  x <- numeric(n)
  y <- numeric(n)
  

  x[1] <- rnorm(1)
  y[1] <- rnorm(1)
  
  
  for (t in 2:n) {
 
    x[t] <- rnorm(1, mean = rho * y[t - 1], sd = sqrt(1 - rho^2))
    
 
    y[t] <- rnorm(1, mean = rho * x[t], sd = sqrt(1 - rho^2))
  }
  
  # Return the bivariate normal chain
  return(data.frame(X = x, Y = y))
}


set.seed(123)


n <- 1000
chain <- gibbs_sampler(n, rho = 0.9)

plot(chain$X, chain$Y, type = "l", col = "blue", xlab = "X", ylab = "Y", main = "Bivariate Normal Chain")



lm_model <- lm(Y ~ X, data = chain)
abline(lm_model, col = "red")
residuals <- residuals(lm_model)


par(mfrow = c(1, 2))
hist(residuals, main = "Histogram of Residuals", col = "lightblue", probability = TRUE)
lines(density(residuals), col = "red")


qqnorm(residuals, main = "Q-Q Plot of Residuals")
qqline(residuals, col = "red")

cat("Intercept:", coef(lm_model)[1], "\n")
cat("beta:", coef(lm_model)[2], "\n")

## ----collapse=TRUE------------------------------------------------------------

library(coda)

f_computed<-function(x,sigma1)
{
  return(  (x)/sigma1^2 *exp(1)^{-x^2/(2*sigma1^2)} )
}

metropolis_hastings_rayleigh <- function(n, sigma_proposal) 
  {
  x <- numeric(n)
  x[1] <- 1 

  for (t in 2:n) {
    x_proposed <- x[t - 1] + rnorm(1, 0, sigma_proposal)
    
    acceptance_ratio <- min(1, f_computed(x_proposed,1)/f_computed(x[t - 1],1)  )
    
    if (runif(1) < acceptance_ratio)
  {
      x[t] <- x_proposed
    } else {
      x[t] <- x[t - 1]
    }
  }
  
  return(x)
}

set.seed(123)

n_samples <- 1000
sigma_proposal <- 0.5
chain <- metropolis_hastings_rayleigh(n_samples, sigma_proposal)

num_chains <- 3

chains <- lapply(1:num_chains, function(i) metropolis_hastings_rayleigh(n_samples, sigma_proposal))

mcmc_chains <- lapply(chains, as.mcmc)

mcmc_list <- mcmc.list(mcmc_chains)

gelman_diag <- gelman.diag(mcmc_list)

print(gelman_diag)

gelman.plot(mcmc_list)


## ----collapse=TRUE------------------------------------------------------------
# 观测数据
data <- matrix(c(11, 12, 8, 9, 27, 28, 13, 14, 16, 17, 0, 1, 23, 24, 10, 11, 24, 25, 2, 3), ncol = 2, byrow = TRUE)

# EM算法
lambda <- 1  

for (iter in 1:100) {  
  # E步
  E_values <- numeric(nrow(data))
   for (i in 1:nrow(data)) {
    integrand <- function(x) {
      x * lambda * exp(-lambda * x) / (exp(-lambda *data[i, 1])-exp(-lambda *data[i, 2]))
    }
    
    E_values[i] <- integrate(integrand, lower = data[i, 1], upper = data[i, 2])$value
  }
  
 
  
  lambda_new <- 1 / mean(E_values)
  

  if (abs(lambda_new - lambda) < 1e-6) {
    break
  } else {
    lambda <- lambda_new
  }
}


print(paste("最大似然估计：", 1 / mean(data)))


print(paste("EM算法估计：", lambda))


## -----------------------------------------------------------------------------
library(boot)


solve.game <- function(A) {
#solve the two player zero-sum game by simplex method
#optimize for player 1, then player 2
#maximize v subject to ...
#let x strategies 1:m, and put v as extra variable
#A1, the <= constraints
#
min.A <- min(A)
A <- A - min.A #so that v >= 0
max.A <- max(A)
A <- A / max(A)
m <- nrow(A)
n <- ncol(A)
it <- n^3
a <- c(rep(0, m), 1) #objective function
A1 <- -cbind(t(A), rep(-1, n)) #constraints <=
b1 <- rep(0, n)
A3 <- t(as.matrix(c(rep(1, m), 0))) #constraints sum(x)=1
b3 <- 1
sx <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=TRUE, n.iter=it)
#the ’solution’ is [x1,x2,...,xm | value of game]
#
#minimize v subject to ...
#let y strategies 1:n, with v as extra variable
a <- c(rep(0, n), 1) #objective function
A1 <- cbind(A, rep(-1, m)) #constraints <=
b1 <- rep(0, m)
A3 <- t(as.matrix(c(rep(1, n), 0))) #constraints sum(y)=1
b3 <- 1
sy <- simplex(a=a, A1=A1, b1=b1, A3=A3, b3=b3,
maxi=FALSE, n.iter=it)
soln <- list("A" = A * max.A + min.A,
"x" = sx$soln[1:m],
"y" = sy$soln[1:n],
"v" = sx$soln[m+1] * max.A + min.A)
soln
}


A <- matrix(c( 0,-2,-2,3,0,0,4,0,0,
2,0,0,0,-3,-3,4,0,0,
2,0,0,3,0,0,0,-4,-4,
-3,0,-3,0,4,0,0,5,0,
0,3,0,-4,0,-4,0,5,0,
0,3,0,0,4,0,-5,0,-5,
-4,-4,0,0,0,5,0,0,6,
0,0,4,-5,-5,0,0,0,6,
0,0,4,0,0,5,-6,-6,0), 9, 9)

s <- solve.game(A+2)

round(cbind(s$x, s$y), 7)


## -----------------------------------------------------------------------------

my_vector <- c(1, 2, 3, 4, 5)

dimensions <- dim(my_vector)

print(dimensions)


## -----------------------------------------------------------------------------
my_df <- data.frame(a = c(1, 2, 3), b = c("apple", "banana", "orange"))

my_matrix <- as.matrix(my_df)

str(my_matrix)

## -----------------------------------------------------------------------------
# 创建一个0行的数据框
empty_df_rows <- data.frame()

# 打印结果
print(empty_df_rows)

## -----------------------------------------------------------------------------

# 创建一个0列的数据框
empty_df_cols <- data.frame(matrix(nrow = 1, ncol = 0))

# 打印结果
print(empty_df_cols)

## -----------------------------------------------------------------------------
scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}


my_df <- data.frame(
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  C = c(7, 8, 9)
)

# 使用 apply 函数应用 scale01 到每一列
scaled_df <- apply(my_df, 2, scale01)


print(scaled_df)

library(dplyr)

my_df <- data.frame(
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  C = c(7, 8, 9),
  D = c("apple", "banana", "orange"),
  E = c(TRUE, FALSE, TRUE)
)


scaled_df <- my_df %>%
  mutate_if(is.numeric, scale01)


print(scaled_df)

## -----------------------------------------------------------------------------
numeric_df <- data.frame(
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  C = c(7, 8, 9)
)

std_devs <- vapply(numeric_df, sd, numeric(1))

print(std_devs)



mixed_df <- data.frame(
  A = c(1, 2, 3),
  B = c(4, 5, 6),
  C = c(7, 8, 9),
  D = c("apple", "banana", "orange"),
  E = c(TRUE, FALSE, TRUE)
)


numeric_std_devs <- vapply(mixed_df, function(x) if(is.numeric(x)) sd(x) else NA, numeric(1))


print(numeric_std_devs)

## -----------------------------------------------------------------------------
# 安装并加载 microbenchmark 包
# install.packages("microbenchmark")
library(microbenchmark)

# 定义R函数
gibbs_sampler_R <- function(n, a, b, iterations) {
  x <- numeric(iterations)
  y <- numeric(iterations)
  
  for (i in 1:iterations) {
    # 从条件分布中采样
    if(i==1)
    {
      s<-0
    }
    else{
      s<-x[i-1]
    }
    y[i] <- rbeta(1, s + a, n - s + b)
    x[i] <- rbinom(1, n, y[i])
  }
  
  return(data.frame(x, y))
}

# 运行 microbenchmark 比较计算时间
result_R <- microbenchmark(gibbs_sampler_R(10, 2, 3, 1000))
print(result_R)

## -----------------------------------------------------------------------------
library(Rcpp)
sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix gibbs_sampler_Rcpp(int n, double a, double b, int iterations)
{
  NumericVector x(iterations);
  NumericVector y(iterations);
  NumericMatrix out_stat(iterations,2);
  double s;

  for (int i = 0; i < iterations; i++) {
    if(i==0)
    {
      s=0;
    }
    else
    {
      s=x(i-1);
    }
    
    y(i) = R::rbeta(s + a, n - s + b);
    x(i) = R::rbinom( n, y(i));
    out_stat(i,0) = x(i);
    out_stat(i,1) = y(i);
  }
  
  return out_stat;
}
')

# 运行 microbenchmark 比较计算时间
result_Rcpp <- microbenchmark::microbenchmark(gibbs_sampler_Rcpp(10, 2, 3, 1000))
print(result_Rcpp)

