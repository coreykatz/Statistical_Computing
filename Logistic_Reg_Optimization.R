library(tidyverse)
library(grid)

gradient_descent_1d <- function(x,f,df,t=1,alpha = .25,beta=.75,tol = 1e-6){
  ### Input: function, dervative, initial guess, ...
  ### Output: Minimum of f
  
  error <- c()
  for (i in 1:1000){
    dx <- -df(x[i]) #delta x is the negative of the gradient
    while (f(x[i]+t*dx)>(f(x[i])+alpha*t*df(x[i])* dx)){#line search
      t <- beta * t 
    }
    x[i+1] <- x[i] + t * dx #update x
    error[i] <- abs(f(x[i+1]) - f(x[i]))#check stopping condition
    if (error[i] < tol){return(list(x=x,f_x = f(x),error=error))}
  }
}


my_gradient_descent<- function(f,df,x,t=1,alpha = .25,beta=.75,tol = 1e-6){
  ### This function works on any low dim function
  ### Input: function, dervative, initial guess, ...
  ### Output: Minimum of f
  error <- c()
  x <- as.matrix(x,ncol=1)
  for (i in 1:1000){
    dx <- -df(x[,i])
    while (f(x[,i]+t*dx)>(f(x[,i])+alpha*t*t(-dx) %*% dx)){
      t <- beta * t
    }
    x <- cbind(x,x[,i] + t * dx)
    error[i] <- abs(f(x[,i+1])-f(x[,i]))
    if (error[i] < tol){return(list(x=matrix(x[,i+1]),error=error))}
  }
}


## Example
f <- function(x) {-log(1-x)-log(1-x^2)}
df <- function(x) {(1/(1-x))+(1/(1-x^2))*2*x}
d2f <- function(x) {-(-3*x^2-2*x-3)/(((x+1)^2)*((x-1)^2))} 
d2f(.5)
my_gradient_descent(f,df,x=-.1)
res_gd<- gradient_descent_1d(f,df,x=-.1)
df(0)^2 * 1/d2f(0)


## Newton's Method 1D
newton_method_1d <- function(x,f,df,d2f,t=1,alpha = .25,beta=.75, tol = 1e-6){
  ### Input: function, dervative,second dervative, initial guess, ...
  ### Output: Minimum of f
  dx <- -(1/d2f(x))*df(x) # change in x
  lambda2 <- (df(x)^2) * (1/d2f(x)) # first Newton decrement
  error <- c()
  j= 1
  while(lambda2/2 > tol){
    
    while (f(x[j]+t*dx)>(f(x[j])+alpha*t*df(x[j])* dx)) {
      t <- beta * t} #line search
    x[j+1] <- x[j] + t * dx #Updating x
    dx <- -(1/d2f(x[j+1]))*df(x[j+1])#calculate change in x
    error[j] <- abs(f(x[j+1]) - f(x[j]))
    lambda2 <- (df(x[j+1])^2) * (1/d2f(x[j+1])) #calculate Newton decrement (stopping criteria)
    j = j+1
  }
  x_star <- x[length(x)]
  return(list(x_star= x_star,f_star = f(x_star),x=x,f_x = f(x),error=error))
}

# Example from above
res1_nm <- newton_method_1d(x=-.1,f,df,d2f)


## Sigmoid Function
sigmoid <- function(x){
  return(1/(1+exp(-x)))
}


## Data example for regression
library(readr)
pima <- read_csv("~/Downloads/pima.csv")
X <- as.matrix(pima[,c(3,6,8)])
X <- scale(X) #scale and center X
X <- cbind(1,X)
colnames(X)[1] <- "Intercept"
y <- as.matrix(pima[,9])

# beta <-matrix(c(0.7460   ,  -0.1765   ,   0.8191   ,   0.5771),nrow=4)
# log(1+exp(X %*% beta))


#logistic score (gradient for logistic regression)
logit_reg_score <- function(X,y,beta){
  SUM <- matrix(0,nrow = length(beta),ncol=1)
  for (i in 1:nrow(X)) {
    SUM <- SUM + as.numeric((-y[i]+sigmoid(X[i,]%*%beta))) * X[i,]
  }
  return(-SUM)
}


#likelihood for logistic regression (nll)
nll_logit_reg <- function(X,y,beta){
  SUM <- 0
  for (i in 1:nrow(X)) {
    SUM <- SUM + y[i]*X[i,]%*%beta + log(1-sigmoid(X[i,]%*%beta))
  }
  return(-SUM)
}


##Observed Information Matrix for Logistic Regression
logit_reg_obs_info <- function(X,y,beta){
  W <- diag(as.numeric(sigmoid(X %*% beta) * (1-sigmoid(X %*% beta))))
  info <- -t(X) %*% W %*% X
  return(-info)
}
  

# nll_logit_reg(X,y,beta)
# a <-c()
# b<- c()
# for (i in 1:nrow(X)) {
#   a[i] <- y[i]*X[i,]%*%beta + log(1-sigmoid(X[i,]%*%beta))
#   b[i] <- y[i] * log(sigmoid(X[i,] %*% beta))+ (1-y[i])*log(1-sigmoid(X[i,]%*%beta)) 
# }
# 
# d <- data.frame(a,b)
# nll_logit_reg(X,y,beta)
# logit_reg_score(X,y,beta)
#logLik(glm(y~X-1,family = binomial))

# Logistic Regression Function using GD or Newton's Method
logit_reg <- function(X,y,beta,t=1,alpha = .1,Beta=.5,x,tol = 1e-6, newton = FALSE) {
  beta_old <- beta
  P <- diag(nrow(beta)) #Identiy P for GD
  for (i in 1:1000){
    if (newton==TRUE){
      P <- logit_reg_obs_info(X,y,beta) #observed information for Newton's Method
    }
    U <- logit_reg_score(X,y,beta) #negative score
    dx <- solve(P) %*% U ## change in x
    while (nll_logit_reg(X,y,beta+t*dx)>(nll_logit_reg(X,y,beta)+alpha*t*t(-dx) %*% dx)){
      t <- Beta * t
    } #line search
    beta <- beta + t * dx #update beta
    error <- norm(beta-beta_old)
    if(error <= tol){return(beta)}
    beta_old <- beta
  }
}

#Example Pima Data
beta_newton <- logit_reg(X,y,matrix(1,4,1),1,newton = T)
beta_gd <- logit_reg(X,y,matrix(1,4,1),1,newton = F)
glm(y~X-1,family = binomial) #glm function
# diag(nrow(U)) %*% U
# 
# 
# x <- c(0,0)
# beta <- matrix(1,4,1)

#Iterative Re-Weighted Least Squares for Logistic Regression
IWLS_logit <- function(X,y,beta,tol = 1e-6){
  for (i in 1:1000){
    eta <- X %*% beta[,i]#find eta-hat
    mu <- sigmoid(eta) #find mu-hat
    s <- (mu*(1-mu)) #variance of bernoulli 
    Z <- eta + ((y-mu)/s)#new response vector
    W <- diag(as.numeric(s))#weight matrix
    beta <- cbind(beta,solve(t(X) %*% W %*% X) %*% (t(X) %*% W %*% Z)) #update beta
    if (norm(matrix(abs(beta[,i+1]-beta[,i]),ncol = 1)) < tol){return(beta=beta[,ncol(beta)])}
 }
}

beta_iwls <- IWLS_logit(X,y,matrix(0,4,1))

# f <- function(x){(x[2]+2)^2 + (x[1]-3)^2}
# df <- function(x){
#   df_x1 <- 2*(x[1]-3)
#   df_x2 <- 2*(x[2]+2)
#   df <- matrix(c(df_x1,df_x2),ncol = 1)
#   return(df)
# }
# # df(x)
# 
# my_gradient_descent(f,df,x=x)
# f(c(0,1))
# f(c(2,1))
# 
# optim(par = c(.5,.5),fn = f)

## Plots of Convergence ####################



# dat <- data.frame(x=seq(-.5,.25,.001),f_x = f(seq(-.5,.25,.001)), g =1)
# dat2 <- data.frame(x=res_gd$x,f_x=res_gd$f_x, g=2)
# dat3 <- data.frame(x=res1_nm$x,f_x = res1_nm$f_x )

# ggplot(dat3,aes(x=x,y=f_x,group=factor(g),color=factor(g)))+
#   geom_path()

# 
# ggplot()+
#   geom_line(data = dat, aes(x=x,y=f_x),color = "blue")+
#   geom_segment(data=dat2, aes(x=x,y=f_x,xend=c(tail(dat2$x, n=-1), NA), yend=c(tail(dat2$f_x, n=-1), NA)),
#              arrow=arrow(length=unit(0.2,"cm")))
# 
# 
# ggplot()+
#   geom_line(data = dat, aes(x=x,y=f_x),color = "blue")+
#   geom_segment(data=dat3, aes(x=x,y=f_x,xend=c(tail(dat3$x, n=-1), NA), yend=c(tail(dat3$f_x, n=-1), NA)),
#                arrow=arrow(length=unit(0.2,"cm")))

