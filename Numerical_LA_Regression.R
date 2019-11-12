#Gram-Schimdt Factorization
gs_factorization <- function(A){
  ### Input: a m-by-n matrix (A)
  ### Output: an m-by-n (Q) and n-by-n (R) matrix where A=QR 
  
  
   ## Set dimensions
  m <- nrow(A) 
  n <- ncol(A)
  
  ## Set up empty Q and R Matrices of the correct dimension
  Q <- matrix(0,m,n)
  R <- matrix(0,n,n)
  
  ### First vector for #Gram-Schimdt Factorization
  v <- A[,1]
  R[1,1] <- sqrt(sum(v^2))
  Q[,1] <- v/R[1,1]

  
  for (j  in 2:n){
    v <- A[,j]
    for (i in 1:(j-1)){
      R[i,j] <- t(Q[,i]) %*% A[,j]
      v <- v- (R[i,j] * Q[,i]) #new basis vector
    }
    R[j,j] <- sqrt(sum(v^2))
    Q[,j] <- v/R[j,j]
  }
  
  return(list(R= R,Q = Q)) ## returns R and Q
} 


## Uses Gram-Schimdt to do linear regression
gs_regression <- function(X,y){
  ### Input: Design Matrix of any dimension (X) and the response variable(y)
  ### Output: Beta estimates and the Q,R from Gram-Schimdt Factorization
  
  
  Q <- gs_factorization(X)$Q #find Q from Gram-Schimdt
  R <- gs_factorization(X)$R #find R from Gram-Schimdt
  
  # Solve for betas
  b <- t(Q) %*% y 
  A <- R
  beta <- backsub(A,b)
  return(list(beta=beta,Q,R))
} 



## does back sub to solve systems of equations where A is upper triangular 
backsub <- function(A,b){
  ### Input: an n-by-n matrix on the left side of the linear system of 
  ###        equations and the vector on the right side of the linear
  ###        system of equation Ax=b
  ### Output: The solution to the system of equations (x)
  
  n <- ncol(A)
  x <- c()
  x[n] <- b[n]/A[n,n]
  
  for (i in (n-1):1){
    p <- A[i,(i+1):n] %*% x[(i+1):n]
    x[i] <- (b[i] - p)/A[i,i] #solves for x element by element
  }
  
  return(x)
}

## Does Back sub for vectors of matrices (defaults to finding inverses of upper triangular matrices)
backsub_matrix <- function(A,B = diag(nrow(A))){
  
  n <- ncol(A)
  if (class(B) == 'numeric'){B <- matrix(B,ncol = 1)}
  x <- B * 0
  for (j in 1:ncol(B)) {
    x[n,j] <- B[n,j]/A[n,n]
    for (i in (n-1):1){
      p <- A[i,(i+1):n] %*% x[(i+1):n,j]
      x[i,j] <- (B[i,j] - p)/A[i,i]
    }
  }
  
  return(x)
} 


# backsub_matrix(R)
# gs_regression(X,y)
# lm(y~X-1)

## QR Factorization using householder algorithm
hh_factorization <- function(X){
  ### Input: a m-by-n matrix (X)
  ### Output: an m-by-m (Q) and m-by-n (R) matrix where A=QR 
  
  # Initalize
  n <- nrow(X)
  p <- ncol(X)
  U <- matrix(0,n,p)
  H <- list()
  
  ## run through the algorithm
  for (k in 1:p) {
    w <- X[k:n,k]
    ## Find perpendicular vector
    w[1] <- w[1] + sign(w[1])* (sqrt(sum(w^2))) 
    u <- w/(sqrt(sum(w^2)))
    U[k:n,k] <- u
    g <- 2* (u %*% t(u)) 
    h <-  diag(length(u)) - g # calcualte the refelctor matrix
    h2 <- diag(n)
    h2[k:n,k:n] <- h
    H[[k]] <- h2 
    X[k:n,k:p] = X[k:n,k:p] - ((2*u) %*% (t(u) %*% X[k:n,k:p]))
    
  }
  R <- round(X, digits=10) #Final X matrix, all reflectors applied to X
  return(list(U = U,R=R,H=H))
} 
hh_factorization_reg <- function(X,Y){
  ### Input: Design Matrix of any dimension (X) and the response variable(y)
  ### Output: Beta estimates and the Q,R from Householder Factorization
  
  ## Same householders algorithm as above
  n <- nrow(X)
  p <- ncol(X)
  U <- matrix(0,n,p)
  for (k in 1:p) {
    w <- X[k:n,k]
    w[1] <- w[1] + sign(w[1])* (sqrt(sum(w^2)))
    u <- w/(sqrt(sum(w^2)))
    U[k:n,k] <- u
    X[k:n,k:p] = X[k:n,k:p] - ((2*u) %*% (t(u) %*% X[k:n,k:p]))
    Y[k:n] <- Y[k:n] - ((2*u) %*% (t(u) %*% Y[k:n])) ## updates Y at each step
  }
  R <- round(X, digits=10)
  beta <- backsub(R,Y) #solves the system to get beta estimates
  return(beta)
} ## Use householder algorithm of QR factorization to do linear regression


jacobi_reg <- function(X,y, beta_0 = rep(0,ncol(X)),tol = 1e-6, Nmax= 100){
  ### Input: Design Matrix of any dimension (X), the response variable(y), initial beta
  ### Output: Beta estimates by Jacobi Method
  
  b <- t(X) %*% y
  beta_old <- beta_0
  for (i in 2:Nmax){
    A <- t(X) %*% X
    Pinv<- diag(1/diag(A))
    I <- diag(nrow(A))
    beta <- ((I-(Pinv %*% A)) %*% beta_old) + Pinv %*% b #iteratively updating beta
    if(norm(beta-beta_old) < tol){return(beta)}
    beta_old <- beta
  }
  
} ## Uses jacobi to do linear regression
jacobi_reg(X,y)



##############################Recursive Least Squares

find_inv_XtX <- function(X){
  ### Input: Square Matrix (X) form of X^TX as in regression
  ### Output: inverse of X
  
  R <- gs_factorization(X)$R ## Does G-S factorization to get R
  Rinv <- backsub_matrix(R) 
  Rtinv <- t(Rinv)
  XtXinv <- Rinv %*% Rtinv #X = (R^TR)^-1
  return(XtXinv)
} ## Finds inverse of matrix of the form XtX

find_inv <- function(X){
  ### Input: Any Square Matrix (X)
  ### Output: inverse of X
  QR <- gs_factorization(X) ## Does G-S factorization to get R and Q
  Q <- QR$Q
  R <- QR$R
  Rinv <- backsub_matrix(R)
  
  # Use R and Q and their respective properties to find X^-1
  Xinv <- Rinv %*% t(Q)
  # Xinv <- solve(R) %*% t(Q)
  return(Xinv)
} ## finds the inverse of any matrix

recursive_reg <- function(X_n,X_k,y_n,y_k,beta_n){
  ### Input: Initial Design Matrix (X_n), New Design Matrix(X_k)
  ###        Initial response vector (y_n), New response vector (y_k)
  ###        Initial beta found through G-S regression
  ### Output: inverse of X
  
  ## Using Block Matrix Inverese to update inverse and calculate beta estimates recursively
  Ainv <- find_inv_XtX(X_n)
  B <- t(X_k)
  C <- diag(nrow(X_k))
  D <- X_k
  E <- (C+(D %*% Ainv %*% B))
  Einv <- find_inv(E)
  Rnkinv <- Ainv - (Ainv %*% B) %*% Einv %*% (D %*% Ainv)
  beta_nk <- beta_n + (Rnkinv %*% t(X_k) %*% (y_k - (X_k %*% beta_n)))
  return(beta_nk)
} ## Recursive Linear Regression

# beta_n <- gs_regression(X[1:200,],y[1:200])$beta
# recursive_reg(X[1:200,],X[201:300,],y[1:200],y[201:300],beta_n)


