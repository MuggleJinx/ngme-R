library(Matrix)

# define M = KW + mu * h
sample_post_V <- function(p, a, b, mu, sigma, h, M) {
  n <- max(length(p), length(a), length(b), length(mu), length(sigma), length(h), length(M))
  p <- p - 0.5
  a <- a + (mu / sigma)^2
  b <- b + ((M / sigma)^2)
  
  ngme2::rgig(n, p, a, b) 
}

# M|V, Y = KW + mu * h | V, Y
# check if it is correct
sample_post_M <- function(
  A, Kinv, mu, sigma, sigma_e, V, Y, h # X, beta
) {
  n <- length(Y)
  Q <- Matrix::diag(1 / (sigma^2 * V))
  Q_e <- Matrix::diag(rep(1 / sigma_e^2, n))
  QQ <- Q + t(Kinv) %*% t(A) %*% Q_e %*% A %*% Kinv
  
  b1 <- t(Kinv) %*% t(A) %*% Q_e %*% (-Y - A %*% Kinv %*% h)
  b2 <- mu * Q %*% V
  b <- - (b1 + b2)

  # M ~ N_c(b, QQ)
  U <- chol(QQ)
  w <- solve(t(U), b)
  mu <- solve(U, w)
  z <- rnorm(n)
  v <- solve(U, z)

  as.numeric(mu + v)
}

# W|V, Y
# KW|V ~ N(mu(V-h), SV)
# sample_post_W <- function(
#   A, K, mu, sigma, sigma_e, V, Y, h # X, beta 
# ) {
#   n <- length(Y)
#   inv_SV <- Matrix::diag(1 / (sigma^2 * V))
#   Q <- t(K) %*% inv_SV %*% K
#   QQ <- Q + t(A) %*% A / sigma_e^2

#   b <- t(K) %*% inv_SV %*% (mu * (V - h)) + t(A) %*% (Y/sigma_e^2)

#   # M ~ N_c(b, QQ)
#   U <- chol(QQ) 
#   w <- solve(t(U), b)
#   mu <- solve(U, w)
#   z <- rnorm(n)
#   v <- solve(U, z)

#   W <- as.numeric(mu + v)
#   W
# }
 


sample_post_W <- function(A, K, mu, sigma, sigma_e, V, Y, h, tol = 1e-8) {
  n <- length(Y)
  inv_SV <- Matrix::diag(1 / (sigma^2 * V))
  Q <- t(K) %*% inv_SV %*% K
  QQ <- Q + t(A) %*% A / sigma_e^2
  b <- t(K) %*% inv_SV %*% (mu * (V - h)) + t(A) %*% (Y/sigma_e^2)
  
  # 尝试稳健的Cholesky分解
  robust_chol_solve <- function(A, b, z, tol = 1e-8) {
    # 首先尝试直接Cholesky
    U <- tryCatch(
      chol(A),
      error = function(e) NULL
    )
    
    # 如果直接Cholesky失败，添加正则化项再试
    if (is.null(U)) {
      U <- tryCatch(
        chol(A + diag(tol, nrow(A))),
        error = function(e) NULL
      )
    }
    
    # 如果还是失败，使用SVD
    if (is.null(U)) {
      svd_result <- svd(A)
      d <- svd_result$d
      tol_svd <- tol * max(d)
      d_inv <- ifelse(d > tol_svd, 1/sqrt(d), 0)
      U <- t(svd_result$v %*% diag(sqrt(d)))
    }
    
    # 求解线性方程组
    w <- solve(t(U), b)
    mu <- solve(U, w)
    v <- solve(U, z)
    
    return(list(mu = mu, v = v))
  }
  
  # 生成随机数
  z <- rnorm(n)
  
  # 使用稳健求解
  result <- robust_chol_solve(QQ, b, z, tol)
  
  # 合并结果
  W <- as.numeric(result$mu + result$v)
  
  return(W)
}