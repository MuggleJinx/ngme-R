solve_with_chol <- function(A, b, tol = 1e-8) {
  chol_result = tryCatch(
    chol(A),
    error = function(e) NULL
  )
  
  if(!is.null(chol_result)) {
    return(backsolve(chol_result, forwardsolve(t(chol_result), b)))
  } else {
    # return(solve(A + diag(tol, nrow(A)), b))
    return(MASS::ginv(A) %*% b)
  }
}


# M|V ~ N(mu V, sigma^2 V)
grad_non_centered <- function(
  RB, noise, A, params, V, M, Y, h, C, G
) {
  library(Matrix)
  n <- length(Y)

  th <- params[1]
  mu <- params[2]
  log_sigma <- params[3]
  log_nu <- params[4]
  log_sigma_e <- params[5]

  rho <- th_to_rho(th)
  sigma <- exp(log_sigma)
  sigma_e <- exp(log_sigma_e)
  nu <- exp(log_nu)
  g_mu <- g_log_sigma <- g_log_nu <- g_log_sigma_e <- g_th <- 0

  # update K
  K <- rho * C + G
  Kinv <- solve(K)

  residual <- Y - A %*% Kinv %*% (M - mu*h)

  if (RB) {
    B <- A %*% Kinv 
    Q <- t(B) %*% B / sigma_e^2 + diag(1/V) / sigma^2
    # b <- 1/sigma_e^2 * (t(B) %*% Y - mu * t(B) %*% B %*% h) +
    #   mu / sigma^2

    b <- 1/sigma_e^2 * t(B) %*% (Y + mu * B %*% h) + mu / sigma^2

    cond_M <- solve(Q, b)
    residual_condW <- Y - A %*% Kinv %*% (cond_M - mu*h)

    # Qinv <- solve(Q)
  }

  # grad wrt mu  
  if (!RB) {
    g_mu1 <- - t(residual) %*% A %*% Kinv %*% h / sigma_e^2
    g_mu2 <- sum(M - mu*V) / sigma^2
    g_mu <- (as.numeric(g_mu1) + g_mu2) 
  } else {
    g_mu1 <- - t(residual_condW) %*% A %*% Kinv %*% h / sigma_e^2
    g_mu2 <- sum(cond_M - mu*V) / sigma^2
    g_mu <- (as.numeric(g_mu1) + g_mu2) 
  }

  # grad wrt sigma
  if (!RB) {
    g_sigma <- sum((M^2/V + mu^2*V - 2*M*mu) / sigma^3)  - n/sigma
  } else {
    g_sigma <- sum((cond_M^2/V + mu^2*V - 2*cond_M*mu) / sigma^3)  - n/sigma
    tr_sigma <- sum(diag(solve(Q, diag(1/V)))) / sigma^3
    g_sigma <- g_sigma + tr_sigma
  }
  g_log_sigma <- g_sigma * sigma # chain rule

  # grad wrt nu
  g_nu <- grad_nu(noise, h, V, nu)
  g_log_nu <- g_nu * nu # chain rule

  # grad wrt sigma_e
  if (!RB) {
    g_sigma_e <- sum(residual^2) / sigma_e^3 - n / sigma_e
  } else {
    g_sigma_e <- sum(residual_condW^2) / sigma_e^3 - n / sigma_e
    tr_sigma_e <- sum(diag(solve(Q, t(B) %*% B))) / sigma_e^3
    g_sigma_e <- g_sigma_e + tr_sigma_e
  }
  g_log_sigma_e <- g_sigma_e * sigma_e # chain rule

  # grad wrt rho
  if (!RB) {
    g_rho <- - t(residual) %*% A %*% Kinv %*% C %*% Kinv %*% (M - mu*h) / sigma_e^2 
  } else {
    g_rho <- - t(residual_condW) %*% A %*% Kinv %*% C %*% Kinv %*% (cond_M - mu*h) / sigma_e^2 
    tr_rho <- - sum(diag(solve(Q, t(B) %*% B %*% C %*% Kinv))) / sigma_e^2
    g_rho <- g_rho + tr_rho
  }
  drho <- 2 * exp(th) / (1 + exp(th))^2
  g_th <- g_rho * drho

  c(
    - as.numeric(g_th),
    - as.numeric(g_mu),
    - as.numeric(g_log_sigma),
    - as.numeric(g_log_nu),
    - as.numeric(g_log_sigma_e)
  )
}



# ngme (centered) version of gradients
# KW|V ~ N(mu V - mu h, sigma^2 V)
grad_centered <- function(
  RB, noise, A, params, V, W, Y, h, C, G
) {
  library(Matrix)
  n <- length(Y)

  th <- params[1]
  mu <- params[2]
  log_sigma <- params[3]
  log_nu <- params[4]
  log_sigma_e <- params[5]

  rho <- th_to_rho(th)
  sigma <- exp(log_sigma)
  sigma_e <- exp(log_sigma_e)
  nu <- exp(log_nu)
  g_mu <- g_log_sigma <- g_log_nu <- g_log_sigma_e <- g_th <- 0

  # update K
  K <- rho * C + G
  Kinv <- solve(K)
  M = as.numeric(K %*% W)

  residual <- as.numeric(Y - A %*% W)

  # tmp = K' * inv(SV) * mu(V-h) + A' * inv(Sigma) * Y
  # cond_W = QQ^-1 * tmp
  tmp = t(K) %*% (1/(sigma^2*V) * mu*(V-h)) + t(A) %*% (Y / sigma_e^2)
  
  Q = t(K) %*% diag(1/(sigma^2*V)) %*% K
  QQ = Q + t(A) %*% A / sigma_e^2
  
  if (RB) {
    cond_W = solve(QQ, tmp)
    cond_M = as.numeric(K %*% cond_W)

    residual_condW <- as.numeric(Y - A %*% cond_W)
  }  

  # grad wrt mu  
  g_mu <- if (!RB) {
    sum((V - h) * (M + mu*(h-V)) / (sigma^2 * V))
  } else {
    sum((V - h) * (K %*% cond_W + mu*(h-V)) / (sigma^2 * V))
  }

  # grad wrt sigma
  g_sigma <- if (!RB) {
    sum((M^2 + 2*mu*M*(h-V) + mu^2 * h^2 - 2*mu^2*h*V + V*(mu^2*V-sigma^2)) / (V*sigma^3))
  } else {
    tmp = sum((cond_M^2 + 2*mu*cond_M*(h-V) + mu^2 * h^2 - 2*mu^2*h*V + V*(mu^2*V-sigma^2)) / (V*sigma^3))
    
    tr_sigma = sum(diag(solve(QQ, Q))) / sigma
    tmp + tr_sigma
  }
  g_log_sigma <- g_sigma * sigma # chain rule

  # grad wrt nu
  g_nu <- grad_nu(noise, h, V, nu)
  g_log_nu <- g_nu * nu # chain rule

  # grad wrt sigma_e
  if (!RB) {
    g_sigma_e <- sum(residual^2) / sigma_e^3 - n / sigma_e
  } else {
    g_sigma_e <- sum(residual_condW^2) / sigma_e^3 - n / sigma_e
    tr_sigma_e <- sum(diag(solve(QQ, t(A) %*% A))) / sigma_e^3
    g_sigma_e <- g_sigma_e + tr_sigma_e
  }
  g_log_sigma_e <- g_sigma_e * sigma_e # chain rule

  # grad wrt rho
  # tr = tr(K^-1 C)
  tr <- sum(diag(C %*% Kinv))
  # tr <- sum(diag(solve(K, C)))
  if (!RB) {
    g_rho <- tr - t(W) %*% t(C) %*% diag(1/V) %*% 
      (K%*%W - mu*(V-h)) / sigma_e^2
  } else {
    g_rho <- tr - t(cond_W) %*% t(C) %*% diag(1/V) %*% 
      (K%*%cond_W - mu*(V-h)) / sigma_e^2
    tr_rho <- - sum(diag(solve(QQ, t(C) %*% diag(1/V) %*% K))) / sigma^2
    g_rho <- g_rho + tr_rho
  }
  drho <- 2 * exp(th) / (1 + exp(th))^2
  g_th <- g_rho * drho

# print(paste0("g_mu: ", g_mu))
# print(paste0("g_sigma: ", g_sigma))
# print(paste0("g_nu: ", g_nu))
# print(paste0("g_log_sigma_e: ", g_log_sigma_e))
# print(paste0("g_rho: ", g_rho))

  c(
    - as.numeric(g_th),
    - as.numeric(g_mu),
    - as.numeric(g_log_sigma),
    - as.numeric(g_log_nu),
    - as.numeric(g_log_sigma_e)
  )
}


grad_nu <- function(noise, h, V, nu) {
  if (noise == "gal") {
    g_nu <- sum(h - V + h * log(V) - h * log(1/nu) - h * digamma(h*nu))
  } else if (noise == "nig") {
    g_nu <- sum(-h^2/(2*V) - V/2 + h + 1/(2*nu))
  }
  return(g_nu)
}

