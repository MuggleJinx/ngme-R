
# compute_grad <- function(A, K, nu, mu, sigma, sigma_e, V, M, Y, h) {
#   g_mu <- g_sigma <- g_nu <- g_sigma_e <- g_rho <- 0
  
#   n <- length(Y)
  
#   tmp1 <- solve(K, (V - h))
#   tmp2 <- solve(K, M)
#   tmp3 <- tmp2 + mu * tmp1 # solve(K, M + mu * (V-h))

#   # residual <- Y - A %*% Kinv %*% (M + mu * (V-h))
#   residual <- Y - A %*% tmp3

#   # grad wrt mu  
#   # g_mu <- t(residual) %*% A %*% Kinv %*% (V-h) / sigma_e^2
#   g_mu <- t(residual) %*% A %*% (tmp1) / sigma_e^2
#   # print(paste0("g_mu: ", g_mu))

#   # grad wrt sigma
#   g_sigma <- sum(M^2 / V / sigma^3) - n / sigma

#   # grad wrt nu
#   g_nu <- sum(h - V + h * log(V) - h * log(1/nu) - h * digamma(h*nu))

#   # grad wrt sigma_e
#   g_sigma_e <- sum(residual^2) / sigma_e^3 - n / sigma_e
#   # print(paste0("g_sigma_e: ", g_sigma_e))

#   # grad wrt rho
#   # g_rho <- - t(residual) %*% A %*% Kinv %*% C %*% Kinv %*% (M + mu * (V-h)) / sigma_e^2
#   g_rho <- - t(residual) %*% A %*% solve(K, C %*% tmp3) / sigma_e^2
#   # print(paste0("g_rho: ", g_rho))

#   list(
#     mu      = - as.numeric(g_mu),
#     sigma   = - as.numeric(g_sigma),
#     nu      = - as.numeric(g_nu), 
#     sigma_e = - as.numeric(g_sigma_e),
#     rho     = - as.numeric(g_rho)
#   )
# }

solve_with_chol <- function(A, b, tol = 1e-8) {
  chol_result = tryCatch(
    chol(A),
    error = function(e) NULL
  )
  
  if(!is.null(chol_result)) {
    return(backsolve(chol_result, forwardsolve(t(chol_result), b)))
  } else {
    # 如果Cholesky分解失败，使用正则化
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

    cond_M <- solve_with_chol(Q, b)
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
    tr_sigma <- sum(diag(solve_with_chol(Q, diag(1/V)))) / sigma^3
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
    tr_sigma_e <- sum(diag(solve_with_chol(Q, t(B) %*% B))) / sigma_e^3
    g_sigma_e <- g_sigma_e + tr_sigma_e
  }
  g_log_sigma_e <- g_sigma_e * sigma_e # chain rule

  # grad wrt rho
  if (!RB) {
    g_rho <- - t(residual) %*% A %*% Kinv %*% C %*% Kinv %*% (M - mu*h) / sigma_e^2 
  } else {
    g_rho <- - t(residual_condW) %*% A %*% Kinv %*% C %*% Kinv %*% (cond_M - mu*h) / sigma_e^2 
    tr_rho <- - sum(diag(solve_with_chol(Q, t(B) %*% B %*% C %*% Kinv))) / sigma_e^2
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
    cond_W = solve_with_chol(QQ, tmp)
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
    
    tr_sigma = sum(diag(solve_with_chol(QQ, Q))) / sigma
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
    tr_sigma_e <- sum(diag(solve_with_chol(QQ, t(A) %*% A))) / sigma_e^3
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
    tr_rho <- - sum(diag(solve_with_chol(QQ, t(C) %*% diag(1/V) %*% K))) / sigma^2
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


# test non-centered gradient
test_grad <- function(
  true_params = list(nu=2, rho=0.5, mu=-3, sigma=4, sigma_e=0.5),
  start_params = list(nu=0.1, rho=0.5, mu=-3, sigma=4, sigma_e=0.5),
  lr = 0.01, 
  n_gibbs = 5,
  mode="non-centered",
  n_obs=100, 
  niter=1000, 
  nburn=100, 
  noise="gal", 
  seed=20,
  sampling=FALSE,
  RB=FALSE,
  n_chain=2
) {
  set.seed(seed)

  sim <- simulation_ar1(
    n_obs = n_obs,
    feff = c(-1, 2),

    rho = true_params$rho,
    mu = true_params$mu,
    sigma = true_params$sigma,
    nu = true_params$nu,
    sigma_e = true_params$sigma_e,

    noise = noise,
    seed = seed
  )

  est <- multi_chain_estimation(
    n_chain = n_chain,
    Y = sim$Y,
    A = sim$A,
    C = sim$C,
    G = sim$G,
    noise = noise,
    mode = mode,
    optimizer = "adam",
    n_gibbs = n_gibbs,
    lr = lr,
    RB = RB,

    # start parameters
    rho = start_params$rho,
    mu = start_params$mu,
    sigma = start_params$sigma,
    sigma_e = start_params$sigma_e,
    nu = start_params$nu,

    iteration = niter,
    burnin = nburn,
    sampling = sampling,
    true_M = if (!sampling) sim$M else NULL,
    true_V = if (!sampling) sim$V else NULL,
    true_W = if (!sampling) sim$W else NULL
  )

  return(est)
}