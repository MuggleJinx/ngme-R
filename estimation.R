# Simulate data from the center model
simulation_ar1 <- function(
  n_obs,
  feff,
  rho,
  mu,
  sigma,
  nu,
  sigma_e,
  noise = "gal",
  seed = 16
) {
  # x1 <- runif(n_obs)
  # x2 <- rexp(n_obs)
  # X <- (model.matrix(~ 0 + x1 + x2))

  if (noise == "gal") {
    noise <- noise_gal(mu = mu, sigma = sigma, nu = nu)
  } else if (noise == "nig") {
    noise <- noise_nig(mu = mu, sigma = sigma, nu = nu)
  }

  day <- 1:n_obs
  ar1_model <- f(
    day,
    model = "ar1", 
    rho = rho,
    noise = noise
  )
  ar1_model

  s <- ngme2::simulate(ar1_model, seed = seed, nsim = 1)
  V <- attr(s, "V")[[1]]
  W <- s[[1]]

  # Y <- as.numeric(X %*% feff) + W + rnorm(n_obs, sd = sigma_e)
  Y <- W + rnorm(n_obs, sd = sigma_e)

  h <- rep(1, n_obs)
  A <- ar1_model$A
  C <- ar1_model$operator$C
  G <- ar1_model$operator$G
  K <- rho * C + G

  M <- as.numeric(K %*% W) + mu * h

  return(list(
    Y = Y, 
    M = M, 
    V = V, 
    W = W,
    A = A,
    C = C,
    G = G
  ))
}


# Estimation procedure, K = rho * C + G
estimation <- function(
  Y,
  A, 
  C, 
  G, 
  noise = "gal", 
  mode = "centered", # "centered" (same as ngme), "non-centered"
  optimizer = "adam",
  RB = FALSE,
  iteration = 1000, 
  burnin = 100,
  n_gibbs = 10,
  rho = 0.5, # AR1
  nu = 1, 
  mu = 0, 
  sigma = 1, 
  sigma_e = 1,
  lr = 0.01,
  true_M = NULL,
  true_V = NULL,
  true_W = NULL,
  sampling = FALSE,
  fix_params = c(rho=FALSE, mu=FALSE, sigma=FALSE, nu=FALSE, sigma_e=FALSE)
) {
  # adam parameters
  beta1 <- 0.9
  beta2 <- 0.999
  epsilon <- 1e-8
  m <- v <- rep(0, 5)

  # Convert parameters, order is important
  th <- rho_to_th(rho)
  mu <- mu 
  log_sigma <- log(sigma)
  log_nu <- log(nu)
  log_sigma_e <- log(sigma_e)
  params <- c(th, mu, log_sigma, log_nu, log_sigma_e)

  n <- length(Y)
  h <- rep(1, n)
  
  # Init V and M
  if (noise == "gal") {
    V <- rgamma(n, shape = h * nu, rate = nu)
  } else if (noise == "nig") {
    V <- rgig(n, p = rep(-0.5, n), a = rep(nu, n), b = nu * h^2)
  }

  # Init GIG parameters  
  noise_param <- convert_noise(noise, nu, h)
  # p, a, b are vectors of length n
  p <- noise_param$p
  a <- noise_param$a
  b <- noise_param$b

  K <- rho * C + G
  # Kinv <- solve(K)

  # Burnin
  for (i in 1:burnin) {
    # if (sample_M) M <- sample_post_M(A, Kinv, mu, sigma, sigma_e, V, Y, h)
    if (sampling) {
      W <- sample_post_W(A, K, mu, sigma, sigma_e, V, Y, h)
      M <- as.numeric(K %*% W) + mu * h
      V <- sample_post_V(p, a, b, mu, sigma, h, M)
    }
  }

  mu_trace <- sigma_trace <- nu_trace <- sigma_e_trace <- rho_trace <- double(iteration+1)
  mu_trace[1] <- mu; sigma_trace[1] <- sigma; nu_trace[1]<-nu;
  sigma_e_trace[1] <- sigma_e; rho_trace[1] <- rho

  # Fix M, V, W
if (!is.null(true_M)) M <- true_M
if (!is.null(true_V)) V <- true_V
if (!is.null(true_W)) W <- true_W
if (!sampling) n_gibbs <- 1

  # SGD
  for (i in 1:iteration) {
    print(paste0("iteration: ", i))
    g <- 0
    
    # Gibbs sampling
    for (j in 1:n_gibbs) {
      # Gibbs sampling
      if (sampling) {
        W <- sample_post_W(A, K, mu, sigma, sigma_e, V, Y, h)
        M <- as.numeric(K %*% W) + mu * h
        V <- sample_post_V(p, a, b, mu, sigma, h, M)
      }
      # if (mode == "non-centered" && sampling) {
      #   M <- sample_post_M(A, Kinv, mu, sigma, sigma_e, V, Y, h)
      #   V <- sample_post_V(p, a, b, mu, sigma, h, M)
      # }

      if (!is.null(true_M)) M <- true_M
      if (!is.null(true_V)) V <- true_V
      if (!is.null(true_W)) W <- true_W

      ###### compute gradient ######
      if (mode == "centered")
        g <- g + grad_centered(RB, noise, A, params, V, W, Y, h, C, G) / n_gibbs
      if (mode == "non-centered")
        g <- g + grad_non_centered(RB, noise, A, params, V, M, Y, h, C, G) / n_gibbs
      
      if (fix_params[1])  g[1] <- 0
      if (fix_params[2])  g[2] <- 0
      if (fix_params[3])  g[3] <- 0
      if (fix_params[4])  g[4] <- 0
      if (fix_params[5])  g[5] <- 0
    }

    # update using adam
    if (optimizer == "adam") {
      m <- beta1 * m + (1 - beta1) * g
      v <- beta2 * v + (1 - beta2) * g^2
      m_hat <- m / (1 - beta1^i)
      v_hat <- v / (1 - beta2^i)
      params <- params - lr * m_hat / (sqrt(v_hat) + epsilon)
    } else if (optimizer == "vanilla") {
      params <- params - lr * g
    } else {
      stop("Unknown optimizer")
    }

    # update parameters
    rho <- th_to_rho(params[1])
    mu <- params[2]
    sigma <- exp(params[3])
    nu <- exp(params[4])
    sigma_e <- exp(params[5])
    # update K
    K <- rho * C + G

    # update noise parameters p, a, b
    noise_param <- convert_noise(noise, nu, h)
    p <- noise_param$p
    a <- noise_param$a
    b <- noise_param$b

    # update trace
    rho_trace[i+1]     <- rho 
    mu_trace[i+1]      <- mu
    sigma_trace[i+1]   <- sigma
    nu_trace[i+1]      <- nu
    sigma_e_trace[i+1] <- sigma_e
  }

  list(
    rho = rho_trace,
    mu = mu_trace,
    sigma = sigma_trace,
    nu = nu_trace,
    sigma_e = sigma_e_trace,
    final = c(rho_trace[iteration+1], mu_trace[iteration+1], sigma_trace[iteration+1], nu_trace[iteration+1], sigma_e_trace[iteration+1])
  )
}

multi_chain_estimation <- function(n_chain, ...) {
  ret <- parallel::mclapply(1:n_chain, function(i) estimation(...), mc.cores = parallel::detectCores())
  
  
  # Combine results
  rho_df     <- do.call(cbind, lapply(ret, function(x) x$rho))
  mu_df      <- do.call(cbind, lapply(ret, function(x) x$mu))
  sigma_df   <- do.call(cbind, lapply(ret, function(x) x$sigma))
  nu_df      <- do.call(cbind, lapply(ret, function(x) x$nu))
  sigma_e_df <- do.call(cbind, lapply(ret, function(x) x$sigma_e))
  n_iter <- nrow(rho_df)

  final_rho <- mean(rho_df[n_iter, ])
  final_mu <- mean(mu_df[n_iter, ])
  final_sigma <- mean(sigma_df[n_iter, ])
  final_nu <- mean(nu_df[n_iter, ])
  final_sigma_e <- mean(sigma_e_df[n_iter, ])

  rho_df <- tidyr::pivot_longer(as.data.frame(rho_df), cols = everything(), names_to = "chain", values_to = "rho")
  rho_df$iteration <- rep(1:n_iter, each=n_chain)

  mu_df <- tidyr::pivot_longer(as.data.frame(mu_df), cols = everything(), names_to = "chain", values_to = "mu")
  mu_df$iteration <- rep(1:n_iter, each=n_chain)

  sigma_df <- tidyr::pivot_longer(as.data.frame(sigma_df), cols = everything(), names_to = "chain", values_to = "sigma")
  sigma_df$iteration <- rep(1:n_iter, each=n_chain)

  nu_df <- tidyr::pivot_longer(as.data.frame(nu_df), cols = everything(), names_to = "chain", values_to = "nu")
  nu_df$iteration <- rep(1:n_iter, each=n_chain)

  sigma_e_df <- tidyr::pivot_longer(as.data.frame(sigma_e_df), cols = everything(), names_to = "chain", values_to = "sigma_e")
  sigma_e_df$iteration <- rep(1:n_iter, each=n_chain)

  return(list(
    rho     = rho_df,
    mu      = mu_df,
    sigma   = sigma_df,
    nu      = nu_df,
    sigma_e = sigma_e_df,
    final = 
      c(
        rho = final_rho,
        mu = final_mu,
        sigma = final_sigma,
        nu = final_nu,
        sigma_e = final_sigma_e
      )
  ))
}
