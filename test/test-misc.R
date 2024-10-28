source_all <- function() {
  library(ngme2)
  library(testthat)
  source("Code/gradient.R")
  source("Code/sampling.R")
  source("Code/estimation.R")
  source("Code/util.R")
}
source_all()

test_that("Test simulation", {
  set.seed(16)
  sim1 <- simulation_ar1(
    n_obs = 500,
    feff = c(-1, 2),
    rho = 0.5,
    mu = -3,
    sigma = 4,
    nu = 0.5,
    sigma_e = 0.5,
    noise = "gal"
  )

  sim2 <- simulation_ar1(
    n_obs = 500,
    feff = c(-1, 2),
    rho = 0.5,
    mu = -3,
    sigma = 4,
    nu = 0.5,
    sigma_e = 0.5,
    noise = "nig"
  )

  rho_est1 <- acf(sim1$W)$acf[2]
  rho_est2 <- acf(sim2$W)$acf[2]
  expect_true(round(rho_est1, 2) < 0.6 & round(rho_est1, 2) > 0.4)
  expect_true(round(rho_est2, 2) < 0.6 & round(rho_est2, 2) > 0.4)
})

test_that("Test single-chain estimation", {
  source_all()
  set.seed(16)
  true_rho     <- 0.5
  true_mu      <- -3
  true_sigma   <- 4
  true_nu      <- 0.5
  true_sigma_e <- 0.2

  sim <- simulation_ar1(
    n_obs = 300,
    feff = c(-1, 2),
    rho = true_rho,
    mu = true_mu,
    sigma = true_sigma,
    nu = true_nu,
    sigma_e = true_sigma_e,
    noise = "nig",
    seed = 20
  )

  source_all()
  system.time({
    est <- estimation(
      Y = sim$Y,
      A = sim$A,
      C = sim$C,
      G = sim$G,
      rho = 0,
      mu = true_mu,
      sigma = true_sigma,
      sigma_e = true_sigma_e,
      noise = "nig",
      mode = "centered",
      optimizer = "adam",
      RB = TRUE,
      n_gibbs = 5,
      lr = 0.01,
      iteration = 3000,
      burnin = 100,
      # true_M = sim$M,
      # true_V = sim$V,
      # true_W = sim$W
      sampling = TRUE
    )
  })

  traceplot_ggplot(est, 
    rho=true_rho, 
    mu=true_mu, 
    sigma=true_sigma, 
    nu=true_nu, 
    sigma_e=true_sigma_e
  )
})

test_that("Test multi-chain estimation", {
  set.seed(16)
  true_rho <- 0.7
  true_mu <- -3
  true_sigma <- 2
  true_nu <- 0.5
  true_sigma_e <- 0.5

  sim <- simulation_ar1(
    n_obs = 300,
    feff = c(-1, 2),
    rho = true_rho,
    mu = true_mu,
    sigma = true_sigma,
    nu = true_nu,
    sigma_e = true_sigma_e,
    noise = "nig"
  )

  source_all()
  system.time({
    est <- multi_chain_estimation(
      n_chain = 4,
      Y = sim$Y,
      A = sim$A,
      C = sim$C,
      G = sim$G,
      rho = 0,
      mu = true_mu,
      sigma = true_sigma,
      sigma_e = true_sigma_e,
      noise = "nig",
      mode = "centered",
      optimizer = "adam",
      # RB = TRUE,
      n_gibbs = 5,
      lr = 0.01,
      iteration = 100,
      burnin = 100,
      # true_M = sim$M,
      # true_V = sim$V,
      # true_W = sim$W
      sampling = TRUE
    )
  })

  source_all()
  traceplot_df(
    est$rho, est$mu, est$sigma, est$nu, est$sigma_e,
    true_rho = true_rho, true_mu = true_mu, true_sigma = true_sigma, true_nu = true_nu, true_sigma_e = true_sigma_e
  )

  ret <- est$final
  names(ret) <- c("rho", "mu", "sigma", "nu", "sigma_e")
  print(ret)
})

test_that("Test non-centered sampling (M, V)", {
  source_all()
  set.seed(16)
  sim <- simulation_ar1(
    n_obs = 100,
    feff = c(-1, 2),
    rho = 0.5,
    mu = -3,
    sigma = 4,
    nu = 0.5,
    sigma_e = 0.5,
    noise = "nig"
  )

  est <- estimation(
    Y = sim$Y,
    A = sim$A,
    C = sim$C,
    G = sim$G,
    noise = "nig",
    mode = "non-centered",
    optimizer = "adam",
    n_gibbs = 5,
    rho = 0,
    lr = 0.01,
    iteration = 3000,
    burnin = 100,
    sampling = TRUE
  )
  ret <- est$final
  names(ret) <- c("rho", "mu", "sigma", "nu", "sigma_e")
  print(ret)

  plot_trace(est, rho=0.5, mu=-3, sigma=4, nu=0.5, sigma_e=0.5)
  ret <- est$final
  names(ret) <- c("rho", "mu", "sigma", "nu", "sigma_e")
  print(ret)
})