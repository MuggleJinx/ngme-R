source_all <- function() {
  library(ngme2)
  library(testthat)
  source("Code/gradient.R")
  source("Code/sampling.R")
  source("Code/estimation.R")
  source("Code/util.R")
}

test_that("Test single-chain estimation using ngme2", {
  source_all()
  set.seed(16)
  true_rho     <- 0.5
  true_mu      <- -3
  true_sigma   <- 4
  true_nu      <- 0.5
  true_sigma_e <- 0.2

  sim <- simulation_ar1(
    n_obs = 500,
    feff = c(-1, 2),
    rho = true_rho,
    mu = true_mu,
    sigma = true_sigma,
    nu = true_nu,
    sigma_e = true_sigma_e,
    noise = "nig",
    seed = 20
  )

  # Estimation using ngme2
  n_obs <- length(sim$Y)
  ngme2_fit <- ngme(
    Y ~ 0 + f(1:n_obs, model = "ar1", noise = noise_nig()),
    family = "normal",
    data = data.frame(Y = sim$Y),
    control_opt = control_opt(
      iteration = 2000,
      rao_blackwellization = TRUE
    )
  )
  ngme2_fit
  ngme2::traceplot(ngme2_fit, "field1")
  ngme2::traceplot(ngme2_fit)
})
