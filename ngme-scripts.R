# This file contains scripts for the experiments for comparing estimation with 
# centered and non-centered parameterization.

# load packages
library(tidyr)
library(ggplot2)
library(cowplot)
library(patchwork)
library(testthat)
library(ngme2)

source("R/gradient.R")
source("R/estimation.R")
source("R/util.R")
source("R/sampling.R")


# main test function
main_test <- function(
  nu = 2, 
  noise = "nig", 
  mode = "centered",
  sampling = TRUE, 
  RB = FALSE,
  niter = 1000, 
  n_chain = 4,
  true_params = list(rho=0.5, mu=-3, sigma=4, nu=nu, sigma_e=1),
  start_params = list(rho=0, mu=0, sigma=1, nu=1, sigma_e=0.5),
  lr = 0.01,
  n_gibbs = 5,
  nburn = 100,
  n_obs = 300,
  seed = 20,
  plot_title = FALSE
) {
  method <- if (RB) "RB" else "MC"

  # simulation
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

  # estimation
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
    iteration = niter,
    burnin = nburn,

    # start parameters
    rho = start_params$rho,
    mu = start_params$mu,
    sigma = start_params$sigma,
    sigma_e = start_params$sigma_e,
    nu = start_params$nu,

    sampling = sampling,
    true_M = if (!sampling) sim$M else NULL,
    true_V = if (!sampling) sim$V else NULL,
    true_W = if (!sampling) sim$W else NULL
  )
  print("after estimation")

  # print results
  print(paste("Noise=", noise, ", nu=", nu, ", mode=", mode))
  final <- est$final 
  names(final) <- c("rho", "mu", "sigma", "nu", "sigma_e")
  print(final)
  print(unlist(true_params))

  # traceplot
  main <- if (plot_title) { 
    paste0("Traceplot of ", noise, " model (nu=", nu, ") with ", mode, " parameterization using ", method, " gradient")
  } else {
    NULL
  }
  print("before traceplot")
  plot <- traceplot_df( 
    true_rho = true_params$rho, 
    true_mu = true_params$mu, 
    true_sigma = true_params$sigma, 
    true_nu = true_params$nu, 
    true_sigma_e = true_params$sigma_e, 

    rho_df = est$rho,
    mu_df = est$mu,
    sigma_df = est$sigma,
    nu_df = est$nu,
    sigma_e_df = est$sigma_e,

    main = main
  )
  print(plot)

  # save plot
  filename <- paste0(noise, "-nu=", nu, "-", mode, ".png")
  filename <- paste0("Figures/", method, "/", filename)

  print("before ggsave")
  ggsave(
    plot=plot, 
    filename=filename,
    width = 10, height = 5
  )
}

#########################
# GAL model
#########################

# MC, gal, nu=0.3
# for (mode in c("centered", "non-centered")) { 
#   main_test(
#     n_obs = 300,
#     nu = 1,
#     noise = "nig",
#     mode = "centered",
#     sampling = FALSE,
#     niter = 1000,
#     RB = FALSE,
#     n_gibbs = 5,
#     plot_title = FALSE,
#     seed = 32
#   )
# }

# MC, gal, nu=2
# for (mode in c("centered", "non-centered")) { 
#   main_test(
#     nu = 2, 
#     noise = "gal", 
#     mode = mode,
#     sampling = TRUE, 
#     niter = 3000, 
#     RB = FALSE, 
#     n_chain = 4,
#     plot_title = FALSE,
#     seed = 32
#   )
# }

# RB, gal, nu=2
# system.time({
#   for (mode in c("centered", "non-centered")) { 
#   main_test(
#     nu = 2, 
#     noise = "gal", 
#     mode = mode,
#     sampling = TRUE, 
#     niter = 3000, 
#     RB = TRUE, 
#     n_chain = 4,
#     plot_title = FALSE,
#     seed = 32
#     )
#   }
# })


# RB, gal, nu=0.5
system.time({
  for (mode in c("centered")) { 
  main_test(
    nu = 0.4, 
    noise = "gal", 
    mode = mode,
    sampling = TRUE, 
    niter = 2000, 
    RB = TRUE, 
    n_chain = 4,
    plot_title = FALSE,
    seed = 32
    )
  }
})


# random test 
# for (mode in c("centered", "non-centered")) { 
#   main_test(
#     nu = 2.1, 
#     noise = "nig", 
#     mode = mode,
#     # sampling = TRUE, 
#     sampling = FALSE, 
#     niter = 1000, 
#     RB = TRUE, 
#     n_chain = 4,
#     n_gibbs = 5,
#     plot_title = TRUE,
#     seed = 32
#   )
# }
