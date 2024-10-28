# following similar structure as test-center-grad.R
source_all <- function() {
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(testthat)
  library(ngme2)

  source("Code/gradient.R")
  source("Code/sampling.R")
  source("Code/estimation.R")
  source("Code/util.R")
}

source_all()

main_test <- function(
  nu=2, 
  noise="nig", 
  mode = "non-centered",
  sampling=TRUE, 
  RB=FALSE,
  niter=1000, 
  n_chain=4,
  seed=20
) {
  method = if (RB) "RB" else "MC"

  test_that(paste("noise=", noise, ", nu=", nu, ", sampling=", sampling), {
    true_params <- list(rho=0.5, mu=-3, sigma=4, nu=nu, sigma_e=1)
    start_params <- list(rho=0, mu=0, sigma=1, nu=1, sigma_e=0.5)

    est <- test_grad(
      true_params = true_params,
      start_params = start_params,
      lr = 0.01, 
      n_gibbs = 5,
      mode = mode,
      n_obs = 300, 
      niter = niter, 
      nburn = 100, 
      noise = noise, 
      seed = seed,
      sampling = sampling,
      n_chain = n_chain,
      RB = RB
    ) 
    
    print(paste("Noise=", noise, ", nu=", nu, ", mode=", mode))
    final <- est$final 
    names(final) <- c("rho", "mu", "sigma", "nu", "sigma_e")
    print(final)
    print(unlist(true_params))

    # plot
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

      main = paste0("Traceplot of ", noise, " model (nu=", nu, ") with ", mode, " parameterization using ", method, " gradient")
    )
    plot

    # save_plot(
    #   plot=plot, 
    #   filename=paste0("Figures/nig-nu=", nu, "-", mode, ".png")
    # )
    filename <- paste0(noise, "-nu=", nu, "-", mode, ".png")
    if (RB) {
      filename <- paste0("Figures/RB/", filename)
    } else {
      filename <- paste0("Figures/MC/", filename)
    }

    ggsave(
      plot=plot, 
      filename=filename,
      width = 10, height = 5
    )

    # expect_true(abs(final["rho"] - true_params$rho) < 0.1)
    # expect_true(abs(final["mu"] - true_params$mu) < 1)
    # expect_true(abs(final["sigma"] - true_params$sigma) < 1)
    # expect_true(abs(final["nu"] - true_params$nu) < 0.5)
    # expect_true(abs(final["sigma_e"] - true_params$sigma_e) < 0.5)
  })
}

# main_test(nu=2, noise="nig", sampling=TRUE, niter=1500)
# main_test(nu=1, noise="nig", sampling=TRUE, niter=1500)
# main_test(nu=0.5, noise="nig", sampling=TRUE, niter=1500)
# main_test(nu=0.3, noise="nig", sampling=TRUE, niter=1500)

# main_test(nu=2, noise="gal", sampling=TRUE, niter=1500)
# main_test(nu=1, noise="gal", sampling=TRUE, niter=1500)
# main_test(nu=0.5, noise="gal", sampling=TRUE, niter=1500)
# main_test(nu=0.4, noise="gal", sampling=TRUE, niter=1500)
# main_test(nu=0.3, noise="gal", sampling=TRUE, niter=1500)
# main_test(nu=0.2, noise="gal", sampling=TRUE, niter=1500)

main_test(nu=0.45, noise="gal", sampling=TRUE, niter=1000)
main_test(nu=0.45, noise="gal", sampling=TRUE, niter=1000, RB=TRUE, seed = 30)

main_test(nu=0.4, noise="gal", sampling=TRUE, niter=500, RB=TRUE, seed = 28)
main_test(nu=0.35, noise="gal", sampling=TRUE, niter=500, RB=TRUE, seed = 29)

main_test(nu=0.5, noise="gal", sampling=TRUE, niter=3000, RB=TRUE)
main_test(nu=2, noise="gal", sampling=TRUE, niter=3000, RB=TRUE)

source_all()
main_test(
  nu=0.3, 
  mode = "non-centered",
  noise="gal", 
  sampling=TRUE, 
  niter=1000, 
  n_chain = 1,
  RB=TRUE, 
  seed = 29
)
