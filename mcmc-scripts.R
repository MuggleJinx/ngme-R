set.seed(20)
library(ngme2)
source("R/MCMC.R")

#########################
# Trace-class
#########################

# Histogram of thining samples with true marginal density
niter = 20000
nburnin = 1000
p1 = compare_MC_MCMC(
  p = -0.5, a = 1, b = 0.1,
  h = 1, mu = 1, sigma = 1,
  niter = niter, nburnin = nburnin,
  xlim = c(0, 3),
  print_var = c("p", "a", "b")
)

p2 = compare_MC_MCMC(
  p = -0.5, a = 1, b = 0.5,
  h = 1, mu = 1, sigma = 1,
  niter = niter, nburnin = nburnin,
  xlim = c(0, 3),
  print_var = c("p", "a", "b")
)

p3 = compare_MC_MCMC(
  p = -0.5, a = 1, b = 1,
  h = 1, mu = 1, sigma = 1,
  niter = niter, nburnin = nburnin,
  xlim = c(0, 3),
  print_var = c("p", "a", "b")
)

p4 = compare_MC_MCMC(
  p = 0.55, a = 1, b = 0,
  h = 1, mu = 1, sigma = 1,
  niter = niter, nburnin = nburnin,
  xlim = c(0, 6), 
  print_var = c("p", "a", "b")
)

p5 = compare_MC_MCMC(
  p = 0.7, a = 1, b = 0,
  h = 1, mu = 1, sigma = 1,
  niter = niter, nburnin = nburnin,
  xlim = c(0, 6), 
  print_var = c("p", "a", "b")
)

p6 = compare_MC_MCMC(
  p = 1, a = 1, b = 0,
  h = 1, mu = 1, sigma = 1,
  niter = niter, nburnin = nburnin,
  xlim = c(0, 6),
  print_var = c("p", "a", "b")
)

library(patchwork)
histogram_tc <- p1 / p2 / p3 / p4 / p5 / p6 + plot_layout(ncol=3)
ggsave(histogram_tc, filename = "Figures/histogram_tc.png", width = 10, height = 5)


p1 = traceplot_mean(
  p = -0.5, a = 1, b = 0.1, 
  print_var = c("p", "a", "b")
)
p2 = traceplot_mean(
  p = -0.5, a = 1, b = 0.5, 
  print_var = c("p", "a", "b")
)
p3 = traceplot_mean(
  p = -0.5, a = 1, b = 1, 
  print_var = c("p", "a", "b")
)
p4 = traceplot_mean(
  p = 0.55, a = 1, b = 0, 
  print_var = c("p", "a", "b")
)
p5 = traceplot_mean(
  p = 0.7, a = 1, b = 0, 
  print_var = c("p", "a", "b")
)
p6 = traceplot_mean(
  p = 1, a = 1, b = 0, 
  print_var = c("p", "a", "b")
)

traceplot_tc <- p1 / p2 / p3 / p4 / p5 / p6 + plot_layout(ncol=3)
ggsave(traceplot_tc, filename = "Figures/traceplot_tc.png", width = 10, height = 5)


#########################
# Non-trace-class
#########################
p1 = compare_MC_MCMC(
  p = 0.2, a = 1, b = 0,
  h = 1, mu = 1, sigma = 1,
  niter = niter, nburnin = nburnin,
  xlim = c(0, 3), print_var = c("p", "a", "b")
)

p2 = compare_MC_MCMC(
  p = 0.3, a = 1, b = 0,
  h = 1, mu = 1, sigma = 1,
  niter = niter, nburnin = nburnin,
  xlim = c(0, 3), print_var = c("p", "a", "b")
)

p3 = compare_MC_MCMC(
  p = 0.5, a = 1, b = 0,
  h = 1, mu = 1, sigma = 1,
  niter = niter, nburnin = nburnin,
  xlim = c(0, 3), print_var = c("p", "a", "b")
)

p4 = compare_mcmc_invGamma(1.1, 1)
p5 = compare_mcmc_invGamma(1.5, 1)
p6 = compare_mcmc_invGamma(2, 1)

library(patchwork)
histogram_ntc <- p1 / p2 / p3 / p4 / p5 / p6 + plot_layout(ncol=3)
ggsave(histogram_ntc, filename = "Figures/histogram_ntc.png", width = 10, height = 5)


p1 = traceplot_mean(p = 0.2, a = 1, b = 0, print_var = c("p", "a", "b"))
p2 = traceplot_mean(p = 0.3, a = 1, b = 0, print_var = c("p", "a", "b"))
p3 = traceplot_mean(p = 0.5, a = 1, b = 0, print_var = c("p", "a", "b"))

p4 = traceplot_mean(p = -1.1, a = 0, b = 2, mu = 0, print_var = c("p", "a", "b")) 
p5 = traceplot_mean(p = -1.5, a = 0, b = 2, mu = 0, print_var = c("p", "a", "b"))
p6 = traceplot_mean(p = -2, a = 0, b = 2, mu = 0, print_var = c("p", "a", "b"))

traceplot_ntc <- p1 / p2 / p3 / p4 / p5 / p6 + plot_layout(ncol=3)
ggsave(traceplot_ntc, filename = "Figures/traceplot_ntc.png", width = 10, height = 5)
