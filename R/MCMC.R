# MCMC sampling
library(ngme2)
library(ggplot2)
library(tidyverse)

### mean, and var of GIG distribution
mean_GIG <- function(p, a, b) {
  m = sqrt(b)/sqrt(a) * besselK(sqrt(a*b), p+1) / besselK(sqrt(a*b), p)
  if (b == 0) {
    # GIG(p, a, 0) = Gamma(p, a/2)
    m = p / (a/2)
  }
  if (a == 0) {
    # GIG(p, 0, b) = InverseGamma(-p, b/2)
    stopifnot(-p > 1)
    m = (b/2) / (-p-1)
  }
  m
}

var_GIG <- function(p, a, b) {
  b/a * (besselK(sqrt(a*b), p+2) / besselK(sqrt(a*b), p) -
    (besselK(sqrt(a*b), p+1) / besselK(sqrt(a*b), p))**2)
}


# Run MCMC (assume K invertible, denote M=KW)
# return samples of V
#
# M|V ~ N(mu*(V-h), sigma*sqrt(V))
# V|M ~ GIG(p-0.5, a+mu^2, (M+mu*h)^2)
MCMC <- function(
  V0, p, a, b, mu, h=1, sigma=1, niter=10000, nburnin=5000, threshold=0.01
){

  if (a == 0 && mu == 0) {
    return (mcmc_invGamma(alpha=-p, beta=b/2, V0=V0, niter=niter, nburnin=nburnin))
  }

  n_chain <- length(V0)
  M <- matrix(0, nrow = niter, ncol = n_chain)
  V <- matrix(0, nrow = niter, ncol = n_chain)
  V[1,] <- V0

  for (iter in 2 : niter){
    M[iter,] <- rnorm(
      n_chain,
      mean = mu * (V[iter-1, ] - h),
      sd = sigma * sqrt(V[iter-1, ])
    )

    if (any(is.na(M[iter, ]))) {
      browser()
V[iter-1, ]
      # NA because V is NA
      print(mu * (V[iter-1, ] - h))

      # M = -mu
      print(M[iter -1,])

      # M = -1, since last V is close to 0
      V[iter-2, ] - h

      # V is close to 0 because b is close to 0
      # b is 0, then V ~ Gamma distribution(p, 2a)
      (M[iter-2, ] + mu*h)**2
    }

# if (b + (M[iter, ] + mu*h)**2 < 1e-6) browser()
    
    V[iter,] <- rgig(
      n_chain,
      p = p-0.5,
      a = a + (mu/sigma)**2,
      b = b + (M[iter, ] + mu*h)**2
    )
  }

  # which_lag acf is less than threshold
  lag = which(acf(V[,1], lag.max = 1000, plot=FALSE)$acf < threshold)[1];

  # Remove the burnins
  V = V[-seq_len(nburnin), ,drop = FALSE]

  V_ind <- V[seq(1, nrow(V), lag), ,drop = FALSE]

  list(
    # M = M[-seq_len(nburnin), ,drop = FALSE],
    V = V,
    V_ind = V_ind
  )
}


### compare samples from MC and MCMC
compare_MC_MCMC <- function(
  p, a=1, b=1, mu=1, 
  h=1, sigma=1,
  niter = 10^6, nburnin=5000,
  xlim = c(0, 5),
  print_var="b"
) {

  # hist(rgig(5000, p, a, b), breaks = 50, freq=F,
  #   main = paste0("p=", p, ", a=", a, ", b=", b))
  # curve(dgig(x, p, a, b), add = T, col = "red")

  # MCMC
  V_samples_ind <- MCMC(
    V0 = 1,
    p = p,
    a = a,
    b = b,
    mu = mu,
    h  = h,
    sigma = sigma,
    niter = niter,
    nburnin = nburnin
  )[["V_ind"]]

  V_filtered <- V_samples_ind[V_samples_ind < xlim[2]]

  # hist(V_filtered, 
  #   xlim = xlim,
  #   breaks = 100, 
  #   freq=F, 
  #   main = paste0("p=", p, ", a=", a, ", b=", b, ", mu=", mu, ", h=", h, ", sigma=", sigma)
  # )
  
  # # add the curve of density of GIG(p, a, b)
  # curve(
  #   dgig(x, p, a, b) / pgig(xlim[2], p, a, b), 
  #   add = T, col = "red"
  # )

  cond_dgig = function(x, p, a, b) {
    dgig(x, p, a, b) / pgig(xlim[2], p, a, b)
  }

  title = ""
  for (var in print_var) {
    title = paste0(title, var, " = ", get(var), ", ")
  }
  title = title %>% substr(1, nchar(title) - 2)

  g = ggplot(data.frame(V=V_filtered), aes(x=V)) +
    geom_histogram(aes(y=..density..), bins=100, color="black", fill="#525150") +
    stat_function(fun=cond_dgig, args=list(p=p, a=a, b=b), color="red") +
    xlim(xlim) +
    # coord_cartesian(xlim = xlim) +
    labs(
      title = title,
      x = NULL,
      y = NULL
    )

  # ll = list(
  #   lag = lag,
  #   n_ind = length(V_samples_ind),
  #   V_samples_ind = V_samples_ind
  # )
  g
}

# InvGamma prior
mcmc_invGamma <- function(
  alpha = 3, 
  beta = 1,
  V0 = 1,
  niter=10^5, 
  nburnin=10^4
) {
  n_chain <- length(V0)
  M <- matrix(0, nrow = niter, ncol = n_chain)
  V <- matrix(0, nrow = niter, ncol = n_chain)
  V[1,] <- V0

  for (iter in 2 : niter){
    M[iter,] <- rnorm(
      n_chain,
      mean = 0,
      sd = sqrt(V[iter-1, ])
    )

    if (any(is.na(M[iter, ]))) {
      browser()
    }

    V[iter,] <- MCMCpack::rinvgamma(
      n_chain,
      shape = alpha+0.5,
      scale = beta + (M[iter, ])**2/2
    )
  }

  # find the lag
  lag = which(acf(V[,1], plot=FALSE, lag.max = 100)$acf < 0.05)[1];

  # keep every lag-th sample
  V_ind <- V[seq(1, nrow(V), lag), ,drop = FALSE]

  if (nburnin > 0)
    return(
      list(M = M[-seq_len(nburnin), ,drop = FALSE],
        V = V[-seq_len(nburnin), ,drop = FALSE],
        V_ind = V_ind))
  else
    return(
      list(M = M, V = V, V_ind = V_ind))
}

compare_mcmc_invGamma <- function(
  alpha, beta
) {
  ret_invGamma <- mcmc_invGamma(alpha, beta)

  # hist(ret_invGamma[["V_ind"]], breaks = 100, freq = FALSE)

  # filter V < 10
  V_filtered <- ret_invGamma[["V"]][ret_invGamma[["V"]] < 10]

  # hist(
  #   V_filtered, breaks = 100, freq = FALSE,
  #   main = paste("alpha = ", alpha, ", beta = ", beta)
  # )
  # curve(
  #   dinvgamma(x, alpha, beta) / ngme2::pigam(10, alpha, beta), 
  #   col = "red",
  #   from = 0, to = 10, add=TRUE
  # )

  cond_dinvgamma = function(x, alpha, beta) {
    MCMCpack::dinvgamma(x, alpha, beta) / ngme2::pigam(10, alpha, beta)
  }

  p = ggplot(data.frame(V=V_filtered), aes(x=V)) +
    geom_histogram(aes(y=..density..), bins=100, color="black", fill="#525150") +
    stat_function(fun=cond_dinvgamma, args=list(alpha=alpha, beta=beta), color="red") +
    xlim(c(0, 10)) +
    labs(
      # title = paste0("alpha = ", alpha, ", beta = ", beta),
      title = paste0("p = ", -alpha, ", a = ", 0, ", b = ", 2*beta),
      x = NULL,
      y = NULL
    )
  p
}



# return a table of sample mean of each chain
experiment_diff_p <- function(
  p_list,
  burnin=10000,
  niter=20000,
  n_chain=5,
  a=1, b=1, mu=5, h=1, sigma=1
) {
  sim_mean_V <- double(length(p_list)); names(sim_mean_V) <- p_list
  mean_V <- double(length(p_list)); names(mean_V) <- p_list
  var_V <- double(length(p_list)); names(var_V) <- p_list
  sim_var_V <- double(length(p_list)); names(sim_var_V) <- p_list

  for (p in p_list){
    V_samples <- MCMC(
      V0 = c(10, 30, 50, 80, 100),
      p = p,
      a = a,
      b = b,
      mu = mu,
      h = h,
      sigma = sigma,
      niter = niter
    )[["V"]]
    df <- as.data.frame(V_samples)
    df[["iter"]] <- 1:niter
    df <- pivot_longer(
      df,
      cols = starts_with("V"),
      names_to = "chain",
      values_to = "V"
    )

    df_filtered <- df %>% filter(iter > niter-1000)
    mean_V[[as.character(p)]] <-
      sqrt(b)/sqrt(a) * besselK(sqrt(a*b), p+1) / besselK(sqrt(a*b), p)
    sim_mean_V[[as.character(p)]] <-
      df_filtered %>% summarize(mean(V)) %>% pull()

    var_V[[as.character(p)]] <- var_GIG(p, a, b)
    sim_var_V[[as.character(p)]] <-
      df_filtered %>% summarize(var(V)) %>% pull()
  }

  return(data.frame(
    mcmc_mean=sim_mean_V,
    true_mean=mean_V,
    rela_error_mean = abs((sim_mean_V - mean_V) / mean_V)
    # mcmc_var=sim_var_V,
    # true_var=var_V,
    # rela_error_var=abs((sim_var_V - var_V) / var_V)
  ))
}

traceplot_V <- function(
  V0 = c(10, 30, 50, 80, 100),
  p = 0.2, a = 1, b = 1, mu=1, h=1, sigma=1,
  niter = 10000,
  nburnin = 1000,
  plot_last = 1000,
  log_scale = FALSE
) {
  V_samples <- MCMC(
    V0 = V0,
    p = p,
    a = a,
    b = b,
    mu = mu,
    h = h,
    sigma = sigma,
    niter = niter,
    nburnin = nburnin
  )[["V"]]
  df <- as.data.frame(V_samples)
  df[["iter"]] <- (nburnin+1) : niter
  df <- pivot_longer(
    df,
    cols = starts_with("V"),
    names_to = "chain",
    values_to = "V"
  )
  mean_V <- mean_GIG(p, a, b)

  # get last Plot_last samples
  df_filtered <- df %>% filter(iter > niter - plot_last)
  transform <- if (log_scale) log else identity

  ggplot(df_filtered, aes(x = iter, y = transform(V), color = chain)) +
    geom_line() +
    geom_hline(yintercept = transform(mean_V),
      linetype = "dashed") +
    geom_hline(yintercept = transform(mean(df_filtered$V)),
      linetype = "dotted") +
    theme_minimal() +
    # show every xlabel every 10 iterations
    # scale_x_continuous(breaks = seq(0, niter - nburnin, by = 10)) +
    labs(
      title = paste0("p=", p, ", a=", a, ", b=", b, ", mu=", mu, ", h=", h, ", sigma=", sigma),
      x = "Iteration",
      y = "V"
    )
}


traceplot_mean <- function(
  V0 = c(1),
  p = 0.2, a = 1, b = 1, mu=1, h=1, sigma=1,
  niter = 5000,
  nburnin = 1000,
  log_scale = FALSE,
  print_var = "b"
) {
  V_samples <- MCMC(
    V0 = V0,
    p = p,
    a = a,
    b = b,
    mu = mu,
    h = h,
    sigma = sigma,
    niter = niter,
    nburnin = nburnin
  )[["V"]][, 1]

  mean_V <- mean_GIG(p, a, b)
  df <- data.frame(V = V_samples)
  df[["iter"]] <- (nburnin+1) : niter

  # accumulated mean of sample mean
  df[["mean_V"]] <- cumsum(df[["V"]]) / (df[["iter"]] - nburnin)
  
  # sd of iterated samples
  df[["sd_V"]] <- sapply(1:nrow(df), function(i) {
    sd(V_samples[1:i])
  })

  transform <- if (log_scale) log else identity

  title = ""
  for (var in print_var) {
    title = paste0(title, var, " = ", get(var), ", ")
  }
  title = title %>% substr(1, nchar(title) - 2)

  # add legend, sample mean and true mean
  ggplot(df, aes(x = iter, y = transform(mean_V))) +
    geom_line() +
    geom_hline(yintercept = transform(mean_V),
      linetype = "dashed") +
    # geom_line(aes(y = transform(mean_V + sd_V)), linetype = "dotted") +
    # geom_line(aes(y = transform(mean_V - sd_V)), linetype = "dotted") +
    theme_minimal() + 
    labs(
      title = title,
      x = "Samples",
      y = NULL
    )
}
  