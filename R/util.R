

traceplot <- function(
  est, 
  rho, 
  mu, 
  sigma, 
  nu, 
  sigma_e,
  main = NULL
) {
  compute_ylim <- function(x) {
    ylow <- if (min(x) > 0) 0.9 * min(x) else 1.1 * min(x)
    yhigh <- if (max(x) > 0) 1.1 * max(x) else 0.9 * max(x)
    c(ylow, yhigh)
  }

  par(mfrow = c(3, 2))
  plot(est$mu, type = "l", main = "mu", xlab = "Iteration", ylab = "Value", 
    ylim = compute_ylim(est$mu))
  abline(h = mu, col = "red")
  
  plot(est$sigma, type = "l", main = "sigma", xlab = "Iteration", ylab = "Value", 
    ylim = compute_ylim(est$sigma))
  abline(h = sigma, col = "red")
  
  plot(est$nu, type = "l", main = "nu", xlab = "Iteration", ylab = "Value", 
    ylim = compute_ylim(est$nu))
  abline(h = nu, col = "red")

  plot(est$sigma_e, type = "l", main = "sigma_e", xlab = "Iteration", ylab = "Value", 
    ylim = compute_ylim(est$sigma_e))
  abline(h = sigma_e, col = "red")
  
  plot(est$rho, type = "l", main = "rho", xlab = "Iteration", ylab = "Value", 
    ylim = compute_ylim(est$rho))
  abline(h = rho, col = "red")

  # Add a single title for the entire plot
  if (!is.null(main)) {
    mtext(main, side = 3, outer = TRUE, line = -2, cex = 1.5)
  }
}


traceplot_ggplot <- function(est, rho, mu, sigma, nu, sigma_e, main = NULL) {
  library(ggplot2)
  library(patchwork)

  # Create data frames for each parameter
  df_mu <- data.frame(iteration = seq_along(est$mu), value = est$mu)
  df_sigma <- data.frame(iteration = seq_along(est$sigma), value = est$sigma)
  df_nu <- data.frame(iteration = seq_along(est$nu), value = est$nu)
  df_sigma_e <- data.frame(iteration = seq_along(est$sigma_e), value = est$sigma_e)
  df_rho <- data.frame(iteration = seq_along(est$rho), value = est$rho)

  # Create individual plots
  plot_mu <- ggplot(df_mu, aes(x = iteration, y = value)) +
    geom_line() +
    geom_hline(yintercept = mu, color = "red") +
    ggtitle("mu") +
    theme_minimal()

  plot_sigma <- ggplot(df_sigma, aes(x = iteration, y = value)) +
    geom_line() +
    geom_hline(yintercept = sigma, color = "red") +
    ggtitle("sigma") +
    theme_minimal()

  plot_nu <- ggplot(df_nu, aes(x = iteration, y = value)) +
    geom_line() +
    geom_hline(yintercept = nu, color = "red") +
    ggtitle("nu") +
    theme_minimal()

  plot_sigma_e <- ggplot(df_sigma_e, aes(x = iteration, y = value)) +
    geom_line() +
    geom_hline(yintercept = sigma_e, color = "red") +
    ggtitle("sigma_e") +
    theme_minimal()

  plot_rho <- ggplot(df_rho, aes(x = iteration, y = value)) +
    geom_line() +
    geom_hline(yintercept = rho, color = "red") +
    ggtitle("rho") +
    theme_minimal()

  # Combine the plots into a 3x2 grid and add a single title
  combined_plot <- (plot_mu | plot_sigma) / (plot_nu | plot_sigma_e) / (plot_rho | plot_spacer())

  if (!is.null(main)) {
    combined_plot <- combined_plot + plot_annotation(title = main)
  }

  # Print the combined plot
  print(combined_plot)
}



rho_to_th <- function(r) {
  log((-1 - r) / (-1 + r))
}

th_to_rho <- function(th) {
  (exp(th) - 1) / (exp(th) + 1)
}

convert_noise <- function(noise, nu, h) {
  n <- length(h)
  if (noise == "gal") {
    p <- h * nu
    a <- rep(2 * nu, n)
    b <- rep(0, n)
  } else if (noise == "nig") {
    p <- rep(-0.5, n)
    a <- rep(nu, n)
    b <- nu * h^2
  }

  return(list(p = p, a = a, b = b))
}


traceplot_df <- function(
  rho_df, mu_df, sigma_df, nu_df, sigma_e_df, 
  true_rho, true_mu, true_sigma, true_nu, true_sigma_e,
  main = NULL
) {
  library(ggplot2)
  library(patchwork)

  plot_one_df <- function(df, true_value, title, value, legend=TRUE) {
    iteration <- nrow(df)
    if (legend) {
      plot <- ggplot(df, aes(x = iteration, y = !!sym(value), color = chain)) +
        geom_line() +
        geom_hline(aes(yintercept = true_value, linetype = "True Value"), color = "black") +
        ggtitle(title) +
        theme_minimal() +
        scale_color_discrete(name = "Chain") +
        scale_linetype_manual(name = "Reference", values = c("True Value" = "dashed")) +
        guides(color = guide_legend(order = 1),
               linetype = guide_legend(order = 2)) +
        theme(legend.position = "bottom",
              legend.box = "vertical")
    } else {
      plot <- ggplot(df, aes(x = iteration, y = !!sym(value), color = chain)) +
        geom_line() +
        geom_hline(yintercept = true_value, linetype = "dashed", color = "black") +
        ggtitle(title) +
        theme_minimal() +
        theme(legend.position = "none",
              axis.title.y = element_blank())
    }
    return(plot)
  }

  plot_rho <- plot_one_df(rho_df, true_rho, "rho", "rho", legend=FALSE)
  plot_mu <- plot_one_df(mu_df, true_mu, "mu", "mu", legend=FALSE)
  plot_sigma <- plot_one_df(sigma_df, true_sigma, "sigma", "sigma", legend=FALSE)
  plot_nu <- plot_one_df(nu_df, true_nu, "nu", "nu", legend=FALSE)
  plot_sigma_e <- plot_one_df(sigma_e_df, true_sigma_e, "sigma_e", "sigma_e", legend=FALSE)

  library(cowplot)
  # Function to extract the legend from a ggplot object
  extract_legend <- function(p) {
    g <- ggplotGrob(p)
    legend <- g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
    return(legend)
  }
  # Extract the legend from plot1
  plot_rho_legend <- plot_one_df(rho_df, true_rho, "rho", "rho", legend=TRUE)
  legend_plot <- extract_legend(plot_rho_legend)
  # Create an empty plot and add the legend
  legend_only_plot <- ggdraw() + draw_grob(legend_plot)

  # Combine the plots into a 3x2 grid and add a single title
  combined_plot <- (plot_rho | plot_mu) / (plot_sigma | plot_nu) / (plot_sigma_e | legend_only_plot) 

  if (!is.null(main)) {
    combined_plot <- combined_plot +
      plot_annotation(
        title = main,
      theme = theme(
        plot.title = element_text(size = 15, hjust = 0.5)  # hjust = 0.5 centers the title
        )
      )
  }

  # Print the combined plot
  return(combined_plot)
}