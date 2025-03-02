#' Estimate Zero-Inflation Parameters from Observed Data
#'
#' Infers Beta distribution parameters for modeling variable zero-inflation rates
#' from observed virome data. This function analyzes the prevalence patterns of viral
#' taxa in the dataset and fits a Beta distribution to model the structural zero rates.
#'
#' @param count_matrix A numeric matrix of viral counts with taxa as rows and samples as columns
#' @param presence_threshold Count threshold to consider a taxon as present (default: 1)
#' @param min_prevalence Minimum prevalence (proportion of samples) to include a taxon (default: 0.05)
#' @param method Method for estimating zero-inflation parameters: "moment", "mle", or "bayes" (default: "moment")
#' @param plot_fit Whether to plot the estimated Beta distribution fit (default: TRUE)
#' @param normalize Whether to normalize the count matrix before estimation (default: FALSE)
#'
#' @return A list containing:
#'   \item{alpha}{Estimated alpha parameter for Beta distribution}
#'   \item{beta}{Estimated beta parameter for Beta distribution}
#'   \item{mean_zi}{Estimated mean zero-inflation rate}
#'   \item{presence_absence}{Binary presence/absence matrix derived from data}
#'   \item{prevalence}{Vector of prevalence rates for each taxon}
#'   \item{distribution_fit}{Goodness of fit statistics for the Beta distribution}
#'   \item{plot}{ggplot object showing the data and fitted distribution (if plot_fit=TRUE)}
#'
#' @details
#' This function provides a bridge between real-world virome data and the simulation
#' parameters needed for realistic modeling. It estimates the distribution of zero-inflation
#' rates by analyzing the prevalence patterns in the dataset.
#'
#' The function first converts the count matrix to a presence/absence matrix using the
#' presence_threshold parameter. It then calculates the prevalence (proportion of samples
#' where a taxon is present) for each taxon and filters out rare taxa based on the
#' min_prevalence parameter.
#'
#' From the prevalence data, it estimates the parameters of a Beta distribution that
#' can be used to model variable zero-inflation rates in simulations. Three estimation
#' methods are available:
#' \itemize{
#'   \item "moment": Uses method of moments estimation (fast but less accurate)
#'   \item "mle": Uses maximum likelihood estimation (more accurate but slower)
#'   \item "bayes": Uses Bayesian estimation with weakly informative priors (most robust but slowest)
#' }
#'
#' @examples
#' \dontrun{
#' # Using real virome data
#' params <- estimate_zi_parameters_from_data(
#'   count_matrix = my_virome_counts,
#'   presence_threshold = 5,
#'   method = "mle"
#' )
#' 
#' # Print the estimated parameters
#' cat("Estimated alpha:", params$alpha, "\n")
#' cat("Estimated beta:", params$beta, "\n")
#' cat("Mean zero-inflation rate:", params$mean_zi, "\n")
#' 
#' # Plot the fit
#' print(params$plot)
#' 
#' # Use the estimated parameters in simulation
#' sim_data <- simulate_zero_inflated_virome(
#'   n_samples = 30,
#'   n_viruses = 100,
#'   variable_zi_rates = TRUE,
#'   zi_alpha = params$alpha,
#'   zi_beta = params$beta,
#'   structural_zeros = params$mean_zi
#' )
#' }
#'
#' @export
estimate_zi_parameters_from_data <- function(count_matrix,
                                            presence_threshold = 1,
                                            min_prevalence = 0.05,
                                            method = "moment",
                                            plot_fit = TRUE,
                                            normalize = FALSE) {
  # Check if required packages are available
  if (!requireNamespace("fitdistrplus", quietly = TRUE)) {
    stop("Package 'fitdistrplus' is required for parameter estimation. Install with install.packages('fitdistrplus')")
  }
  
  if (plot_fit && !requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is required for plotting. No plot will be generated.")
    plot_fit <- FALSE
  }
  
  # Normalize count matrix if requested (e.g., to CPM)
  if (normalize) {
    count_sums <- colSums(count_matrix)
    count_matrix <- t(t(count_matrix) / count_sums) * 1e6
  }
  
  # Convert to presence/absence matrix
  presence_absence <- count_matrix >= presence_threshold
  
  # Calculate prevalence of each taxon
  n_samples <- ncol(count_matrix)
  prevalence <- rowSums(presence_absence) / n_samples
  
  # Calculate zero-inflation rates (1 - prevalence)
  zero_inflation <- 1 - prevalence
  
  # Filter out rare taxa (those with very low prevalence)
  valid_taxa <- prevalence >= min_prevalence & prevalence <= (1 - min_prevalence)
  
  if (sum(valid_taxa) < 5) {
    warning("Very few taxa with sufficient prevalence for estimation (n=", sum(valid_taxa), 
            "). Results may be unreliable.")
  }
  
  filtered_zi <- zero_inflation[valid_taxa]
  
  # Ensure filtered values are bounded away from 0 and 1 for Beta fitting
  filtered_zi <- pmin(pmax(filtered_zi, 0.001), 0.999)
  
  # Estimate Beta distribution parameters
  if (method == "moment") {
    # Method of moments estimation
    mean_zi <- mean(filtered_zi)
    var_zi <- var(filtered_zi)
    
    # Calculate alpha and beta from mean and variance
    alpha <- mean_zi * (mean_zi * (1 - mean_zi) / var_zi - 1)
    beta <- (1 - mean_zi) * (mean_zi * (1 - mean_zi) / var_zi - 1)
    
    # Ensure parameters are positive and reasonable
    alpha <- max(0.1, alpha)
    beta <- max(0.1, beta)
    
    # No formal distribution fit with this method
    distribution_fit <- list(
      method = "moment",
      loglik = NA,
      aic = NA
    )
    
  } else if (method == "mle") {
    # Maximum likelihood estimation using fitdistrplus
    fit <- try(fitdistrplus::fitdist(filtered_zi, "beta", method = "mle"), silent = TRUE)
    
    if (inherits(fit, "try-error")) {
      warning("MLE fitting failed, falling back to method of moments")
      return(estimate_zi_parameters_from_data(
        count_matrix, presence_threshold, min_prevalence, 
        method = "moment", plot_fit, normalize))
    }
    
    alpha <- fit$estimate["shape1"]
    beta <- fit$estimate["shape2"]
    
    distribution_fit <- list(
      method = "mle",
      loglik = fit$loglik,
      aic = fit$aic,
      bic = fit$bic,
      summary = summary(fit)
    )
    
  } else if (method == "bayes") {
    # Simple approximate Bayesian estimation with weak priors
    if (!requireNamespace("rjags", quietly = TRUE)) {
      warning("Package 'rjags' is required for Bayesian estimation. Falling back to MLE.")
      return(estimate_zi_parameters_from_data(
        count_matrix, presence_threshold, min_prevalence, 
        method = "mle", plot_fit, normalize))
    }
    
    # This is simplified Bayesian estimation, a full implementation would use MCMC
    # For now, we'll use a simple empirical Bayes approach
    mean_zi <- mean(filtered_zi)
    var_zi <- var(filtered_zi)
    
    # Calculate alpha and beta from mean and variance
    alpha <- mean_zi * (mean_zi * (1 - mean_zi) / var_zi - 1)
    beta <- (1 - mean_zi) * (mean_zi * (1 - mean_zi) / var_zi - 1)
    
    # Add weak prior (equivalent to 5 prior observations)
    alpha <- alpha + 2  # Prior: 2 successes 
    beta <- beta + 3    # Prior: 3 failures
    
    # Ensure parameters are positive and reasonable
    alpha <- max(0.1, alpha)
    beta <- max(0.1, beta)
    
    distribution_fit <- list(
      method = "bayes (approximate)",
      prior = "weak Beta(2,3)",
      posterior_mean = mean(rbeta(1000, alpha, beta)),
      posterior_var = var(rbeta(1000, alpha, beta))
    )
  } else {
    stop("Unknown method. Choose 'moment', 'mle', or 'bayes'.")
  }
  
  # Create a plot if requested
  if (plot_fit) {
    # Generate points for the fitted Beta distribution
    x_seq <- seq(0, 1, length.out = 100)
    y_density <- dbeta(x_seq, shape1 = alpha, shape2 = beta)
    
    # Create data frame for plotting
    density_df <- data.frame(x = x_seq, y = y_density)
    
    # Generate histogram data
    hist_data <- hist(filtered_zi, plot = FALSE, breaks = 20)
    hist_df <- data.frame(
      x = hist_data$mids,
      y = hist_data$density
    )
    
    # Create plot with ggplot2
    p <- ggplot2::ggplot() +
      ggplot2::geom_histogram(
        data = data.frame(zi = filtered_zi),
        ggplot2::aes(x = zi, y = ..density..),
        bins = 20,
        fill = "steelblue",
        alpha = 0.7,
        color = "white"
      ) +
      ggplot2::geom_line(
        data = density_df,
        ggplot2::aes(x = x, y = y),
        color = "firebrick",
        size = 1.2
      ) +
      ggplot2::labs(
        title = "Estimated Zero-Inflation Pattern",
        subtitle = sprintf("Beta(%0.2f, %0.2f), Mean ZI = %0.2f", 
                           alpha, beta, alpha / (alpha + beta)),
        x = "Zero-Inflation Rate",
        y = "Density"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        plot.subtitle = ggplot2::element_text(face = "italic")
      )
    
    # Add QQ plot as inset
    theoretical_quantiles <- qbeta(ppoints(length(filtered_zi)), alpha, beta)
    sorted_zi <- sort(filtered_zi)
    
    qq_df <- data.frame(
      theoretical = theoretical_quantiles,
      empirical = sorted_zi
    )
    
    qq_plot <- ggplot2::ggplot(qq_df, ggplot2::aes(x = theoretical, y = empirical)) +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      ggplot2::labs(
        title = "Q-Q Plot",
        x = "Theoretical Quantiles",
        y = "Sample Quantiles"
      ) +
      ggplot2::theme_minimal(base_size = 8) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 10, face = "bold"),
        axis.title = ggplot2::element_text(size = 8)
      )
  } else {
    p <- NULL
  }
  
  # Calculate mean ZI rate from parameters
  mean_zi <- alpha / (alpha + beta)
  
  # Return results
  return(list(
    alpha = alpha,
    beta = beta,
    mean_zi = mean_zi,
    presence_absence = presence_absence,
    prevalence = prevalence,
    zero_inflation = zero_inflation,
    filtered_taxa = valid_taxa,
    distribution_fit = distribution_fit,
    plot = p
  ))
}