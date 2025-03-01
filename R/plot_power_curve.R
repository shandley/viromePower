#' Plot Power Curve for Virome Studies
#'
#' Generates a power curve plot showing the relationship between sample size and
#' statistical power for detecting differences in viral abundances.
#'
#' @param effect_size Expected effect size (fold change)
#' @param n_viruses Number of viral taxa tested
#' @param sample_sizes Vector of sample sizes to evaluate
#' @param alpha Significance level (default: 0.05)
#' @param sparsity Proportion of zeros in the data (default: 0.8)
#' @param dispersion Dispersion parameter for viral abundance (default: 2)
#' @param method Test method to use: "wilcoxon", "t.test", or "deseq" (default: "wilcoxon")
#'
#' @return A ggplot object showing the power curve
#' @export
#'
#' @examples
#' power_curve <- plot_power_curve(effect_size = 1.5, n_viruses = 100, 
#'                                sample_sizes = seq(5, 30, by = 5))
plot_power_curve <- function(effect_size, n_viruses, 
                           sample_sizes = seq(5, 30, by = 5),
                           alpha = 0.05, sparsity = 0.8, dispersion = 2,
                           method = "wilcoxon") {
  # Check that ggplot2 is installed
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The 'ggplot2' package is required for this function")
  }
  
  # Validate input parameters
  if (length(sample_sizes) < 2) {
    stop("At least 2 sample sizes are needed to create a power curve")
  }
  if (effect_size <= 0) {
    stop("effect_size must be positive")
  }
  
  # Sort sample sizes
  sample_sizes <- sort(sample_sizes)
  
  # Calculate power for each sample size with error handling
  powers <- numeric(length(sample_sizes))
  
  # For very small sample sizes, use more simulations to improve accuracy
  n_sim_base <- 20
  
  for (i in seq_along(sample_sizes)) {
    # Calculate appropriate number of simulations based on sample size
    # More simulations for smaller sample sizes which tend to be more variable
    n_sim_factor <- max(0.5, min(2, 5 / sample_sizes[i]))
    n_sim_reduced <- max(5, min(30, round(n_sim_base * n_sim_factor)))
    
    # Catch and handle errors for each sample size
    result <- try({
      # Calculate power for this sample size
      calc_virome_power(
        n_samples = sample_sizes[i],
        effect_size = effect_size,
        n_viruses = n_viruses,
        alpha = alpha,
        sparsity = sparsity,
        dispersion = dispersion,
        method = method,
        n_sim = n_sim_reduced
      )
    }, silent = TRUE)
    
    # If calculation succeeded, use the power value, otherwise use estimated value
    if (!inherits(result, "try-error")) {
      powers[i] <- result$power
    } else {
      # If calculation failed, estimate power based on sample size and effect size
      warning("Power calculation failed for sample size ", sample_sizes[i], 
              ". Using estimated power instead.")
      
      # Simple logistic growth model for power estimation
      k <- 20 / effect_size  # Scale factor based on effect size
      powers[i] <- 1 / (1 + exp(-0.3 * (sample_sizes[i] - k)))
    }
  }
  
  # Smooth the power curve to reduce noise
  if (length(sample_sizes) >= 5) {
    # Simple smoothing
    smooth_powers <- powers
    for (i in 2:(length(powers)-1)) {
      smooth_powers[i] <- mean(powers[(i-1):(i+1)])
    }
    powers <- smooth_powers
  }
  
  # Create data frame for plotting
  plot_data <- data.frame(
    sample_size = sample_sizes,
    power = powers
  )
  
  # Find the sample size needed for 80% power (common threshold)
  power_80 <- 0.8
  sample_size_80 <- NA
  
  if (max(powers) >= power_80) {
    if (min(powers) <= power_80) {
      # Interpolate to find the sample size for 80% power
      sample_size_80 <- stats::approx(
        x = powers,
        y = sample_sizes,
        xout = power_80
      )$y
    } else {
      # All powers are above 80%, use the smallest sample size
      sample_size_80 <- min(sample_sizes)
    }
  }
  
  # Create the plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = sample_size, y = power)) +
    ggplot2::geom_line(linewidth = 1.2, color = "#3366cc") +
    ggplot2::geom_point(size = 3, color = "#3366cc") +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "#cc3366") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
    ggplot2::labs(
      x = "Sample Size per Group",
      y = "Statistical Power",
      title = "Power Curve for Virome Study",
      subtitle = paste("Effect size:", effect_size, 
                      "| Number of viruses:", n_viruses, 
                      "| Method:", method)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10)
    )
  
  # Add a vertical line at the 80% power point if applicable
  if (!is.na(sample_size_80)) {
    p <- p + ggplot2::geom_vline(
      xintercept = sample_size_80, 
      linetype = "dashed", 
      color = "#cc3366"
    ) +
    ggplot2::annotate(
      "text",
      x = sample_size_80 + 0.5,
      y = 0.5,
      label = paste("n =", round(sample_size_80, 1)),
      hjust = 0,
      color = "#cc3366"
    )
  }
  
  # Add power and sample size info as a data attribute
  attr(p, "power_data") <- plot_data
  attr(p, "sample_size_80") <- sample_size_80
  
  return(p)
}
