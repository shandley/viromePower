#' Plot Observed vs Expected Zeros
#'
#' Creates a visualization comparing observed zero counts with expected zero counts
#' under standard and zero-inflated negative binomial models.
#'
#' @param obs_data A matrix or data frame containing the observed count data with
#'   samples as rows and taxa as columns.
#' @param standard_model_params List containing standard model parameters (mu, size).
#' @param zinb_model_params List containing zero-inflated model parameters (mu, size, zi_prob).
#' @param group_factor Optional factor indicating group membership for samples.
#'   If provided, points will be colored by group.
#' @param title Optional custom title for the plot.
#'
#' @return A ggplot2 object displaying the comparison between observed and expected zeros.
#'
#' @examples
#' \dontrun{
#' # Simulate zero-inflated data
#' sim_data <- simulate_zero_inflated_virome(
#'   n_samples = 50,
#'   n_taxa = 100,
#'   mu = 10,
#'   size = 0.5,
#'   zi_prob = 0.3,
#'   group_effect = 1.5,
#'   group_factor = rep(c("A", "B"), each = 25)
#' )
#'
#' # Plot observed vs expected zeros
#' plot_observed_vs_expected_zeros(
#'   obs_data = sim_data$counts,
#'   standard_model_params = list(mu = 10, size = 0.5),
#'   zinb_model_params = list(mu = 10, size = 0.5, zi_prob = 0.3),
#'   group_factor = sim_data$groups
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline theme_minimal labs
#' @importFrom stats rnbinom
#' @export
plot_observed_vs_expected_zeros <- function(obs_data, standard_model_params, zinb_model_params, 
                                            group_factor = NULL, title = "Observed vs. Expected Zeros") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.", call. = FALSE)
  }
  
  # Validate inputs
  if (!is.matrix(obs_data) && !is.data.frame(obs_data)) {
    stop("obs_data must be a matrix or data frame")
  }
  
  if (!is.list(standard_model_params) || !all(c("mu", "size") %in% names(standard_model_params))) {
    stop("standard_model_params must be a list containing 'mu' and 'size'")
  }
  
  if (!is.list(zinb_model_params) || !all(c("mu", "size", "zi_prob") %in% names(zinb_model_params))) {
    stop("zinb_model_params must be a list containing 'mu', 'size', and 'zi_prob'")
  }
  
  # Calculate observed zeros per sample
  observed_zeros <- rowSums(obs_data == 0)
  n_taxa <- ncol(obs_data)
  
  # Calculate expected zeros under standard negative binomial
  p_zero_nb <- (standard_model_params$size / (standard_model_params$size + standard_model_params$mu))^standard_model_params$size
  expected_zeros_nb <- n_taxa * p_zero_nb
  
  # Calculate expected zeros under zero-inflated negative binomial
  p_zero_zinb <- zinb_model_params$zi_prob + (1 - zinb_model_params$zi_prob) * 
    (zinb_model_params$size / (zinb_model_params$size + zinb_model_params$mu))^zinb_model_params$size
  expected_zeros_zinb <- n_taxa * p_zero_zinb
  
  # Create data frame for plotting
  plot_data <- data.frame(
    sample_id = 1:nrow(obs_data),
    observed = observed_zeros,
    expected_nb = expected_zeros_nb,
    expected_zinb = expected_zeros_zinb
  )
  
  # Add group information if provided
  if (!is.null(group_factor)) {
    if (length(group_factor) != nrow(obs_data)) {
      stop("Length of group_factor must match number of samples in obs_data")
    }
    plot_data$group <- group_factor
  }
  
  # Create plots
  p1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = expected_nb, y = observed)) +
    ggplot2::geom_point(ggplot2::aes(color = if (!is.null(group_factor)) group else NULL), 
                        alpha = 0.7, size = 3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Standard Negative Binomial Model",
      x = "Expected Number of Zeros",
      y = "Observed Number of Zeros"
    )
  
  p2 <- ggplot2::ggplot(plot_data, ggplot2::aes(x = expected_zinb, y = observed)) +
    ggplot2::geom_point(ggplot2::aes(color = if (!is.null(group_factor)) group else NULL), 
                        alpha = 0.7, size = 3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Zero-Inflated Negative Binomial Model",
      x = "Expected Number of Zeros",
      y = "Observed Number of Zeros"
    )
  
  # Combine plots if gridExtra is available
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined_plot <- gridExtra::grid.arrange(
      p1, p2, 
      ncol = 2, 
      top = grid::textGrob(title, gp = grid::gpar(fontsize = 14, fontface = "bold"))
    )
    return(combined_plot)
  } else {
    message("For combined plots, install the 'gridExtra' package. Returning standard model plot.")
    return(p1)
  }
}

#' Plot Zero-Inflation Parameter Distribution
#'
#' Creates visualizations of the zero-inflation parameter distribution across samples or taxa.
#'
#' @param zinb_model_fit A fitted zero-inflated model object or a list containing
#'   estimated zero-inflation parameters.
#' @param by_taxa Logical indicating whether to show distribution by taxa (TRUE) or
#'   by samples (FALSE, default).
#' @param group_factor Optional factor indicating group membership for samples or taxa.
#'   If provided, distributions will be faceted by group.
#' @param highlight_threshold Optional numeric value. Zero-inflation parameters
#'   above this threshold will be highlighted.
#'
#' @return A ggplot2 object visualizing the distribution of zero-inflation parameters.
#'
#' @examples
#' \dontrun{
#' # Simulate zero-inflated data
#' sim_data <- simulate_zero_inflated_virome(
#'   n_samples = 50,
#'   n_taxa = 100,
#'   mu = 10,
#'   size = 0.5,
#'   zi_prob = 0.3,
#'   group_effect = 1.5,
#'   group_factor = rep(c("A", "B"), each = 25)
#' )
#'
#' # Assume we've fit a ZINB model and have parameter estimates
#' zi_params <- list(
#'   pi = matrix(rbeta(50*100, 2, 5), nrow = 50),  # Example parameter estimates
#'   groups = sim_data$groups,
#'   taxa_names = paste0("Taxon_", 1:100)
#' )
#'
#' # Plot zero-inflation distribution
#' plot_zero_inflation_distribution(
#'   zinb_model_fit = zi_params,
#'   by_taxa = FALSE,
#'   group_factor = sim_data$groups,
#'   highlight_threshold = 0.5
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_density facet_wrap theme_minimal labs geom_vline
#' @export
plot_zero_inflation_distribution <- function(zinb_model_fit, by_taxa = FALSE, 
                                             group_factor = NULL, highlight_threshold = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.", call. = FALSE)
  }
  
  # Extract zero-inflation parameters from model or list
  if (is.list(zinb_model_fit) && "pi" %in% names(zinb_model_fit)) {
    zi_params <- zinb_model_fit$pi
  } else {
    stop("zinb_model_fit must be a list containing 'pi' (zero-inflation parameters)")
  }
  
  # Convert to data frame for plotting
  if (by_taxa) {
    # Calculate mean ZI parameter for each taxon
    if (is.matrix(zi_params)) {
      zi_values <- colMeans(zi_params)
      taxa_names <- if (!is.null(zinb_model_fit$taxa_names)) zinb_model_fit$taxa_names else paste0("Taxon_", 1:length(zi_values))
      plot_data <- data.frame(
        taxon = taxa_names,
        zi_prob = zi_values
      )
      if (!is.null(group_factor) && length(group_factor) == length(zi_values)) {
        plot_data$group <- group_factor
      }
      x_lab <- "Zero-Inflation Probability"
      y_lab <- "Count"
      title <- "Distribution of Zero-Inflation Parameters by Taxon"
    } else {
      stop("Zero-inflation parameters must be in matrix format for by_taxa=TRUE")
    }
  } else {
    # Calculate mean ZI parameter for each sample
    if (is.matrix(zi_params)) {
      zi_values <- rowMeans(zi_params)
      sample_ids <- 1:length(zi_values)
      plot_data <- data.frame(
        sample = sample_ids,
        zi_prob = zi_values
      )
      if (!is.null(group_factor) && length(group_factor) == length(zi_values)) {
        plot_data$group <- group_factor
      }
      x_lab <- "Zero-Inflation Probability"
      y_lab <- "Count"
      title <- "Distribution of Zero-Inflation Parameters by Sample"
    } else {
      stop("Zero-inflation parameters must be in matrix format for by_sample calculation")
    }
  }
  
  # Create distribution plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = zi_prob)) +
    ggplot2::geom_histogram(ggplot2::aes(fill = if (!is.null(group_factor)) group else NULL),
                           alpha = 0.7, bins = 30, position = "identity") +
    ggplot2::geom_density(ggplot2::aes(color = if (!is.null(group_factor)) group else NULL),
                        alpha = 0.5, linewidth = 1) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = x_lab,
      y = y_lab
    )
  
  # Add highlight threshold if provided
  if (!is.null(highlight_threshold)) {
    p <- p + ggplot2::geom_vline(xintercept = highlight_threshold, 
                               linetype = "dashed", color = "red", linewidth = 1)
  }
  
  # Add faceting by group if provided
  if (!is.null(group_factor) && requireNamespace("dplyr", quietly = TRUE)) {
    p <- p + ggplot2::facet_wrap(~ group, scales = "free_y")
  }
  
  return(p)
}

#' Compare Power Curves Between Standard and Zero-Inflated Models
#'
#' Visualizes and compares statistical power between standard and zero-inflated
#' negative binomial models across different sample sizes or effect sizes.
#'
#' @param standard_power_results A data frame or list containing power analysis results
#'   from standard model calculations.
#' @param zinb_power_results A data frame or list containing power analysis results
#'   from zero-inflated model calculations.
#' @param x_variable Character specifying which variable to use on x-axis: "n_samples"
#'   (default) or "effect_size".
#' @param alpha Significance level used in power calculations (default: 0.05).
#' @param highlight_thresh Optional numeric threshold to highlight on the plot (e.g., 0.8 for 80% power).
#'
#' @return A ggplot2 object comparing power curves between the two models.
#'
#' @examples
#' \dontrun{
#' # Calculate power using standard model
#' std_power <- calc_virome_power(
#'   n_samples = seq(10, 100, by = 10),
#'   effect_size = 1.5,
#'   mu = 10,
#'   size = 0.5
#' )
#'
#' # Calculate power using zero-inflated model
#' zinb_power <- calc_zinb_bayesian_power(
#'   n_samples = seq(10, 100, by = 10),
#'   effect_size = 1.5,
#'   mu = 10,
#'   size = 0.5,
#'   zi_prob = 0.3
#' )
#'
#' # Compare power curves
#' compare_power_curves(
#'   standard_power_results = std_power,
#'   zinb_power_results = zinb_power,
#'   x_variable = "n_samples",
#'   highlight_thresh = 0.8
#' )
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_point geom_hline theme_minimal labs
#' @export
compare_power_curves <- function(standard_power_results, zinb_power_results,
                                 x_variable = "n_samples", alpha = 0.05,
                                 highlight_thresh = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.", call. = FALSE)
  }
  
  # Validate inputs
  valid_x_vars <- c("n_samples", "effect_size")
  if (!x_variable %in% valid_x_vars) {
    stop(paste0("x_variable must be one of: ", paste(valid_x_vars, collapse = ", ")))
  }
  
  # Prepare standard model results
  if (is.list(standard_power_results) && !is.data.frame(standard_power_results)) {
    std_data <- data.frame(
      x_val = standard_power_results[[x_variable]],
      power = standard_power_results$power,
      model = "Standard NB"
    )
  } else if (is.data.frame(standard_power_results)) {
    if (!all(c(x_variable, "power") %in% names(standard_power_results))) {
      stop(paste0("standard_power_results must contain columns: ", x_variable, " and power"))
    }
    std_data <- data.frame(
      x_val = standard_power_results[[x_variable]],
      power = standard_power_results$power,
      model = "Standard NB"
    )
  } else {
    stop("standard_power_results must be a data frame or list containing power analysis results")
  }
  
  # Prepare zero-inflated model results
  if (is.list(zinb_power_results) && !is.data.frame(zinb_power_results)) {
    zinb_data <- data.frame(
      x_val = zinb_power_results[[x_variable]],
      power = zinb_power_results$power,
      model = "Zero-Inflated NB"
    )
  } else if (is.data.frame(zinb_power_results)) {
    if (!all(c(x_variable, "power") %in% names(zinb_power_results))) {
      stop(paste0("zinb_power_results must contain columns: ", x_variable, " and power"))
    }
    zinb_data <- data.frame(
      x_val = zinb_power_results[[x_variable]],
      power = zinb_power_results$power,
      model = "Zero-Inflated NB"
    )
  } else {
    stop("zinb_power_results must be a data frame or list containing power analysis results")
  }
  
  # Combine data
  plot_data <- rbind(std_data, zinb_data)
  
  # Set appropriate axis labels
  x_lab <- if (x_variable == "n_samples") "Number of Samples" else "Effect Size"
  
  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x_val, y = power, color = model, linetype = model)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::labs(
      title = "Power Comparison: Standard vs. Zero-Inflated Models",
      x = x_lab,
      y = "Statistical Power",
      color = "Model Type",
      linetype = "Model Type"
    ) +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2))
  
  # Add significance level line
  p <- p + ggplot2::geom_hline(yintercept = alpha, linetype = "dotted", 
                              color = "darkgray", linewidth = 0.8)
  
  # Add highlight threshold if provided
  if (!is.null(highlight_thresh)) {
    p <- p + ggplot2::geom_hline(yintercept = highlight_thresh, linetype = "dashed", 
                               color = "red", linewidth = 1)
  }
  
  # Add power difference annotation if both models have same x values
  if (all(std_data$x_val == zinb_data$x_val) && requireNamespace("dplyr", quietly = TRUE)) {
    diff_data <- data.frame(
      x_val = std_data$x_val,
      power_diff = zinb_data$power - std_data$power
    )
    
    # Create secondary plot showing difference
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      p2 <- ggplot2::ggplot(diff_data, ggplot2::aes(x = x_val, y = power_diff)) +
        ggplot2::geom_bar(stat = "identity", fill = "steelblue", alpha = 0.7) +
        ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "black") +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = "Power Difference (ZINB - Standard)",
          x = x_lab,
          y = "Power Difference"
        )
      
      combined_plot <- gridExtra::grid.arrange(p, p2, ncol = 1, heights = c(2, 1))
      return(combined_plot)
    }
  }
  
  return(p)
}

#' Plot Zero-Inflated Negative Binomial Model Diagnostics
#'
#' Creates diagnostic plots for assessing the fit of zero-inflated negative binomial models
#' to observed virome data.
#'
#' @param obs_data A matrix or data frame containing the observed count data with
#'   samples as rows and taxa as columns.
#' @param zinb_model_fit A fitted zero-inflated model object or a list containing
#'   model parameters (mu, size, zi_prob).
#' @param n_taxa_to_plot Integer specifying the number of randomly selected taxa to
#'   show in detailed plots (default: 4).
#' @param seed Optional random seed for reproducible taxa selection.
#'
#' @return A list of ggplot2 objects with different diagnostic plots.
#'
#' @examples
#' \dontrun{
#' # Simulate zero-inflated data
#' sim_data <- simulate_zero_inflated_virome(
#'   n_samples = 50,
#'   n_taxa = 100,
#'   mu = 10,
#'   size = 0.5,
#'   zi_prob = 0.3,
#'   group_effect = 1.5,
#'   group_factor = rep(c("A", "B"), each = 25)
#' )
#'
#' # Assume we've fit a ZINB model and have parameter estimates
#' zinb_fit <- list(
#'   mu = matrix(rgamma(50*100, 5, 0.5), nrow = 50),
#'   size = 0.5,
#'   zi_prob = matrix(rbeta(50*100, 2, 5), nrow = 50)
#' )
#'
#' # Create diagnostic plots
#' diagnostics <- plot_zinb_diagnostics(
#'   obs_data = sim_data$counts,
#'   zinb_model_fit = zinb_fit,
#'   n_taxa_to_plot = 6,
#'   seed = 123
#' )
#'
#' # Display the plots
#' diagnostics$zero_proportion_plot
#' diagnostics$qq_plot
#' diagnostics$taxa_fit_plots
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_abline theme_minimal labs facet_wrap
#' @importFrom stats rnbinom quantile
#' @export
plot_zinb_diagnostics <- function(obs_data, zinb_model_fit, n_taxa_to_plot = 4, seed = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is needed for this function to work. Please install it.", call. = FALSE)
  }
  
  # Validate inputs
  if (!is.matrix(obs_data) && !is.data.frame(obs_data)) {
    stop("obs_data must be a matrix or data frame")
  }
  
  if (!is.list(zinb_model_fit)) {
    stop("zinb_model_fit must be a list containing model parameters")
  }
  
  # Extract model parameters
  if (all(c("mu", "size", "zi_prob") %in% names(zinb_model_fit))) {
    mu <- zinb_model_fit$mu
    size <- zinb_model_fit$size
    zi_prob <- zinb_model_fit$zi_prob
  } else {
    stop("zinb_model_fit must contain 'mu', 'size', and 'zi_prob' parameters")
  }
  
  # Ensure compatible dimensions
  if (is.matrix(mu) && nrow(mu) != nrow(obs_data) || ncol(mu) != ncol(obs_data)) {
    stop("Dimensions of mu matrix must match obs_data dimensions")
  }
  
  # If mu is not a matrix, create one with the same value repeated
  if (!is.matrix(mu)) {
    mu <- matrix(mu, nrow = nrow(obs_data), ncol = ncol(obs_data))
  }
  
  # Same for zi_prob
  if (!is.matrix(zi_prob)) {
    zi_prob <- matrix(zi_prob, nrow = nrow(obs_data), ncol = ncol(obs_data))
  }
  
  # Calculate proportion of zeros by taxon
  zero_prop_obs <- colMeans(obs_data == 0)
  
  # Calculate expected proportion of zeros under ZINB model
  calc_exp_zeros <- function(mu_vec, size_val, zi_prob_vec) {
    p_zero_nb <- (size_val / (size_val + mu_vec))^size_val
    p_zero_zinb <- zi_prob_vec + (1 - zi_prob_vec) * p_zero_nb
    return(p_zero_zinb)
  }
  
  zero_prop_exp <- sapply(1:ncol(obs_data), function(j) {
    mean(calc_exp_zeros(mu[, j], size, zi_prob[, j]))
  })
  
  # Create zero proportion comparison plot
  zero_plot_data <- data.frame(
    taxon = 1:ncol(obs_data),
    observed = zero_prop_obs,
    expected = zero_prop_exp
  )
  
  zero_prop_plot <- ggplot2::ggplot(zero_plot_data, ggplot2::aes(x = expected, y = observed)) +
    ggplot2::geom_point(alpha = 0.7, size = 3) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Observed vs. Expected Zero Proportions by Taxon",
      x = "Expected Proportion of Zeros (ZINB Model)",
      y = "Observed Proportion of Zeros"
    )
  
  # Create QQ plot for model fit
  # Generate expected counts under ZINB model
  generate_zinb_counts <- function(n, mu_val, size_val, zi_prob_val) {
    is_zero <- runif(n) < zi_prob_val
    counts <- rnbinom(n, mu = mu_val, size = size_val)
    counts[is_zero] <- 0
    return(counts)
  }
  
  # Collect observed and expected counts
  all_obs_counts <- as.vector(obs_data)
  all_exp_counts <- numeric(length(all_obs_counts))
  
  idx <- 1
  for (i in 1:nrow(obs_data)) {
    for (j in 1:ncol(obs_data)) {
      all_exp_counts[idx] <- generate_zinb_counts(1, mu[i, j], size, zi_prob[i, j])
      idx <- idx + 1
    }
  }
  
  # Create QQ plot data
  qq_data <- data.frame(
    observed_quantiles = quantile(all_obs_counts, probs = seq(0, 1, 0.01)),
    expected_quantiles = quantile(all_exp_counts, probs = seq(0, 1, 0.01))
  )
  
  qq_plot <- ggplot2::ggplot(qq_data, ggplot2::aes(x = expected_quantiles, y = observed_quantiles)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Quantile-Quantile Plot for ZINB Model Fit",
      x = "Expected Quantiles",
      y = "Observed Quantiles"
    )
  
  # Create detailed plots for a few random taxa
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  n_taxa <- min(n_taxa_to_plot, ncol(obs_data))
  selected_taxa <- sample(1:ncol(obs_data), n_taxa)
  
  # Create taxa fit plot data
  taxa_fit_data <- data.frame()
  
  for (t in selected_taxa) {
    # Observed counts for this taxon
    obs_counts <- obs_data[, t]
    
    # Expected counts for this taxon
    exp_counts <- sapply(1:nrow(obs_data), function(i) {
      generate_zinb_counts(1, mu[i, t], size, zi_prob[i, t])
    })
    
    # Add to plot data
    taxa_fit_data <- rbind(taxa_fit_data, data.frame(
      sample = 1:nrow(obs_data),
      observed = obs_counts,
      expected = exp_counts,
      taxon = paste0("Taxon_", t)
    ))
  }
  
  taxa_fit_plots <- ggplot2::ggplot(taxa_fit_data, ggplot2::aes(x = sample)) +
    ggplot2::geom_point(ggplot2::aes(y = observed), color = "blue", alpha = 0.7) +
    ggplot2::geom_point(ggplot2::aes(y = expected), color = "red", alpha = 0.7, shape = 4) +
    ggplot2::facet_wrap(~ taxon, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "Observed vs. Expected Counts for Selected Taxa",
      x = "Sample",
      y = "Count"
    )
  
  # Return all diagnostic plots
  return(list(
    zero_proportion_plot = zero_prop_plot,
    qq_plot = qq_plot,
    taxa_fit_plots = taxa_fit_plots
  ))
}