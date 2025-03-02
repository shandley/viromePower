#' Plot Power Analysis Results Across Different Zero-Inflation Patterns
#'
#' Creates visualizations comparing statistical power across different zero-inflation
#' patterns, helping researchers understand how various zero-inflation distributions
#' affect study power for virome analysis.
#'
#' @param sample_sizes Vector of sample sizes to evaluate (default: 10, 20, 30, 40, 50)
#' @param effect_sizes Vector of effect sizes to evaluate (default: 1.5, 2.0, 2.5)
#' @param n_viruses Number of viral taxa to simulate (default: 100)
#' @param zi_patterns Named list of alpha/beta parameter pairs for different zero-inflation patterns
#'   (default provides right-skewed, left-skewed, uniform, and symmetric patterns)
#' @param n_sim Number of simulations per parameter combination (default: 10)
#' @param fixed_zi_rate Mean zero-inflation rate for all patterns (default: 0.7)
#' @param sampling_zeros Proportion of sampling zeros (default: 0.2)
#' @param plot_type Type of plot to generate: "heatmap", "line", or "combined" (default: "combined")
#' @param parallel Whether to use parallel processing (default: FALSE)
#' @param n_cores Number of CPU cores for parallel processing (default: detectCores() - 1)
#'
#' @return A list containing:
#'   \item{plots}{List of generated plots}
#'   \item{power_results}{Data frame with power results for all parameter combinations}
#'   \item{summary}{Summary statistics comparing patterns}
#'
#' @details
#' This function systematically evaluates statistical power across different zero-inflation
#' patterns to help researchers understand how the distribution of structural zeros affects
#' the ability to detect differences between groups in virome data.
#'
#' The function:
#' 1. Simulates data for each parameter combination using different zero-inflation patterns
#' 2. Calculates power for each scenario using ZINB models
#' 3. Creates visualizations showing how power varies with sample size, effect size, and ZI pattern
#'
#' The default zero-inflation patterns provided are:
#' \itemize{
#'   \item "right_skewed": Beta(2, 5) - More taxa with lower ZI rates
#'   \item "left_skewed": Beta(5, 2) - More taxa with higher ZI rates
#'   \item "uniform": Beta(1, 1) - Uniform distribution of ZI rates
#'   \item "symmetric": Beta(2, 2) - Bell-shaped distribution centered at 0.5
#'   \item "fixed": No variability, all taxa have the same ZI rate
#' }
#'
#' These patterns maintain the same mean ZI rate for fair comparison, but differ in their
#' distribution shapes, representing different prevalence patterns in virome data.
#'
#' @examples
#' \dontrun{
#' # Basic comparison of different ZI patterns
#' power_comparison <- plot_power_by_zi_pattern(
#'   sample_sizes = c(20, 40),
#'   effect_sizes = c(1.5, 2.5),
#'   n_sim = 5  # Use higher for actual analysis (e.g., 20-50)
#' )
#' 
#' # Print the first plot
#' print(power_comparison$plots$combined)
#' 
#' # Custom ZI patterns
#' custom_patterns <- list(
#'   highly_variable = c(alpha = 0.5, beta = 0.5),
#'   mostly_present = c(alpha = 1, beta = 5),
#'   mostly_absent = c(alpha = 5, beta = 1)
#' )
#' 
#' power_custom <- plot_power_by_zi_pattern(
#'   sample_sizes = c(15, 30, 45),
#'   effect_sizes = 2.0,
#'   zi_patterns = custom_patterns,
#'   plot_type = "line"
#' )
#' }
#'
#' @export
plot_power_by_zi_pattern <- function(sample_sizes = c(10, 20, 30, 40, 50),
                                    effect_sizes = c(1.5, 2.0, 2.5),
                                    n_viruses = 100,
                                    zi_patterns = NULL,
                                    n_sim = 10,
                                    fixed_zi_rate = 0.7,
                                    sampling_zeros = 0.2,
                                    plot_type = "combined",
                                    parallel = FALSE,
                                    n_cores = NULL) {
  
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for plotting")
  }
  
  # Initialize parallel processing if requested
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      warning("Package 'parallel' is required for parallel processing. Falling back to sequential.")
      parallel <- FALSE
    } else {
      if (is.null(n_cores)) {
        n_cores <- max(1, parallel::detectCores() - 1)
      }
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl))
      
      # Export required functions to the cluster
      parallel::clusterExport(cl, c("simulate_zero_inflated_virome", "calc_zinb_bayesian_power"),
                             envir = environment())
    }
  }
  
  # Define default zero-inflation patterns if not provided
  if (is.null(zi_patterns)) {
    zi_patterns <- list(
      right_skewed = c(alpha = 2, beta = 5),    # More taxa with lower ZI rates
      left_skewed = c(alpha = 5, beta = 2),     # More taxa with higher ZI rates
      uniform = c(alpha = 1, beta = 1),         # Uniform distribution
      symmetric = c(alpha = 2, beta = 2),       # Bell-shaped, centered at 0.5
      fixed = c(alpha = NA, beta = NA)          # Fixed rate (traditional approach)
    )
  }
  
  # Create all parameter combinations
  param_grid <- expand.grid(
    n_samples = sample_sizes,
    effect_size = effect_sizes,
    pattern_name = names(zi_patterns),
    stringsAsFactors = FALSE
  )
  
  # Initialize results data frame
  results <- data.frame(
    n_samples = param_grid$n_samples,
    effect_size = param_grid$effect_size,
    pattern_name = param_grid$pattern_name,
    power = NA,
    ci_lower = NA,
    ci_upper = NA,
    zi_alpha = NA,
    zi_beta = NA,
    stringsAsFactors = FALSE
  )
  
  # Add alpha/beta parameters
  for (i in 1:nrow(results)) {
    pattern <- results$pattern_name[i]
    if (pattern != "fixed") {
      results$zi_alpha[i] <- zi_patterns[[pattern]]["alpha"]
      results$zi_beta[i] <- zi_patterns[[pattern]]["beta"]
    }
  }
  
  # Calculate power for each combination
  calculate_power <- function(i) {
    row <- results[i, ]
    pattern_name <- row$pattern_name
    
    # Fixed or variable ZI rates
    if (pattern_name == "fixed") {
      variable_zi <- FALSE
      zi_alpha <- NA
      zi_beta <- NA
    } else {
      variable_zi <- TRUE
      zi_alpha <- zi_patterns[[pattern_name]]["alpha"]
      zi_beta <- zi_patterns[[pattern_name]]["beta"]
    }
    
    # Run power analysis
    power_result <- calc_zinb_bayesian_power(
      n_samples = row$n_samples,
      effect_size = row$effect_size,
      n_viruses = n_viruses,
      structural_zeros = fixed_zi_rate,
      sampling_zeros = sampling_zeros,
      n_sim = n_sim,
      variable_zi_rates = variable_zi,
      zi_alpha = zi_alpha,
      zi_beta = zi_beta
    )
    
    # Return power result
    c(
      power = power_result$power,
      ci_lower = power_result$ci_lower,
      ci_upper = power_result$ci_upper
    )
  }
  
  # Apply the calculation function to each row
  if (parallel) {
    power_results <- parallel::parApply(cl, matrix(1:nrow(results)), 1, calculate_power)
    power_results <- t(power_results)  # Transpose to get the correct structure
  } else {
    power_results <- t(sapply(1:nrow(results), calculate_power))
  }
  
  # Update results with power values
  results$power <- power_results[, "power"]
  results$ci_lower <- power_results[, "ci_lower"]
  results$ci_upper <- power_results[, "ci_upper"]
  
  # Add pattern descriptions
  pattern_descriptions <- list(
    right_skewed = "Right-skewed (more taxa with lower ZI)",
    left_skewed = "Left-skewed (more taxa with higher ZI)",
    uniform = "Uniform distribution",
    symmetric = "Symmetric (bell-shaped)",
    fixed = "Fixed (all taxa same ZI rate)"
  )
  
  results$pattern_label <- sapply(results$pattern_name, function(p) {
    if (p %in% names(pattern_descriptions)) {
      return(pattern_descriptions[[p]])
    } else {
      return(p)
    }
  })
  
  # Create factor with ordered levels for plotting
  results$pattern_label <- factor(results$pattern_label, 
                                levels = unlist(pattern_descriptions[names(zi_patterns)]))
  
  # Prepare plots
  plots <- list()
  
  # Heatmap: Power by sample size and effect size for each pattern
  if (plot_type %in% c("heatmap", "combined")) {
    # Create a heatmap for each pattern
    for (pattern in unique(results$pattern_name)) {
      pattern_data <- results[results$pattern_name == pattern, ]
      
      p <- ggplot2::ggplot(pattern_data, 
                          ggplot2::aes(x = n_samples, y = effect_size, fill = power)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
          low = "white", mid = "steelblue", high = "darkblue",
          midpoint = 0.5, limits = c(0, 1), name = "Power"
        ) +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", power)), 
                          color = "white", size = 3) +
        ggplot2::labs(
          title = pattern_descriptions[[pattern]],
          x = "Sample Size per Group",
          y = "Effect Size"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          legend.position = "right"
        )
      
      plots[[paste0("heatmap_", pattern)]] <- p
    }
    
    # Combined heatmap with facets
    p_facet <- ggplot2::ggplot(results, 
                              ggplot2::aes(x = n_samples, y = effect_size, fill = power)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "white", mid = "steelblue", high = "darkblue",
        midpoint = 0.5, limits = c(0, 1), name = "Power"
      ) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", power)), 
                        color = "white", size = 3) +
      ggplot2::facet_wrap(~pattern_label) +
      ggplot2::labs(
        title = "Power by Sample Size and Effect Size Across Zero-Inflation Patterns",
        x = "Sample Size per Group",
        y = "Effect Size"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "right",
        strip.text = ggplot2::element_text(face = "bold")
      )
    
    plots[["heatmap_facet"]] <- p_facet
  }
  
  # Line plot: Power by sample size for each pattern and effect size
  if (plot_type %in% c("line", "combined")) {
    # For each effect size, create a line plot
    for (es in unique(results$effect_size)) {
      es_data <- results[results$effect_size == es, ]
      
      p <- ggplot2::ggplot(es_data, 
                          ggplot2::aes(x = n_samples, y = power, 
                                       color = pattern_label, group = pattern_label)) +
        ggplot2::geom_line(size = 1) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
          width = 0.2, alpha = 0.7
        ) +
        ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50") +
        ggplot2::labs(
          title = sprintf("Power by Sample Size for Effect Size %.1f", es),
          subtitle = "Comparing different zero-inflation patterns",
          x = "Sample Size per Group",
          y = "Statistical Power",
          color = "ZI Pattern"
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = ggplot2::element_text(face = "bold")
        )
      
      plots[[paste0("line_es", es)]] <- p
    }
    
    # For each sample size, create a line plot
    for (ss in unique(results$n_samples)) {
      ss_data <- results[results$n_samples == ss, ]
      
      p <- ggplot2::ggplot(ss_data, 
                          ggplot2::aes(x = effect_size, y = power, 
                                       color = pattern_label, group = pattern_label)) +
        ggplot2::geom_line(size = 1) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_errorbar(
          ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
          width = 0.05, alpha = 0.7
        ) +
        ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50") +
        ggplot2::labs(
          title = sprintf("Power by Effect Size for N=%d per Group", ss),
          subtitle = "Comparing different zero-inflation patterns",
          x = "Effect Size",
          y = "Statistical Power",
          color = "ZI Pattern"
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 1)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          legend.position = "bottom",
          legend.title = ggplot2::element_text(face = "bold")
        )
      
      plots[[paste0("line_ss", ss)]] <- p
    }
  }
  
  # Create a combined plot showing the most important comparisons
  if (plot_type == "combined") {
    # Choose a representative effect size and sample size for focused plots
    mid_effect <- median(effect_sizes)
    mid_sample <- median(sample_sizes)
    
    # Subset data for representative plots
    es_data <- results[results$effect_size == mid_effect, ]
    ss_data <- results[results$n_samples == mid_sample, ]
    
    # Create panels
    p1 <- ggplot2::ggplot(es_data, 
                         ggplot2::aes(x = n_samples, y = power, 
                                      color = pattern_label, group = pattern_label)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        width = 0.2, alpha = 0.7
      ) +
      ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50") +
      ggplot2::labs(
        title = sprintf("Effect Size = %.1f", mid_effect),
        x = "Sample Size per Group",
        y = "Statistical Power"
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "none"
      )
    
    p2 <- ggplot2::ggplot(ss_data, 
                         ggplot2::aes(x = effect_size, y = power, 
                                      color = pattern_label, group = pattern_label)) +
      ggplot2::geom_line(size = 1) +
      ggplot2::geom_point(size = 2) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = ci_lower, ymax = ci_upper),
        width = 0.05, alpha = 0.7
      ) +
      ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray50") +
      ggplot2::labs(
        title = sprintf("Sample Size = %d", mid_sample),
        x = "Effect Size",
        y = ""
      ) +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        legend.position = "none"
      )
    
    # Create a small version of the heatmap facet plot
    p3 <- plots[["heatmap_facet"]] +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 10),
        axis.title = ggplot2::element_text(size = 8),
        axis.text = ggplot2::element_text(size = 6),
        strip.text = ggplot2::element_text(size = 8),
        legend.title = ggplot2::element_text(size = 8),
        legend.text = ggplot2::element_text(size = 6)
      )
    
    # Add external legend for the combined plot
    legend_data <- ggplot2::ggplot(results, 
                                  ggplot2::aes(x = n_samples, y = power, 
                                               color = pattern_label)) +
      ggplot2::geom_line() +
      ggplot2::theme(legend.position = "bottom") +
      ggplot2::guides(color = ggplot2::guide_legend(title = "Zero-Inflation Pattern", 
                                                   nrow = 2))
    
    legend <- ggplot2::get_legend(legend_data)
    
    # Combine plots
    if (requireNamespace("gridExtra", quietly = TRUE)) {
      combined_plot <- gridExtra::grid.arrange(
        gridExtra::arrangeGrob(p1, p2, ncol = 2),
        p3,
        legend,
        heights = c(0.4, 0.4, 0.2),
        top = "Comparison of Power Across Zero-Inflation Patterns"
      )
      plots[["combined"]] <- combined_plot
    } else {
      warning("Package 'gridExtra' is required for combined plots. Individual plots will still be generated.")
    }
  }
  
  # Create summary statistics
  summary <- list(
    avg_power_by_pattern = aggregate(power ~ pattern_name, data = results, FUN = mean),
    avg_power_by_pattern_and_es = aggregate(
      power ~ pattern_name + effect_size, 
      data = results, 
      FUN = mean
    ),
    avg_power_by_pattern_and_ss = aggregate(
      power ~ pattern_name + n_samples, 
      data = results, 
      FUN = mean
    ),
    min_sample_size_80pct = NULL
  )
  
  # Calculate minimum sample size for 80% power by pattern and effect size
  min_n_for_power <- function(pattern, es, target_power = 0.8) {
    pattern_es_data <- results[results$pattern_name == pattern & results$effect_size == es, ]
    if (all(pattern_es_data$power < target_power)) {
      return(NA)  # Can't achieve target power
    }
    min_idx <- min(which(pattern_es_data$power >= target_power))
    return(pattern_es_data$n_samples[min_idx])
  }
  
  min_n_df <- expand.grid(
    pattern_name = unique(results$pattern_name),
    effect_size = unique(results$effect_size),
    min_n_for_80pct = NA
  )
  
  for (i in 1:nrow(min_n_df)) {
    min_n_df$min_n_for_80pct[i] <- min_n_for_power(
      min_n_df$pattern_name[i], 
      min_n_df$effect_size[i]
    )
  }
  
  summary$min_sample_size_80pct <- min_n_df
  
  # Return results
  return(list(
    plots = plots,
    power_results = results,
    summary = summary
  ))
}