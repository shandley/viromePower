#' Plot Virus-Specific Zero Inflation Rates
#'
#' Visualizes the distribution of zero-inflation rates across viral taxa, which is useful
#' for understanding the heterogeneity of structural zeros in virome data. Different viral
#' taxa often have different prevalence patterns, leading to varying rates of structural zeros.
#'
#' @param sim_data A list output from simulate_zero_inflated_virome() function
#' @param binwidth Width of bins for histogram (default: 0.05)
#' @param include_density Whether to overlay density curve (default: TRUE)
#' @param color_palette Vector of colors for plot (default: c("navy", "firebrick"))
#'
#' @return A ggplot object showing the distribution of zero-inflation rates
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate simulated data with variable zero-inflation rates
#' sim_data <- simulate_zero_inflated_virome(
#'   n_samples = 30, 
#'   n_viruses = 100,
#'   variable_zi_rates = TRUE
#' )
#' 
#' # Plot distribution of zero-inflation rates
#' plot_virus_specific_zi_rates(sim_data)
#' }
plot_virus_specific_zi_rates <- function(sim_data, 
                                        binwidth = 0.05,
                                        include_density = TRUE,
                                        color_palette = c("navy", "firebrick")) {
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function")
  }
  
  # Verify input is from simulate_zero_inflated_virome
  if (!is.list(sim_data) || !"virus_specific_zi" %in% names(sim_data)) {
    stop("Input must be output from simulate_zero_inflated_virome() function")
  }
  
  # Extract zero-inflation rates
  zi_rates <- sim_data$virus_specific_zi$rates
  
  # Create data frame for plotting
  plot_data <- data.frame(
    virus_index = 1:length(zi_rates),
    zero_inflation_rate = zi_rates
  )
  
  # Create distribution plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = zero_inflation_rate)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      binwidth = binwidth,
      fill = color_palette[1],
      alpha = 0.7,
      color = "white"
    ) +
    ggplot2::labs(
      title = "Distribution of Zero-Inflation Rates Across Viral Taxa",
      subtitle = paste0(
        "Mean: ", round(sim_data$virus_specific_zi$mean, 2),
        ", Median: ", round(sim_data$virus_specific_zi$median, 2),
        ", Range: [", round(sim_data$virus_specific_zi$min, 2), 
        " - ", round(sim_data$virus_specific_zi$max, 2), "]"
      ),
      x = "Structural Zero Rate",
      y = "Density"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Add density line if requested
  if (include_density) {
    p <- p + 
      ggplot2::geom_density(
        color = color_palette[2],
        size = 1.2,
        alpha = 0.8
      )
  }
  
  # Add reference lines for mean and median
  p <- p +
    ggplot2::geom_vline(
      xintercept = sim_data$virus_specific_zi$mean,
      linetype = "dashed",
      color = color_palette[2],
      size = 0.8
    ) +
    ggplot2::geom_vline(
      xintercept = sim_data$virus_specific_zi$median,
      linetype = "dotted",
      color = color_palette[1],
      size = 0.8
    ) +
    ggplot2::annotate(
      "text",
      x = sim_data$virus_specific_zi$mean + 0.02,
      y = 0.8 * max(ggplot2::ggplot_build(p)$data[[1]]$density, na.rm = TRUE),
      label = "Mean",
      color = color_palette[2],
      hjust = 0
    ) +
    ggplot2::annotate(
      "text",
      x = sim_data$virus_specific_zi$median - 0.02,
      y = 0.7 * max(ggplot2::ggplot_build(p)$data[[1]]$density, na.rm = TRUE),
      label = "Median",
      color = color_palette[1],
      hjust = 1
    )
  
  return(p)
}

#' Plot Zero Inflation by Viral Abundance
#'
#' Visualizes the relationship between viral abundance and zero-inflation rates. This helps
#' understand whether rare viral taxa have higher rates of structural zeros, which is a common
#' pattern in virome data.
#'
#' @param sim_data A list output from simulate_zero_inflated_virome() function
#' @param abundance_scale Scale for abundance values ("log" or "linear", default: "log")
#' @param add_regression Whether to add regression line (default: TRUE)
#' @param color_palette Vector of colors for plot (default: c("navy", "firebrick"))
#'
#' @return A ggplot object showing relationship between abundance and zero-inflation rates
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate simulated data with variable zero-inflation rates
#' sim_data <- simulate_zero_inflated_virome(
#'   n_samples = 30, 
#'   n_viruses = 100,
#'   variable_zi_rates = TRUE
#' )
#' 
#' # Plot relationship between abundance and zero-inflation rates
#' plot_zi_by_abundance(sim_data)
#' }
plot_zi_by_abundance <- function(sim_data,
                                abundance_scale = "log",
                                add_regression = TRUE,
                                color_palette = c("navy", "firebrick")) {
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function")
  }
  
  # Verify input is from simulate_zero_inflated_virome
  if (!is.list(sim_data) || !"virus_specific_zi" %in% names(sim_data)) {
    stop("Input must be output from simulate_zero_inflated_virome() function")
  }
  
  # Extract zero-inflation rates
  zi_rates <- sim_data$virus_specific_zi$rates
  
  # Calculate mean abundance for each viral taxon
  abundance <- rowMeans(sim_data$counts)
  
  # Create data frame for plotting
  plot_data <- data.frame(
    virus_index = 1:length(zi_rates),
    zero_inflation_rate = zi_rates,
    abundance = abundance
  )
  
  # Create scatter plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = abundance, y = zero_inflation_rate)) +
    ggplot2::geom_point(
      color = color_palette[1],
      alpha = 0.7,
      size = 3
    ) +
    ggplot2::labs(
      title = "Relationship Between Viral Abundance and Zero-Inflation Rates",
      subtitle = "Higher zero-inflation rates typically occur in less abundant viral taxa",
      x = if (abundance_scale == "log") "Mean Abundance (log scale)" else "Mean Abundance",
      y = "Structural Zero Rate"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Apply log scale if requested
  if (abundance_scale == "log") {
    p <- p + ggplot2::scale_x_log10()
  }
  
  # Add regression line if requested
  if (add_regression) {
    if (abundance_scale == "log") {
      # For log scale, use log-transformed abundance for regression
      plot_data$log_abundance <- log10(plot_data$abundance + 1)
      
      p <- p + 
        ggplot2::geom_smooth(
          ggplot2::aes(x = abundance, y = zero_inflation_rate),
          method = "lm",
          formula = y ~ log(x + 1),
          color = color_palette[2],
          size = 1.2,
          alpha = 0.3
        )
    } else {
      p <- p + 
        ggplot2::geom_smooth(
          method = "lm",
          color = color_palette[2],
          size = 1.2,
          alpha = 0.3
        )
    }
  }
  
  return(p)
}

#' Compare Group-Specific Zero Inflation Patterns
#'
#' Visualizes and compares zero-inflation patterns between experimental groups, which is 
#' useful for understanding differential detection rates across groups.
#'
#' @param sim_data A list output from simulate_zero_inflated_virome() function
#' @param plot_type Type of comparison to plot: "boxplot" or "scatter" (default: "boxplot")
#' @param color_palette Vector of colors for plot (default: c("steelblue", "firebrick"))
#'
#' @return A ggplot object showing comparison of zero-inflation rates between groups
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate simulated data with variable zero-inflation rates and group differences
#' sim_data <- simulate_zero_inflated_virome(
#'   n_samples = 30, 
#'   n_viruses = 100,
#'   variable_zi_rates = TRUE,
#'   zero_inflation_difference = TRUE
#' )
#' 
#' # Compare zero-inflation patterns between groups
#' compare_group_zi_patterns(sim_data)
#' }
compare_group_zi_patterns <- function(sim_data,
                                     plot_type = "boxplot",
                                     color_palette = c("steelblue", "firebrick")) {
  # Check if required packages are available
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function")
  }
  
  # Verify input is from simulate_zero_inflated_virome
  if (!is.list(sim_data) || !"virus_specific_zi" %in% names(sim_data)) {
    stop("Input must be output from simulate_zero_inflated_virome() function")
  }
  
  # Extract structural zero probabilities for both groups
  if (!"structural_zero_probs" %in% names(sim_data$virus_specific_zi)) {
    stop("Group-specific structural zero probabilities not found in sim_data")
  }
  
  zi_probs <- sim_data$virus_specific_zi$structural_zero_probs
  
  # Create data frame for plotting
  plot_data <- data.frame(
    virus_index = rep(1:nrow(zi_probs), 2),
    group = rep(c("Group A", "Group B"), each = nrow(zi_probs)),
    structural_zero_rate = c(zi_probs[, 1], zi_probs[, 2])
  )
  
  # Check if zero inflation difference is present
  zi_diff_mean <- mean(zi_probs[, 2] - zi_probs[, 1])
  has_zi_diff <- abs(zi_diff_mean) > 0.05
  
  # Create plot based on requested type
  if (plot_type == "boxplot") {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = group, y = structural_zero_rate, fill = group)) +
      ggplot2::geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
      ggplot2::labs(
        title = "Comparison of Structural Zero Rates Between Groups",
        subtitle = if (has_zi_diff) 
          paste0("Mean difference: ", round(abs(zi_diff_mean) * 100, 1), "% ", 
                ifelse(zi_diff_mean < 0, "higher", "lower"), " in Group B") 
        else 
          "No substantial difference in structural zero rates between groups",
        x = "",
        y = "Structural Zero Rate"
      ) +
      ggplot2::scale_fill_manual(values = color_palette) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank(),
        legend.position = "none"
      )
    
    # Add individual points if not too many
    if (nrow(zi_probs) <= 100) {
      p <- p + ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 2)
    }
    
  } else if (plot_type == "scatter") {
    # Create scatter plot comparing Group A vs Group B rates
    p <- ggplot2::ggplot(data.frame(group_a = zi_probs[, 1], group_b = zi_probs[, 2]), 
                       ggplot2::aes(x = group_a, y = group_b)) +
      ggplot2::geom_point(alpha = 0.7, color = color_palette[1], size = 3) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
      ggplot2::labs(
        title = "Structural Zero Rates: Group A vs Group B",
        subtitle = if (has_zi_diff) 
          paste0("Group B has ", round(abs(zi_diff_mean) * 100, 1), "% ", 
                ifelse(zi_diff_mean < 0, "higher", "lower"), " structural zero rates on average") 
        else 
          "Similar structural zero rates between groups",
        x = "Group A Structural Zero Rate",
        y = "Group B Structural Zero Rate"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold"),
        panel.grid.minor = ggplot2::element_blank()
      ) +
      ggplot2::coord_equal() +
      ggplot2::xlim(0, 1) +
      ggplot2::ylim(0, 1)
    
    # Add loess smoothing line
    p <- p + ggplot2::geom_smooth(method = "loess", se = TRUE, color = color_palette[2], alpha = 0.2)
    
  } else {
    stop("plot_type must be 'boxplot' or 'scatter'")
  }
  
  return(p)
}