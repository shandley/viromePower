#' Plot Viral Diversity Power Curve
#'
#' Creates a power curve visualization for viral diversity analysis, showing how power
#' changes with sample size or effect size.
#'
#' @param param_range Vector of parameter values to test (sample sizes or effect sizes)
#' @param param_type Type of parameter being varied: "n_samples" or "effect_size" (default: "n_samples")
#' @param effect_size Effect size to use when varying sample size (default: 1.0)
#' @param n_samples Number of samples to use when varying effect size (default: 15)
#' @param n_viruses Number of viral taxa to simulate (default: 200)
#' @param diversity_measure Type of diversity to analyze: "shannon", "simpson", "richness", "evenness", "bray", "jaccard", or "unifrac" (default: "shannon")
#' @param alpha Significance level (default: 0.05)
#' @param sparsity Proportion of zeros in the data (default: 0.8)
#' @param dispersion Dispersion parameter for viral abundance (default: 2)
#' @param n_sim Number of simulations for each point (default: 50)
#' @param show_individual_sims Logical indicating whether to show individual simulation points (default: FALSE)
#'
#' @return A ggplot object with the power curve visualization
#' @export
#'
#' @examples
#' # Power curve for varying sample sizes with Shannon diversity
#' power_curve <- plot_diversity_power_curve(param_range = seq(5, 30, by = 5), 
#'                                          diversity_measure = "shannon")
#'                                          
#' # Power curve for varying effect sizes with Bray-Curtis dissimilarity
#' effect_curve <- plot_diversity_power_curve(param_range = seq(0.05, 0.3, by = 0.05),
#'                                           param_type = "effect_size",
#'                                           n_samples = 20,
#'                                           diversity_measure = "bray")
plot_diversity_power_curve <- function(param_range, param_type = "n_samples", 
                                     effect_size = 1.0, n_samples = 15, 
                                     n_viruses = 200, diversity_measure = "shannon",
                                     alpha = 0.05, sparsity = 0.8, dispersion = 2,
                                     n_sim = 50, show_individual_sims = FALSE) {
  
  # Validate input parameters
  if (length(param_range) < 2) {
    stop("param_range must contain at least 2 values")
  }
  if (!(param_type %in% c("n_samples", "effect_size"))) {
    stop("param_type must be either 'n_samples' or 'effect_size'")
  }
  
  # Initialize results storage
  results <- data.frame(
    param_value = numeric(0),
    power = numeric(0),
    sim_id = numeric(0)
  )
  
  # For each parameter value, calculate power
  for (param_value in param_range) {
    if (param_type == "n_samples") {
      # Vary sample size, keep effect size constant
      power_result <- calc_viral_diversity_power(
        n_samples = param_value,
        effect_size = effect_size,
        n_viruses = n_viruses,
        alpha = alpha,
        sparsity = sparsity,
        dispersion = dispersion,
        diversity_measure = diversity_measure,
        n_sim = n_sim
      )
    } else {
      # Vary effect size, keep sample size constant
      power_result <- calc_viral_diversity_power(
        n_samples = n_samples,
        effect_size = param_value,
        n_viruses = n_viruses,
        alpha = alpha,
        sparsity = sparsity,
        dispersion = dispersion,
        diversity_measure = diversity_measure,
        n_sim = n_sim
      )
    }
    
    # Extract power and individual simulation results
    if (show_individual_sims) {
      # Add individual simulation results
      sim_results <- data.frame(
        param_value = rep(param_value, n_sim),
        power = as.numeric(power_result$significant_tests),
        sim_id = 1:n_sim
      )
      results <- rbind(results, sim_results)
    } else {
      # Just add the average power
      results <- rbind(results, data.frame(
        param_value = param_value,
        power = power_result$power,
        sim_id = 0
      ))
    }
  }
  
  # Create x-axis label based on parameter type
  if (param_type == "n_samples") {
    x_label <- "Number of Samples per Group"
  } else {
    x_label <- "Effect Size"
  }
  
  # Format diversity measure name for title
  diversity_title <- switch(diversity_measure,
                            shannon = "Shannon Diversity",
                            simpson = "Simpson Diversity",
                            richness = "Species Richness",
                            evenness = "Pielou's Evenness",
                            bray = "Bray-Curtis Dissimilarity",
                            jaccard = "Jaccard Distance",
                            unifrac = "UniFrac Distance",
                            "Diversity")
  
  # Create plot title
  if (param_type == "n_samples") {
    plot_title <- paste0("Power Analysis for ", diversity_title, 
                         " (Effect Size = ", effect_size, ")")
  } else {
    plot_title <- paste0("Power Analysis for ", diversity_title, 
                         " (n = ", n_samples, " per group)")
  }
  
  # Create power curve plot using ggplot2
  power_plot <- ggplot2::ggplot(results, ggplot2::aes(x = param_value, y = power)) +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "darkgray") +
    ggplot2::annotate("text", x = max(param_range) * 0.9, y = 0.82, 
                     label = "80% Power", color = "darkgray") +
    ggplot2::labs(
      title = plot_title,
      x = x_label,
      y = "Statistical Power",
      caption = paste0("Based on ", n_sim, " simulations per point, ", 
                      n_viruses, " viral taxa, alpha = ", alpha)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2))
  
  if (show_individual_sims) {
    # Add individual simulation points (jittered for visibility)
    power_plot <- power_plot +
      ggplot2::geom_point(alpha = 0.2, position = ggplot2::position_jitter(width = 0.1, height = 0)) +
      ggplot2::stat_summary(fun = mean, geom = "line", size = 1, color = "blue") +
      ggplot2::stat_summary(fun = mean, geom = "point", size = 3, color = "blue")
  } else {
    # Just show the average power curve
    power_plot <- power_plot +
      ggplot2::geom_line(size = 1, color = "blue") +
      ggplot2::geom_point(size = 3, color = "blue")
  }
  
  return(power_plot)
}