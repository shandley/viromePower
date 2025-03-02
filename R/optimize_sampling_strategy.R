#' Optimize Sampling Strategy for Virome Studies with Variable Zero-Inflation
#'
#' Determines the optimal combination of sample number, sequencing depth, and
#' statistical approach to maximize power given specific zero-inflation patterns
#' and study constraints.
#'
#' @param target_power Desired statistical power (default: 0.8)
#' @param effect_size Expected effect size to detect (default: 2.0)
#' @param zi_alpha Shape parameter alpha for Beta distribution of zero-inflation rates (default: 2)
#' @param zi_beta Shape parameter beta for Beta distribution of zero-inflation rates (default: 5)
#' @param total_sequencing_effort Total sequencing capacity in millions of reads (default: 100)
#' @param mean_zi Expected mean zero-inflation rate (default: 0.7)
#' @param sampling_zeros Expected sampling zero rate (default: 0.2)
#' @param budget_constraint Whether to apply a budget constraint (default: TRUE)
#' @param depth_range Range of sequencing depths to evaluate in reads per sample (default: 0.5-5 million)
#' @param sample_range Range of sample sizes to evaluate (default: 10-100)
#' @param test_variable_zi Whether to test variable zero-inflation modeling (default: TRUE)
#' @param visualize Whether to create visualization of results (default: TRUE)
#'
#' @return A list containing:
#'   \item{optimal_strategy}{The recommended sampling strategy}
#'   \item{power_results}{Data frame with power estimates for different strategies}
#'   \item{contour_plot}{ggplot object showing power by sample size and depth}
#'   \item{comparison_plot}{ggplot object comparing fixed vs variable ZI approaches}
#'   \item{recommendation}{Text string summarizing the recommendation}
#'
#' @details
#' This function helps researchers design optimal sampling strategies for virome studies
#' by systematically evaluating the trade-off between sample size and sequencing depth
#' under different zero-inflation patterns.
#'
#' The function:
#' 1. Simulates power across a grid of sample sizes and sequencing depths
#' 2. Considers the total sequencing effort constraint (sample size × depth)
#' 3. Evaluates both fixed and variable zero-inflation modeling approaches
#' 4. Recommends the most efficient sampling strategy to achieve the target power
#'
#' For virome studies, there is often a trade-off between:
#' - Fewer samples with higher sequencing depth (better for rare virus detection)
#' - More samples with lower sequencing depth (better for prevalence patterns)
#'
#' This function helps navigate this trade-off based on the expected zero-inflation 
#' pattern in the dataset. The Beta distribution parameters (zi_alpha, zi_beta) control
#' the expected distribution of structural zero rates across viral taxa.
#'
#' @examples
#' \dontrun{
#' # Basic optimization with default parameters
#' opt_result <- optimize_sampling_strategy(
#'   target_power = 0.8,
#'   effect_size = 2.0
#' )
#' 
#' # Print the recommendation
#' cat(opt_result$recommendation)
#' 
#' # Plot the results
#' print(opt_result$contour_plot)
#' 
#' # Custom optimization for a specific scenario
#' custom_opt <- optimize_sampling_strategy(
#'   target_power = 0.9,
#'   effect_size = 1.8,
#'   zi_alpha = 1,
#'   zi_beta = 3,
#'   total_sequencing_effort = 200,
#'   sample_range = c(20, 150)
#' )
#' }
#'
#' @export
optimize_sampling_strategy <- function(target_power = 0.8,
                                      effect_size = 2.0,
                                      zi_alpha = 2,
                                      zi_beta = 5,
                                      total_sequencing_effort = 100,
                                      mean_zi = 0.7,
                                      sampling_zeros = 0.2,
                                      budget_constraint = TRUE,
                                      depth_range = c(0.5, 5),
                                      sample_range = c(10, 100),
                                      test_variable_zi = TRUE,
                                      visualize = TRUE) {
  
  # Check if required packages are available
  if (visualize && !requireNamespace("ggplot2", quietly = TRUE)) {
    warning("Package 'ggplot2' is required for visualization. No plots will be generated.")
    visualize <- FALSE
  }
  
  # Define grid of parameters to evaluate
  sample_sizes <- floor(seq(sample_range[1], sample_range[2], length.out = 10))
  
  if (budget_constraint) {
    # Calculate depths based on total sequencing effort and sample sizes
    # Constraint: n_samples * depth = total_sequencing_effort
    depths <- total_sequencing_effort / sample_sizes
    
    # Filter out depths outside the specified range
    valid_idx <- depths >= depth_range[1] & depths <= depth_range[2]
    
    if (sum(valid_idx) < 3) {
      warning("Very few valid combinations within constraints. Consider adjusting depth_range or total_sequencing_effort.")
      
      # Adjust ranges to ensure at least some valid combinations
      sample_sizes <- floor(seq(max(5, sample_range[1]), 
                             min(sample_range[2], total_sequencing_effort / depth_range[1]), 
                             length.out = 8))
      depths <- total_sequencing_effort / sample_sizes
      valid_idx <- depths >= depth_range[1] & depths <= depth_range[2]
    }
    
    sample_sizes <- sample_sizes[valid_idx]
    depths <- depths[valid_idx]
    
    # Create parameter grid
    param_grid <- data.frame(
      n_samples = sample_sizes,
      depth = depths,
      total_effort = sample_sizes * depths,
      stringsAsFactors = FALSE
    )
    
  } else {
    # Without budget constraint, evaluate a grid of sample sizes and depths
    depths <- seq(depth_range[1], depth_range[2], length.out = 8)
    
    # Create all combinations
    param_grid <- expand.grid(
      n_samples = sample_sizes,
      depth = depths,
      stringsAsFactors = FALSE
    )
    
    param_grid$total_effort <- param_grid$n_samples * param_grid$depth
  }
  
  # Add variable/fixed ZI columns if testing both approaches
  if (test_variable_zi) {
    param_grid <- rbind(
      cbind(param_grid, zi_modeling = "fixed"),
      cbind(param_grid, zi_modeling = "variable")
    )
  } else {
    param_grid$zi_modeling <- "fixed"
  }
  
  # Initialize results
  param_grid$power <- NA
  param_grid$reads_per_virus <- NA
  param_grid$cost_efficiency <- NA
  
  # Define number of viruses based on depth
  base_n_viruses <- 1000  # Base number at 1M depth
  
  # Function to estimate power for a given parameter set
  estimate_power <- function(row) {
    # Estimate number of viruses based on sequencing depth
    # This is a simplified model: more depth → more viruses detected
    # In reality, this relationship is non-linear and plateaus
    row$reads_per_virus <- row$depth * 1e6 / base_n_viruses
    depth_factor <- min(3, 0.5 + log10(row$depth) / log10(5))
    n_viruses <- ceiling(base_n_viruses * depth_factor)
    
    # Set variable ZI parameters based on modeling approach
    if (row$zi_modeling == "variable") {
      variable_zi <- TRUE
      local_zi_alpha <- zi_alpha
      local_zi_beta <- zi_beta
    } else {
      variable_zi <- FALSE
      local_zi_alpha <- NA
      local_zi_beta <- NA
    }
    
    # Calculate power
    power_result <- calc_zinb_bayesian_power(
      n_samples = row$n_samples,
      effect_size = effect_size,
      n_viruses = n_viruses,
      structural_zeros = mean_zi,
      sampling_zeros = sampling_zeros,
      variable_zi_rates = variable_zi,
      zi_alpha = local_zi_alpha,
      zi_beta = local_zi_beta,
      n_sim = 10  # Low number for speed; increase for more precision
    )
    
    # Calculate cost efficiency (power per unit of sequencing effort)
    cost_efficiency <- power_result$power / row$total_effort
    
    c(
      power = power_result$power,
      reads_per_virus = row$reads_per_virus,
      cost_efficiency = cost_efficiency
    )
  }
  
  # Apply the function to each row
  power_results <- t(apply(param_grid, 1, estimate_power))
  
  # Update the results dataframe
  param_grid$power <- power_results[, "power"]
  param_grid$reads_per_virus <- power_results[, "reads_per_virus"]
  param_grid$cost_efficiency <- power_results[, "cost_efficiency"]
  
  # Find optimal strategies
  
  # 1. Strategy that achieves target power with minimum sequencing effort
  target_power_rows <- param_grid$power >= target_power
  
  if (sum(target_power_rows) > 0) {
    min_effort_strategy <- param_grid[target_power_rows, ][
      which.min(param_grid$total_effort[target_power_rows]), 
    ]
  } else {
    # If no strategy achieves target power, find the one with highest power
    min_effort_strategy <- param_grid[which.max(param_grid$power), ]
  }
  
  # 2. Most cost-efficient strategy (highest power per sequencing effort)
  cost_efficient_strategy <- param_grid[which.max(param_grid$cost_efficiency), ]
  
  # 3. Strategy with highest absolute power
  max_power_strategy <- param_grid[which.max(param_grid$power), ]
  
  # Prepare the recommendation
  if (min_effort_strategy$power >= target_power) {
    primary_recommendation <- min_effort_strategy
    recommendation_type <- "minimum_effort"
  } else {
    # If target power not achievable, recommend most powerful strategy
    primary_recommendation <- max_power_strategy
    recommendation_type <- "maximum_power"
  }
  
  # Create visualizations if requested
  if (visualize) {
    plots <- list()
    
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # For contour plot, we need to pivot and prepare the data
      if (test_variable_zi) {
        # Separate fixed and variable approaches for comparison
        fixed_data <- param_grid[param_grid$zi_modeling == "fixed", ]
        var_data <- param_grid[param_grid$zi_modeling == "variable", ]
        
        # Create comparison plot
        comparison_df <- data.frame(
          n_samples = rep(unique(param_grid$n_samples), 2),
          fixed_power = sapply(unique(param_grid$n_samples), function(n) {
            max(fixed_data$power[fixed_data$n_samples == n])
          }),
          variable_power = sapply(unique(param_grid$n_samples), function(n) {
            max(var_data$power[var_data$n_samples == n])
          }),
          stringsAsFactors = FALSE
        )
        
        # Long format for plotting
        comparison_long <- rbind(
          data.frame(
            n_samples = comparison_df$n_samples,
            power = comparison_df$fixed_power,
            approach = "Fixed ZI",
            stringsAsFactors = FALSE
          ),
          data.frame(
            n_samples = comparison_df$n_samples,
            power = comparison_df$variable_power,
            approach = "Variable ZI",
            stringsAsFactors = FALSE
          )
        )
        
        p_comparison <- ggplot2::ggplot(comparison_long, 
                                       ggplot2::aes(x = n_samples, y = power, 
                                                   color = approach, group = approach)) +
          ggplot2::geom_line(size = 1.2) +
          ggplot2::geom_point(size = 3) +
          ggplot2::geom_hline(yintercept = target_power, linetype = "dashed", color = "gray50") +
          ggplot2::labs(
            title = "Power Comparison: Fixed vs. Variable Zero-Inflation Modeling",
            subtitle = sprintf("Effect size = %.1f, Target power = %.2f", effect_size, target_power),
            x = "Number of Samples",
            y = "Statistical Power",
            color = "Modeling Approach"
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(
            plot.title = ggplot2::element_text(face = "bold"),
            legend.position = "bottom"
          )
        
        plots$comparison_plot <- p_comparison
        
        # Create separate contour plots for each approach
        for (approach in c("fixed", "variable")) {
          approach_data <- param_grid[param_grid$zi_modeling == approach, ]
          
          if (budget_constraint) {
            # With budget constraint, we have paired samples/depths
            # Create a scatter plot instead of contour
            p <- ggplot2::ggplot(approach_data, 
                               ggplot2::aes(x = n_samples, y = depth, color = power, size = power)) +
              ggplot2::geom_point(alpha = 0.8) +
              ggplot2::scale_color_gradient2(
                low = "white", mid = "steelblue", high = "darkblue",
                midpoint = 0.5, limits = c(0, 1), name = "Power"
              ) +
              ggplot2::scale_size_continuous(range = c(3, 8), guide = "none") +
              ggplot2::labs(
                title = sprintf("Power by Sample Size and Depth (%s ZI)", 
                               ifelse(approach == "fixed", "Fixed", "Variable")),
                subtitle = sprintf("Total sequencing effort: %d million reads", total_sequencing_effort),
                x = "Number of Samples",
                y = "Sequencing Depth (million reads/sample)"
              ) +
              ggplot2::theme_minimal() +
              ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold")
              )
          } else {
            # Without budget constraint, create a proper contour plot
            p <- ggplot2::ggplot(approach_data, 
                               ggplot2::aes(x = n_samples, y = depth, z = power)) +
              ggplot2::geom_contour_filled(binwidth = 0.1) +
              ggplot2::scale_fill_viridis_d(name = "Power") +
              ggplot2::labs(
                title = sprintf("Power Contours (%s ZI)", 
                               ifelse(approach == "fixed", "Fixed", "Variable")),
                x = "Number of Samples",
                y = "Sequencing Depth (million reads/sample)"
              ) +
              ggplot2::theme_minimal() +
              ggplot2::theme(
                plot.title = ggplot2::element_text(face = "bold"),
                legend.position = "right"
              )
          }
          
          # Add marker for optimal strategy
          opt_row <- primary_recommendation
          if (opt_row$zi_modeling == approach) {
            p <- p + 
              ggplot2::geom_point(
                data = opt_row,
                ggplot2::aes(x = n_samples, y = depth),
                color = "red", size = 5, shape = 8
              ) +
              ggplot2::annotate(
                "text",
                x = opt_row$n_samples,
                y = opt_row$depth * 1.1,
                label = "Optimal",
                color = "red",
                fontface = "bold"
              )
          }
          
          plots[[paste0("contour_", approach)]] <- p
        }
        
      } else {
        # Single contour plot for fixed ZI approach
        if (budget_constraint) {
          p <- ggplot2::ggplot(param_grid, 
                             ggplot2::aes(x = n_samples, y = depth, color = power, size = power)) +
            ggplot2::geom_point(alpha = 0.8) +
            ggplot2::scale_color_gradient2(
              low = "white", mid = "steelblue", high = "darkblue",
              midpoint = 0.5, limits = c(0, 1), name = "Power"
            ) +
            ggplot2::scale_size_continuous(range = c(3, 8), guide = "none") +
            ggplot2::labs(
              title = "Power by Sample Size and Depth",
              subtitle = sprintf("Total sequencing effort: %d million reads", total_sequencing_effort),
              x = "Number of Samples",
              y = "Sequencing Depth (million reads/sample)"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold")
            )
        } else {
          p <- ggplot2::ggplot(param_grid, 
                             ggplot2::aes(x = n_samples, y = depth, z = power)) +
            ggplot2::geom_contour_filled(binwidth = 0.1) +
            ggplot2::scale_fill_viridis_d(name = "Power") +
            ggplot2::labs(
              title = "Power Contours",
              x = "Number of Samples",
              y = "Sequencing Depth (million reads/sample)"
            ) +
            ggplot2::theme_minimal() +
            ggplot2::theme(
              plot.title = ggplot2::element_text(face = "bold"),
              legend.position = "right"
            )
        }
        
        # Add marker for optimal strategy
        opt_row <- primary_recommendation
        p <- p + 
          ggplot2::geom_point(
            data = opt_row,
            ggplot2::aes(x = n_samples, y = depth),
            color = "red", size = 5, shape = 8
          ) +
          ggplot2::annotate(
            "text",
            x = opt_row$n_samples,
            y = opt_row$depth * 1.1,
            label = "Optimal",
            color = "red",
            fontface = "bold"
          )
        
        plots$contour_plot <- p
      }
      
      # Efficiency plot
      p_efficiency <- ggplot2::ggplot(param_grid, 
                                     ggplot2::aes(x = n_samples, y = cost_efficiency, 
                                                 color = zi_modeling, group = interaction(zi_modeling, depth))) +
        ggplot2::geom_line(alpha = 0.7) +
        ggplot2::geom_point(ggplot2::aes(size = depth), alpha = 0.7) +
        ggplot2::labs(
          title = "Cost Efficiency by Sample Size",
          subtitle = "Power per million sequencing reads",
          x = "Number of Samples",
          y = "Efficiency (power/million reads)",
          color = "ZI Modeling",
          size = "Depth"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(face = "bold"),
          legend.position = "bottom"
        )
      
      plots$efficiency_plot <- p_efficiency
    }
  } else {
    plots <- NULL
  }
  
  # Format the recommendation text
  recommendation_text <- sprintf(
    "Recommended sampling strategy: %d samples with %.1f million reads per sample\n",
    primary_recommendation$n_samples, primary_recommendation$depth
  )
  
  recommendation_text <- paste0(
    recommendation_text,
    sprintf("Total sequencing effort: %.1f million reads\n", primary_recommendation$total_effort)
  )
  
  recommendation_text <- paste0(
    recommendation_text,
    sprintf("Expected power: %.2f (target: %.2f)\n", primary_recommendation$power, target_power)
  )
  
  recommendation_text <- paste0(
    recommendation_text,
    sprintf("Zero-inflation modeling: %s\n", 
           ifelse(primary_recommendation$zi_modeling == "variable", 
                  "Variable (Beta distribution)", "Fixed (uniform rate)"))
  )
  
  if (recommendation_type == "minimum_effort") {
    recommendation_text <- paste0(
      recommendation_text,
      "This strategy achieves the target power with minimum sequencing effort."
    )
  } else {
    recommendation_text <- paste0(
      recommendation_text,
      "Note: Target power cannot be achieved with the given constraints.\n",
      "This strategy provides the highest achievable power."
    )
  }
  
  # Return the results
  list(
    optimal_strategy = primary_recommendation,
    cost_efficient_strategy = cost_efficient_strategy,
    max_power_strategy = max_power_strategy,
    power_results = param_grid,
    plots = plots,
    recommendation = recommendation_text,
    parameters = list(
      target_power = target_power,
      effect_size = effect_size,
      zi_alpha = zi_alpha,
      zi_beta = zi_beta,
      total_sequencing_effort = total_sequencing_effort,
      mean_zi = mean_zi,
      sampling_zeros = sampling_zeros,
      budget_constraint = budget_constraint,
      depth_range = depth_range,
      sample_range = sample_range,
      test_variable_zi = test_variable_zi
    )
  )
}