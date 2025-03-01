#' Estimate Required Sample Size for Virome Studies
#'
#' Estimates the number of samples needed to achieve a desired statistical power
#' for detecting differences in viral abundances between groups.
#'
#' @param power Desired statistical power (default: 0.8)
#' @param effect_size Expected effect size (fold change)
#' @param n_viruses Number of viral taxa tested
#' @param alpha Significance level (default: 0.05)
#' @param sparsity Proportion of zeros in the data (default: 0.8)
#' @param dispersion Dispersion parameter for viral abundance (default: 2)
#' @param method Test method to use: "wilcoxon", "t.test", or "deseq" (default: "wilcoxon")
#'
#' @return A list containing estimated sample size and parameters
#' @export
#'
#' @examples
#' sample_size <- estimate_sample_size(power = 0.8, effect_size = 1.5, n_viruses = 100)
estimate_sample_size <- function(power = 0.8, effect_size, n_viruses,
                               alpha = 0.05, sparsity = 0.8, dispersion = 2,
                               method = "wilcoxon") {
  # Validate input parameters
  if (power <= 0 || power >= 1) {
    stop("power must be between 0 and 1")
  }
  if (effect_size <= 0) {
    stop("effect_size must be positive")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }
  
  # Define search space for sample sizes
  # Start with a reasonable range based on effect size
  min_n <- max(3, ceiling(2 / effect_size))  # Minimum n depends on effect size
  max_n <- ceiling(50 / sqrt(effect_size))   # Maximum n (smaller effect = larger max)
  
  # Initial search grid (coarse)
  sample_sizes <- seq(min_n, max_n, by = max(1, floor((max_n - min_n) / 10)))
  
  # For very small effect sizes, adjust the grid
  if (effect_size < 1.2) {
    sample_sizes <- c(sample_sizes, seq(max_n, max_n * 2, by = 5))
  }
  
  # Initialize storage for power results
  powers <- numeric(length(sample_sizes))
  
  # Calculate power for each sample size
  for (i in seq_along(sample_sizes)) {
    # Use a reduced number of simulations for efficiency
    n_sim_reduced <- min(50, max(10, floor(100 / length(sample_sizes))))
    
    # Calculate power for this sample size
    power_result <- calc_virome_power(
      n_samples = sample_sizes[i],
      effect_size = effect_size,
      n_viruses = n_viruses,
      alpha = alpha,
      sparsity = sparsity,
      dispersion = dispersion,
      method = method,
      n_sim = n_sim_reduced
    )
    
    powers[i] <- power_result$power
  }
  
  # Find smallest sample size that achieves desired power
  sufficient_power <- which(powers >= power)
  
  if (length(sufficient_power) == 0) {
    # If no sample size achieved desired power, use interpolation or extrapolation
    if (max(powers) < power) {
      # Need to extrapolate beyond our grid
      warning("Could not achieve desired power with sample sizes up to ", max(sample_sizes), 
              ". Extrapolating, but results may be unreliable.")
      
      # Fit a simple model to extrapolate
      power_model <- stats::lm(powers ~ sqrt(sample_sizes))
      
      # Calculate required n using the model (with some safety buffer)
      required_sample_size <- ceiling((power - coef(power_model)[1])^2 / coef(power_model)[2]^2 * 1.2)
      
      # Refine with a smaller number of simulations
      power_check <- calc_virome_power(
        n_samples = required_sample_size,
        effect_size = effect_size,
        n_viruses = n_viruses,
        alpha = alpha,
        sparsity = sparsity,
        dispersion = dispersion,
        method = method,
        n_sim = 20
      )
      
      achieved_power <- power_check$power
    } else {
      # Use interpolation to find the sample size
      required_sample_size <- ceiling(stats::approx(
        x = powers, 
        y = sample_sizes, 
        xout = power, 
        rule = 2
      )$y)
      
      achieved_power <- power
    }
  } else {
    # Take the smallest sample size that achieves desired power
    required_sample_size <- sample_sizes[min(sufficient_power)]
    achieved_power <- powers[min(sufficient_power)]
    
    # If we can further refine with binary search, do so
    if (min(sufficient_power) > 1) {
      lower_bound <- sample_sizes[min(sufficient_power) - 1]
      upper_bound <- required_sample_size
      
      # Do a few steps of binary search to refine
      for (i in 1:3) {
        mid <- ceiling((lower_bound + upper_bound) / 2)
        
        # Skip if the midpoint equals either bound
        if (mid == lower_bound || mid == upper_bound) {
          break
        }
        
        # Calculate power at midpoint
        mid_power <- calc_virome_power(
          n_samples = mid,
          effect_size = effect_size,
          n_viruses = n_viruses,
          alpha = alpha,
          sparsity = sparsity,
          dispersion = dispersion,
          method = method,
          n_sim = 20
        )$power
        
        if (mid_power >= power) {
          upper_bound <- mid
          achieved_power <- mid_power
          required_sample_size <- mid
        } else {
          lower_bound <- mid
        }
      }
    }
  }
  
  # Return the results
  list(
    sample_size = required_sample_size,
    achieved_power = achieved_power,
    power_by_n = data.frame(
      sample_size = sample_sizes,
      power = powers
    ),
    parameters = list(
      power = power,
      effect_size = effect_size,
      n_viruses = n_viruses,
      alpha = alpha,
      sparsity = sparsity,
      dispersion = dispersion,
      method = method
    )
  )
}