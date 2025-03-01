#' Calculate Bayesian Power for Virome Studies
#'
#' This function performs Bayesian power analysis for virome studies comparing two groups.
#' It uses simulation-based methods to estimate the probability of detecting true effects
#' under a Bayesian framework, which can be more appropriate for sparse, high-dimensional
#' virome data than traditional frequentist approaches.
#'
#' @param n_samples Number of samples per group
#' @param effect_size Expected effect size (fold change between groups)
#' @param n_viruses Number of viral taxa to simulate
#' @param sparsity Proportion of zeros in the count matrix (default: 0.7)
#' @param dispersion Dispersion parameter for negative binomial distribution (default: 1.2)
#' @param prior_strength Prior strength (concentration) parameter (default: 2)
#' @param credible_interval Width of the credible interval for decision threshold (default: 0.95)
#' @param fold_change_threshold Minimum fold change to consider biologically relevant (default: 1.5)
#' @param posterior_prob_threshold Posterior probability threshold for significance (default: 0.95)
#' @param n_sim Number of simulation iterations (default: 100)
#' @param seed Random seed for reproducibility (default: NULL)
#'
#' @return A list containing:
#' \itemize{
#'   \item power: Estimated Bayesian power (probability of detecting true effect)
#'   \item expected_discoveries: Expected number of discoveries (true positives)
#'   \item false_discovery_proportion: Expected proportion of false discoveries
#'   \item bayes_factor_summary: Summary statistics of Bayes factors from simulations
#'   \item posterior_prob_summary: Summary of posterior probabilities from simulations
#'   \item simulation_results: Detailed simulation results for each iteration
#'   \item parameters: List of input parameters
#' }
#'
#' @examples
#' # Calculate Bayesian power for a virome study
#' power_result <- calc_bayesian_power(
#'   n_samples = 20,         # 20 samples per group
#'   effect_size = 3.0,      # 3-fold difference between groups
#'   n_viruses = 50,         # 50 viral taxa
#'   sparsity = 0.6,         # 60% zeros in the data
#'   prior_strength = 2.5,   # Moderate prior strength
#'   n_sim = 50              # 50 simulation iterations
#' )
#'
#' # Print the Bayesian power estimate
#' print(paste("Bayesian power:", round(power_result$power * 100, 1), "%"))
#'
#' # Print the expected number of discoveries
#' print(paste("Expected discoveries:", round(power_result$expected_discoveries, 1)))
#'
#' # Generate a Bayesian power report
#' report_path <- generate_bayesian_power_report(
#'   bayesian_power_results = power_result,
#'   output_file = "bayesian_power_report.html",
#'   title = "Bayesian Power Analysis for Virome Study"
#' )
#'
#' @export
calc_bayesian_power <- function(n_samples,
                               effect_size,
                               n_viruses,
                               sparsity = 0.7,
                               dispersion = 1.2,
                               prior_strength = 2,
                               credible_interval = 0.95,
                               fold_change_threshold = 1.5,
                               posterior_prob_threshold = 0.95,
                               n_sim = 100,
                               seed = NULL) {
  
  # Validate inputs
  if (n_samples < 3) {
    stop("n_samples must be at least 3 per group")
  }
  
  if (effect_size <= 1) {
    stop("effect_size must be greater than 1")
  }
  
  if (n_viruses < 1) {
    stop("n_viruses must be at least 1")
  }
  
  if (sparsity < 0 || sparsity > 0.99) {
    stop("sparsity must be between 0 and 0.99")
  }
  
  if (dispersion <= 0) {
    stop("dispersion must be positive")
  }
  
  if (prior_strength <= 0) {
    stop("prior_strength must be positive")
  }
  
  if (credible_interval <= 0 || credible_interval >= 1) {
    stop("credible_interval must be between 0 and 1")
  }
  
  if (posterior_prob_threshold <= 0 || posterior_prob_threshold >= 1) {
    stop("posterior_prob_threshold must be between 0 and 1")
  }
  
  if (n_sim < 10) {
    stop("n_sim should be at least 10 for reliable results")
  }
  
  # Set random seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Store parameters for return
  parameters <- list(
    n_samples = n_samples,
    effect_size = effect_size,
    n_viruses = n_viruses,
    sparsity = sparsity,
    dispersion = dispersion,
    prior_strength = prior_strength,
    credible_interval = credible_interval,
    fold_change_threshold = fold_change_threshold,
    posterior_prob_threshold = posterior_prob_threshold,
    n_sim = n_sim,
    seed = seed
  )
  
  # Define number of viruses with true effects (50% by default)
  n_true_effects <- floor(n_viruses / 2)
  effect_indices <- sample(1:n_viruses, n_true_effects)
  
  # Initialize results storage
  simulation_results <- list()
  
  # Main simulation loop
  for (i in 1:n_sim) {
    # Generate simulated data
    sim_data <- simulate_virome_data(
      n_samples = n_samples * 2,  # Total samples across both groups
      n_viruses = n_viruses,
      group_effect = TRUE,
      effect_size = effect_size,
      effect_indices = effect_indices,
      sparsity = sparsity,
      dispersion = dispersion
    )
    
    # Perform Bayesian analysis
    bayes_result <- analyze_bayesian(
      counts = sim_data$counts,
      groups = sim_data$groups,
      prior_strength = prior_strength,
      credible_interval = credible_interval
    )
    
    # Evaluate results against true effects
    eval_result <- evaluate_bayesian_results(
      bayes_result = bayes_result,
      effect_indices = effect_indices,
      posterior_prob_threshold = posterior_prob_threshold,
      fold_change_threshold = fold_change_threshold
    )
    
    # Store results for this iteration
    simulation_results[[i]] <- list(
      true_positives = eval_result$true_positives,
      false_positives = eval_result$false_positives,
      true_negatives = eval_result$true_negatives,
      false_negatives = eval_result$false_negatives,
      bayes_factors = bayes_result$bayes_factors,
      posterior_probs = bayes_result$posterior_probs,
      effect_sizes = bayes_result$effect_sizes,
      credible_intervals = bayes_result$credible_intervals
    )
  }
  
  # Calculate summary statistics
  true_positives <- sapply(simulation_results, function(x) length(x$true_positives))
  false_positives <- sapply(simulation_results, function(x) length(x$false_positives))
  true_negatives <- sapply(simulation_results, function(x) length(x$true_negatives))
  false_negatives <- sapply(simulation_results, function(x) length(x$false_negatives))
  
  # Calculate power metrics
  power <- mean(true_positives / n_true_effects)
  expected_discoveries <- mean(true_positives + false_positives)
  false_discovery_proportion <- mean(false_positives / pmax(true_positives + false_positives, 1))
  
  # Extract Bayes factors and posterior probabilities for summary
  bayes_factors <- unlist(lapply(simulation_results, function(x) x$bayes_factors))
  posterior_probs <- unlist(lapply(simulation_results, function(x) x$posterior_probs))
  
  # Summarize Bayes factors and posterior probabilities
  bayes_factor_summary <- list(
    mean = mean(bayes_factors),
    median = median(bayes_factors),
    q25 = quantile(bayes_factors, 0.25),
    q75 = quantile(bayes_factors, 0.75)
  )
  
  posterior_prob_summary <- list(
    mean = mean(posterior_probs),
    median = median(posterior_probs),
    q25 = quantile(posterior_probs, 0.25),
    q75 = quantile(posterior_probs, 0.75)
  )
  
  # Create summary of simulation results
  sim_summary <- list(
    true_positives = mean(true_positives),
    false_positives = mean(false_positives),
    true_negatives = mean(true_negatives),
    false_negatives = mean(false_negatives)
  )
  
  # Prepare results
  results <- list(
    power = power,
    expected_discoveries = expected_discoveries,
    false_discovery_proportion = false_discovery_proportion,
    bayes_factor_summary = bayes_factor_summary,
    posterior_prob_summary = posterior_prob_summary,
    simulation_results = simulation_results,
    sim_summary = sim_summary,
    parameters = parameters
  )
  
  return(results)
}

#' Analyze Virome Data Using Bayesian Methods
#'
#' Internal function to perform Bayesian analysis on virome count data
#'
#' @param counts Matrix of virome counts
#' @param groups Vector of group assignments
#' @param prior_strength Prior strength parameter
#' @param credible_interval Width of credible interval
#'
#' @return List with Bayesian analysis results
#'
#' @keywords internal
analyze_bayesian <- function(counts, groups, prior_strength, credible_interval) {
  # Get dimensions
  n_viruses <- nrow(counts)
  
  # Initialize result containers
  bayes_factors <- numeric(n_viruses)
  posterior_probs <- numeric(n_viruses)
  effect_sizes <- numeric(n_viruses)
  credible_intervals <- matrix(0, nrow = n_viruses, ncol = 2)
  
  # Split data by group
  group_indices <- split(1:length(groups), groups)
  counts_g1 <- counts[, group_indices[[1]], drop = FALSE]
  counts_g2 <- counts[, group_indices[[2]], drop = FALSE]
  
  # For each virus, perform Bayesian analysis
  for (i in 1:n_viruses) {
    # Extract counts for this virus
    y1 <- counts_g1[i, ]
    y2 <- counts_g2[i, ]
    
    # Calculate summary statistics
    n1 <- length(y1)
    n2 <- length(y2)
    sum_y1 <- sum(y1)
    sum_y2 <- sum(y2)
    
    # Add prior pseudocounts (Dirichlet-multinomial model)
    alpha1 <- sum_y1 + prior_strength
    alpha2 <- sum_y2 + prior_strength
    
    # Calculate posterior for null hypothesis (no difference)
    log_prob_null <- lbeta(sum_y1 + sum_y2 + 2 * prior_strength, 
                          n1 + n2 - sum_y1 - sum_y2 + 2 * prior_strength) -
                     lbeta(2 * prior_strength, 2 * prior_strength)
    
    # Calculate posterior for alternative hypothesis (difference exists)
    log_prob_alt <- lbeta(sum_y1 + prior_strength, n1 - sum_y1 + prior_strength) +
                   lbeta(sum_y2 + prior_strength, n2 - sum_y2 + prior_strength) -
                   lbeta(prior_strength, prior_strength) -
                   lbeta(prior_strength, prior_strength)
    
    # Calculate Bayes factor (exp to convert from log scale)
    bf <- exp(log_prob_alt - log_prob_null)
    bayes_factors[i] <- bf
    
    # Calculate posterior probability of alternative hypothesis
    # Using equal prior probabilities for null and alternative
    posterior_probs[i] <- bf / (1 + bf)
    
    # Calculate posterior effect size (fold change) using posterior means
    # Add small pseudocount to avoid division by zero
    mean1 <- (sum_y1 + prior_strength) / (n1 + 2 * prior_strength)
    mean2 <- (sum_y2 + prior_strength) / (n2 + 2 * prior_strength)
    effect_sizes[i] <- mean2 / max(mean1, 1e-6)
    
    # Calculate credible interval for effect size using beta distribution properties
    # This is a simplification; for real implementation, consider MCMC
    alpha1_post <- sum_y1 + prior_strength
    beta1_post <- n1 - sum_y1 + prior_strength
    alpha2_post <- sum_y2 + prior_strength
    beta2_post <- n2 - sum_y2 + prior_strength
    
    # Simplified credible interval calculation
    # In practice, would use posterior samples for more accurate intervals
    lower_quantile <- (1 - credible_interval) / 2
    upper_quantile <- 1 - lower_quantile
    
    # Approximate credible interval for effect size
    ci_lower <- qbeta(lower_quantile, alpha2_post, beta2_post) / 
               qbeta(upper_quantile, alpha1_post, beta1_post)
    ci_upper <- qbeta(upper_quantile, alpha2_post, beta2_post) / 
               qbeta(lower_quantile, alpha1_post, beta1_post)
    
    credible_intervals[i, ] <- c(ci_lower, ci_upper)
  }
  
  return(list(
    bayes_factors = bayes_factors,
    posterior_probs = posterior_probs,
    effect_sizes = effect_sizes,
    credible_intervals = credible_intervals
  ))
}

#' Evaluate Bayesian Analysis Results
#'
#' Internal function to evaluate Bayesian analysis results against ground truth
#'
#' @param bayes_result Results from Bayesian analysis
#' @param effect_indices Indices of viruses with true effects
#' @param posterior_prob_threshold Threshold for posterior probability
#' @param fold_change_threshold Threshold for fold change
#'
#' @return List with evaluation metrics
#'
#' @keywords internal
evaluate_bayesian_results <- function(bayes_result, effect_indices, 
                                     posterior_prob_threshold, fold_change_threshold) {
  # Identify significant findings
  # Criteria: high posterior probability AND biologically relevant effect size
  significant_indices <- which(bayes_result$posterior_probs >= posterior_prob_threshold & 
                              abs(log2(bayes_result$effect_sizes)) >= log2(fold_change_threshold))
  
  # Create logical vector for true effects
  n_viruses <- length(bayes_result$posterior_probs)
  true_effect <- rep(FALSE, n_viruses)
  true_effect[effect_indices] <- TRUE
  
  # Calculate performance metrics
  true_positives <- intersect(significant_indices, effect_indices)
  false_positives <- setdiff(significant_indices, effect_indices)
  true_negatives <- setdiff(which(!true_effect), significant_indices)
  false_negatives <- setdiff(effect_indices, significant_indices)
  
  return(list(
    true_positives = true_positives,
    false_positives = false_positives,
    true_negatives = true_negatives,
    false_negatives = false_negatives
  ))
}