#' Calculate Power for Virome Studies
#'
#' Calculates statistical power for detecting differences in viral taxa between groups
#' in a virome study, accounting for multiple testing correction and sparsity.
#'
#' @param n_samples Number of samples per group
#' @param effect_size Expected effect size (fold change)
#' @param n_viruses Number of viral taxa tested
#' @param alpha Significance level (default: 0.05)
#' @param sparsity Proportion of zeros in the data (default: 0.8)
#' @param dispersion Dispersion parameter for viral abundance (default: 2)
#' @param method Test method to use: "wilcoxon", "t.test", or "deseq" (default: "wilcoxon")
#' @param n_sim Number of simulations for Monte Carlo estimation (default: 100)
#'
#' @return A list containing power estimates, simulation results, and parameters
#' @export
#'
#' @examples
#' power_result <- calc_virome_power(n_samples = 10, effect_size = 1.5, n_viruses = 100)
calc_virome_power <- function(n_samples, effect_size, n_viruses,
                             alpha = 0.05, sparsity = 0.8, dispersion = 2,
                             method = "wilcoxon", n_sim = 100) {
  # Validate input parameters
  if (n_samples < 3) {
    stop("n_samples must be at least 3 per group")
  }
  if (effect_size <= 0) {
    stop("effect_size must be positive")
  }
  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }
  if (!(method %in% c("wilcoxon", "t.test", "deseq"))) {
    stop("method must be one of: wilcoxon, t.test, deseq")
  }
  
  # Initialize result storage
  detected_diff <- numeric(n_sim)
  true_positives <- numeric(n_sim)
  false_positives <- numeric(n_sim)
  
  # Run simulations
  for (i in 1:n_sim) {
    # Simulate virome data
    sim_data <- simulate_virome_data(
      n_samples = n_samples * 2,  # Total samples (both groups)
      n_viruses = n_viruses,
      sparsity = sparsity,
      dispersion = dispersion,
      effect_size = effect_size
    )
    
    # Extract data
    counts <- sim_data$counts
    metadata <- sim_data$metadata
    diff_taxa <- sim_data$diff_taxa  # True differentially abundant taxa
    
    # Apply specified statistical test
    p_values <- numeric(n_viruses)
    for (j in 1:n_viruses) {
      group_a <- counts[j, metadata$group == "A"]
      group_b <- counts[j, metadata$group == "B"]
      
      # Skip taxa with all zeros (common in virome data)
      if (all(group_a == 0) && all(group_b == 0)) {
        p_values[j] <- 1
        next
      }
      
      # Apply the selected statistical test
      if (method == "wilcoxon") {
        # Wilcoxon rank-sum test (non-parametric, robust to non-normality)
        test_result <- tryCatch({
          stats::wilcox.test(group_a, group_b)
        }, error = function(e) {
          # Handle ties or other issues
          return(list(p.value = 1))
        })
        p_values[j] <- test_result$p.value
        
      } else if (method == "t.test") {
        # t-test (parametric)
        test_result <- tryCatch({
          stats::t.test(group_a, group_b, var.equal = FALSE)
        }, error = function(e) {
          return(list(p.value = 1))
        })
        p_values[j] <- test_result$p.value
        
      } else if (method == "deseq") {
        # Simplified DESeq-like negative binomial test
        mean_a <- mean(group_a)
        mean_b <- mean(group_b)
        var_a <- var(group_a)
        var_b <- var(group_b)
        
        # Skip if means are zero
        if (mean_a == 0 || mean_b == 0) {
          p_values[j] <- 1
          next
        }
        
        # Estimate dispersion
        disp_a <- max(1e-6, (var_a - mean_a) / (mean_a^2))
        disp_b <- max(1e-6, (var_b - mean_b) / (mean_b^2))
        disp <- mean(c(disp_a, disp_b))
        
        # Simple likelihood ratio test
        log_fc <- log2(mean_b / mean_a)
        se <- sqrt(disp * (1/length(group_a) + 1/length(group_b)))
        z <- log_fc / se
        p_values[j] <- 2 * stats::pnorm(-abs(z))
      }
    }
    
    # Apply multiple testing correction
    adjusted_p <- stats::p.adjust(p_values, method = "BH")  # Benjamini-Hochberg
    
    # Determine significant taxa
    sig_taxa <- which(adjusted_p < alpha)
    
    # Calculate true positives and false positives
    true_positives[i] <- sum(sig_taxa %in% diff_taxa)
    false_positives[i] <- sum(!sig_taxa %in% diff_taxa)
    
    # Store number of detected differences
    detected_diff[i] <- length(sig_taxa)
  }
  
  # Calculate power (proportion of true positives detected)
  n_true_diff <- length(sim_data$diff_taxa)  # Number of truly diff taxa
  power <- mean(true_positives) / n_true_diff
  
  # Calculate false discovery rate
  fdr <- sum(false_positives) / max(1, sum(detected_diff))
  
  # Return results
  list(
    power = power,
    fdr = fdr,
    avg_detected = mean(detected_diff),
    sim_summary = list(
      true_positives = true_positives,
      false_positives = false_positives,
      detected_diff = detected_diff
    ),
    parameters = list(
      n_samples = n_samples,
      effect_size = effect_size,
      n_viruses = n_viruses,
      alpha = alpha,
      sparsity = sparsity,
      dispersion = dispersion,
      method = method,
      n_sim = n_sim
    )
  )
}