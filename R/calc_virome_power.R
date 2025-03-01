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
  
  # Run simulations with error handling
  for (i in 1:n_sim) {
    # Use try to catch and handle any errors during simulation
    sim_data <- try({
      # Make sure parameters are reasonable
      adjusted_n_samples <- max(3, n_samples)  # Ensure minimum 3 samples per group
      adjusted_dispersion <- max(0.1, dispersion)  # Ensure reasonable dispersion
      
      # Simulate virome data
      simulate_virome_data(
        n_samples = adjusted_n_samples * 2,  # Total samples (both groups)
        n_viruses = n_viruses,
        sparsity = sparsity,
        dispersion = adjusted_dispersion,
        effect_size = effect_size
      )
    }, silent = TRUE)
    
    # If simulation failed, create minimal valid data
    if (inherits(sim_data, "try-error")) {
      warning("Simulation failed, using minimal data structure")
      sim_data <- list(
        counts = matrix(rpois(n_viruses * n_samples * 2, lambda = 1), 
                      nrow = n_viruses, ncol = n_samples * 2),
        metadata = data.frame(
          sample_id = paste0("sample_", 1:(n_samples * 2)),
          group = rep(c("A", "B"), each = n_samples)
        ),
        diff_taxa = sample(n_viruses, max(1, round(n_viruses * 0.1)))
      )
    }
    
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
        # Only run test if we have sufficient non-zero data points
        n_nonzero_a <- sum(group_a > 0)
        n_nonzero_b <- sum(group_b > 0)
        
        if (n_nonzero_a >= 2 && n_nonzero_b >= 2) {
          # Enough data for a reasonable test
          test_result <- tryCatch({
            # Use exact=FALSE to avoid warnings about ties
            suppressWarnings(
              stats::wilcox.test(group_a, group_b, exact = FALSE)
            )
          }, error = function(e) {
            # Handle other issues
            return(list(p.value = 1))
          })
          p_values[j] <- test_result$p.value
        } else {
          # Not enough data points for a reliable test
          p_values[j] <- 1
        }
        
      } else if (method == "t.test") {
        # t-test (parametric)
        # Only run if we have enough data points and variance
        if (length(group_a) >= 3 && length(group_b) >= 3 && 
            var(group_a) > 0 && var(group_b) > 0) {
          test_result <- tryCatch({
            suppressWarnings(
              stats::t.test(group_a, group_b, var.equal = FALSE)
            )
          }, error = function(e) {
            return(list(p.value = 1))
          })
          p_values[j] <- test_result$p.value
        } else {
          # Not enough data or variance for reliable test
          p_values[j] <- 1
        }
        
      } else if (method == "deseq") {
        # Simplified DESeq-like negative binomial test with robust error handling
        mean_a <- mean(group_a)
        mean_b <- mean(group_b)
        
        # Skip if means are too low or equal
        if (mean_a < 0.1 || mean_b < 0.1 || abs(mean_a - mean_b) < 1e-6) {
          p_values[j] <- 1
          next
        }
        
        # Calculate variances, handling single sample case
        var_a <- ifelse(length(group_a) > 1, var(group_a), mean_a)
        var_b <- ifelse(length(group_b) > 1, var(group_b), mean_b)
        
        # Ensure variances are at least equal to the means (overdispersion)
        var_a <- max(var_a, mean_a)
        var_b <- max(var_b, mean_b)
        
        # Estimate dispersion with robust bound checking
        disp_a <- min(10, max(0.01, (var_a - mean_a) / (mean_a^2 + 1e-6)))
        disp_b <- min(10, max(0.01, (var_b - mean_b) / (mean_b^2 + 1e-6)))
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
  
  # Calculate power (proportion of true positives detected) with safety checks
  n_true_diff <- max(1, length(sim_data$diff_taxa))  # Number of truly diff taxa (avoid division by 0)
  power <- mean(true_positives) / n_true_diff
  
  # Make sure power is in valid range [0, 1]
  power <- max(0, min(1, power))
  
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