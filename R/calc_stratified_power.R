#' Calculate Power for Stratified Virome Studies
#'
#' Calculates statistical power for detecting differences in viral taxa between groups
#' in a stratified virome study design, accounting for multiple testing correction, 
#' stratification factors, and complex sampling designs.
#'
#' @param strata_sizes Vector of sample sizes for each stratum per group
#' @param effect_sizes Effect size (fold change) per stratum or single value for all strata
#' @param strata_weights Importance/weighting of each stratum (must sum to 1)
#' @param n_viruses Number of viral taxa tested
#' @param clustering_factor Intra-class correlation coefficient for clustered designs (0-1)
#' @param stratification_vars Character describing the stratification basis (informational only)
#' @param allocation_method Method for sample allocation: "proportional", "optimal", or "equal"
#' @param alpha Significance level (default: 0.05)
#' @param sparsity Proportion of zeros in the data, can be vector by stratum (default: 0.8)
#' @param dispersion Dispersion parameter for viral abundance, can be vector by stratum (default: 2)
#' @param method Test method to use: "stratified_test", "mixed_effects", or "survey_adjusted" (default: "stratified_test")
#' @param design_effect Adjustment for complex survey design, if known (default: NULL for calculation)
#' @param n_sim Number of simulations for Monte Carlo estimation (default: 100)
#'
#' @return A list containing power estimates, stratum-specific results, and parameters
#' @export
#'
#' @examples
#' # Simple stratification with two strata of different sizes
#' power_result <- calc_stratified_power(
#'   strata_sizes = c(10, 20),
#'   effect_sizes = 2.5,
#'   strata_weights = c(0.4, 0.6),
#'   n_viruses = 30,
#'   clustering_factor = 0.1,
#'   stratification_vars = "host_factors",
#'   allocation_method = "proportional"
#' )
#'
#' # Example with higher statistical power
#' power_result <- calc_stratified_power(
#'   strata_sizes = c(15, 20),
#'   effect_sizes = c(3.0, 3.5),
#'   strata_weights = c(0.4, 0.6),
#'   n_viruses = 20,
#'   clustering_factor = 0.05,
#'   sparsity = 0.6,
#'   dispersion = 1.5,
#'   stratification_vars = "geography",
#'   method = "mixed_effects"
#' )
#'
#' # Stratification with three strata and different effect sizes
#' power_result <- calc_stratified_power(
#'   strata_sizes = c(15, 15, 15),
#'   effect_sizes = c(2.5, 3.0, 3.5),
#'   strata_weights = c(0.3, 0.3, 0.4),
#'   n_viruses = 30,
#'   clustering_factor = 0.1,
#'   sparsity = 0.7,
#'   stratification_vars = "geography",
#'   allocation_method = "equal",
#'   method = "mixed_effects"
#' )
calc_stratified_power <- function(strata_sizes, 
                                effect_sizes,
                                strata_weights = NULL,
                                n_viruses,
                                clustering_factor = 0,
                                stratification_vars = "host_factors",
                                allocation_method = "proportional",
                                alpha = 0.05, 
                                sparsity = 0.8, 
                                dispersion = 2,
                                method = "stratified_test", 
                                design_effect = NULL,
                                n_sim = 100) {
  
  # Validate input parameters
  if (any(strata_sizes < 3)) {
    stop("Each stratum must have at least 3 samples per group")
  }
  
  n_strata <- length(strata_sizes)
  
  # Validate and standardize effect_sizes
  if (length(effect_sizes) == 1) {
    effect_sizes <- rep(effect_sizes, n_strata)
  } else if (length(effect_sizes) != n_strata) {
    stop("effect_sizes must be either a single value or match the number of strata")
  }
  
  if (any(effect_sizes <= 0)) {
    stop("effect_sizes must be positive")
  }
  
  # Validate and standardize strata_weights
  if (is.null(strata_weights)) {
    strata_weights <- rep(1/n_strata, n_strata)
  } else if (length(strata_weights) != n_strata) {
    stop("strata_weights must match the number of strata")
  } else if (abs(sum(strata_weights) - 1) > 0.001) {
    stop("strata_weights must sum to 1")
  }
  
  # Validate clustering_factor
  if (clustering_factor < 0 || clustering_factor > 1) {
    stop("clustering_factor must be between 0 and 1")
  }
  
  # Validate and standardize sparsity and dispersion
  if (length(sparsity) == 1) {
    sparsity <- rep(sparsity, n_strata)
  } else if (length(sparsity) != n_strata) {
    stop("sparsity must be either a single value or match the number of strata")
  }
  
  if (length(dispersion) == 1) {
    dispersion <- rep(dispersion, n_strata)
  } else if (length(dispersion) != n_strata) {
    stop("dispersion must be either a single value or match the number of strata")
  }
  
  # Validate method
  valid_methods <- c("stratified_test", "mixed_effects", "survey_adjusted")
  if (!(method %in% valid_methods)) {
    stop(paste("method must be one of:", paste(valid_methods, collapse = ", ")))
  }
  
  # Calculate design effect if not provided
  if (is.null(design_effect)) {
    # Basic design effect calculation for stratified sample
    # Assumes simple random sampling within strata
    design_effect <- 1 + (mean(strata_sizes) - 1) * clustering_factor
  }
  
  # Initialize storage for simulation results
  true_positives_overall <- numeric(n_sim)
  false_positives_overall <- numeric(n_sim)
  detected_diff_overall <- numeric(n_sim)
  
  # Create storage for stratum-specific results
  stratum_true_positives <- matrix(0, nrow = n_sim, ncol = n_strata)
  stratum_power <- numeric(n_strata)
  
  # Set proportion of differentially abundant viruses (10% is common in virome studies)
  diff_prop <- 0.1
  
  # Run simulations
  for (i in 1:n_sim) {
    # Create list to store data for each stratum
    strata_data <- list()
    
    # True differentially abundant viruses (same across strata for consistency)
    diff_taxa <- sample(n_viruses, max(1, round(n_viruses * diff_prop)))
    
    # Generate data for each stratum
    for (s in 1:n_strata) {
      # For each stratum, generate appropriate data with stratum-specific parameters
      stratum_data <- simulate_stratified_data(
        n_samples = strata_sizes[s] * 2,  # Total samples (both groups)
        n_viruses = n_viruses,
        sparsity = sparsity[s],
        dispersion = dispersion[s],
        effect_size = effect_sizes[s],
        diff_taxa = diff_taxa,
        clustering_factor = clustering_factor
      )
      
      strata_data[[s]] <- stratum_data
    }
    
    # Apply appropriate statistical test based on method
    if (method == "stratified_test") {
      # Perform stratified analysis
      results <- analyze_stratified_data(strata_data, strata_weights, alpha)
    } else if (method == "mixed_effects") {
      # Perform mixed effects analysis
      results <- analyze_mixed_effects(strata_data, alpha, clustering_factor)
    } else if (method == "survey_adjusted") {
      # Perform survey-adjusted analysis
      results <- analyze_survey_data(strata_data, strata_weights, design_effect, alpha)
    }
    
    # Extract results
    sig_taxa <- results$sig_taxa
    stratum_sig <- results$stratum_sig
    
    # Calculate true positives and false positives overall
    true_positives_overall[i] <- sum(sig_taxa %in% diff_taxa)
    false_positives_overall[i] <- sum(!sig_taxa %in% diff_taxa)
    detected_diff_overall[i] <- length(sig_taxa)
    
    # Calculate stratum-specific true positives
    for (s in 1:n_strata) {
      stratum_true_positives[i, s] <- sum(stratum_sig[[s]] %in% diff_taxa)
    }
  }
  
  # Calculate overall power
  n_true_diff <- length(diff_taxa)  # Number of truly diff taxa
  overall_power <- mean(true_positives_overall) / n_true_diff
  
  # Calculate stratum-specific power
  for (s in 1:n_strata) {
    stratum_power[s] <- mean(stratum_true_positives[, s]) / n_true_diff
  }
  
  # Calculate effective sample size
  total_samples <- sum(strata_sizes)
  effective_sample_size <- total_samples / design_effect
  
  # Calculate false discovery rate
  fdr <- sum(false_positives_overall) / max(1, sum(detected_diff_overall))
  
  # Prepare return values
  list(
    overall_power = overall_power,
    stratum_specific_power = setNames(as.list(stratum_power), paste0("stratum_", 1:n_strata)),
    effective_sample_size = effective_sample_size,
    design_effect_actual = design_effect,
    fdr = fdr,
    avg_detected = mean(detected_diff_overall),
    sim_summary = list(
      true_positives = true_positives_overall,
      false_positives = false_positives_overall,
      detected_diff = detected_diff_overall,
      stratum_detections = stratum_true_positives
    ),
    parameters = list(
      strata_sizes = strata_sizes,
      effect_sizes = effect_sizes,
      strata_weights = strata_weights,
      n_viruses = n_viruses,
      clustering_factor = clustering_factor,
      stratification_vars = stratification_vars,
      allocation_method = allocation_method,
      alpha = alpha,
      sparsity = sparsity,
      dispersion = dispersion,
      method = method,
      design_effect = design_effect,
      n_sim = n_sim
    )
  )
}

# Helper function to simulate stratified data
simulate_stratified_data <- function(n_samples, n_viruses, sparsity, dispersion, 
                                   effect_size, diff_taxa, clustering_factor) {
  # Create counts matrix (viruses x samples)
  counts <- matrix(0, nrow = n_viruses, ncol = n_samples)
  
  # Create metadata
  metadata <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    group = rep(c("A", "B"), each = n_samples/2),
    stratum_id = rep(1, n_samples)  # Single stratum for this function
  )
  
  # Define baseline abundance (log-normal distribution is common in microbiome)
  baseline_abundance <- rlnorm(n_viruses, meanlog = 1, sdlog = 2)
  
  # Define group membership
  group_a <- metadata$group == "A"
  group_b <- metadata$group == "B"
  
  # Simulate counts for each viral taxon
  for (v in 1:n_viruses) {
    # Determine if this virus is differentially abundant
    is_diff <- v %in% diff_taxa
    
    # Base abundance for this virus
    lambda_a <- baseline_abundance[v]
    
    # Effect size only applied to differential taxa in group B
    lambda_b <- if (is_diff) lambda_a * effect_size else lambda_a
    
    # Apply dispersion (negative binomial is appropriate for overdispersed count data)
    size_param <- 1/dispersion  # Size parameter for negative binomial
    
    # Generate counts for group A
    counts_a <- rnbinom(sum(group_a), mu = lambda_a, size = size_param)
    
    # Generate counts for group B
    counts_b <- rnbinom(sum(group_b), mu = lambda_b, size = size_param)
    
    # Apply clustering effect if specified
    if (clustering_factor > 0) {
      # Simple approach: add correlated noise based on clustering factor
      cluster_effect_a <- rnorm(1, 0, lambda_a * clustering_factor)
      cluster_effect_b <- rnorm(1, 0, lambda_b * clustering_factor)
      
      counts_a <- pmax(0, counts_a + cluster_effect_a)
      counts_b <- pmax(0, counts_b + cluster_effect_b)
    }
    
    # Apply sparsity (zeros)
    zero_indices_a <- sample(which(group_a), size = round(sum(group_a) * sparsity))
    zero_indices_b <- sample(which(group_b), size = round(sum(group_b) * sparsity))
    
    # Create the final count vector
    counts_combined <- numeric(n_samples)
    counts_combined[group_a] <- counts_a
    counts_combined[group_b] <- counts_b
    
    # Apply zeros
    counts_combined[zero_indices_a] <- 0
    counts_combined[zero_indices_b] <- 0
    
    # Store in counts matrix
    counts[v, ] <- counts_combined
  }
  
  # Return the simulated data
  list(
    counts = counts,
    metadata = metadata,
    diff_taxa = diff_taxa
  )
}

# Helper function to analyze stratified data
analyze_stratified_data <- function(strata_data, strata_weights, alpha) {
  n_strata <- length(strata_data)
  n_viruses <- nrow(strata_data[[1]]$counts)
  
  # Initialize p-values
  p_values <- numeric(n_viruses)
  stratum_sig <- list()
  
  # Analyze each virus across all strata
  for (v in 1:n_viruses) {
    # Initialize weighted test statistic components
    weighted_stats <- 0
    weighted_vars <- 0
    
    # Store stratum-specific p-values
    stratum_p_values <- numeric(n_strata)
    
    # Analyze each stratum
    for (s in 1:n_strata) {
      stratum <- strata_data[[s]]
      counts <- stratum$counts
      metadata <- stratum$metadata
      
      # Extract data for this virus in this stratum
      group_a <- counts[v, metadata$group == "A"]
      group_b <- counts[v, metadata$group == "B"]
      
      # Skip taxa with all zeros (common in virome data)
      if (all(group_a == 0) && all(group_b == 0)) {
        stratum_p_values[s] <- 1
        next
      }
      
      # Simple Wilcoxon test for demonstration
      # In practice, would use more sophisticated methods accounting for sparsity
      test_result <- tryCatch({
        suppressWarnings(
          stats::wilcox.test(group_a, group_b, exact = FALSE)
        )
      }, error = function(e) {
        return(list(p.value = 1))
      })
      
      stratum_p_values[s] <- test_result$p.value
      
      # For weighted combination, convert to Z-score
      z_score <- stats::qnorm(1 - test_result$p.value/2) * sign(mean(group_b) - mean(group_a))
      if (is.nan(z_score)) z_score <- 0
      
      # Apply stratum weight
      weighted_stats <- weighted_stats + strata_weights[s] * z_score
      weighted_vars <- weighted_vars + strata_weights[s]^2
    }
    
    # Combined p-value using weighted Z-method
    # This is a simple demonstration; more sophisticated methods exist
    combined_z <- weighted_stats / sqrt(weighted_vars)
    p_values[v] <- 2 * (1 - stats::pnorm(abs(combined_z)))
    
    # Store stratum-specific significant taxa
    for (s in 1:n_strata) {
      if (is.null(stratum_sig[[s]])) stratum_sig[[s]] <- numeric(0)
      if (stratum_p_values[s] < alpha) {
        stratum_sig[[s]] <- c(stratum_sig[[s]], v)
      }
    }
  }
  
  # Apply multiple testing correction
  adjusted_p <- stats::p.adjust(p_values, method = "BH")  # Benjamini-Hochberg
  
  # Determine significant taxa
  sig_taxa <- which(adjusted_p < alpha)
  
  # Return results
  list(
    sig_taxa = sig_taxa,
    stratum_sig = stratum_sig,
    p_values = p_values,
    adjusted_p = adjusted_p
  )
}

# Helper function for mixed effects analysis
analyze_mixed_effects <- function(strata_data, alpha, clustering_factor) {
  # In a real implementation, this would use lme4 or other mixed-effects packages
  # Here we'll simulate the results for demonstration
  
  n_strata <- length(strata_data)
  n_viruses <- nrow(strata_data[[1]]$counts)
  
  # Initialize results
  p_values <- numeric(n_viruses)
  stratum_sig <- list()
  for (s in 1:n_strata) stratum_sig[[s]] <- numeric(0)
  
  # Combine data across strata
  combined_counts <- matrix(0, nrow = n_viruses, ncol = 0)
  combined_metadata <- data.frame(sample_id = character(0), group = character(0), stratum_id = numeric(0))
  
  for (s in 1:n_strata) {
    stratum <- strata_data[[s]]
    combined_counts <- cbind(combined_counts, stratum$counts)
    
    stratum_metadata <- stratum$metadata
    stratum_metadata$stratum_id <- s
    combined_metadata <- rbind(combined_metadata, stratum_metadata)
  }
  
  # For each virus, perform mixed-effects analysis
  for (v in 1:n_viruses) {
    # Extract counts for this virus
    virus_counts <- combined_counts[v, ]
    
    # In real implementation: fit mixed model with stratum as random effect
    # lmer(counts ~ group + (1|stratum_id), data=...)
    
    # For demonstration, simulate p-value based on data properties
    mean_diff <- mean(virus_counts[combined_metadata$group == "B"]) - 
                mean(virus_counts[combined_metadata$group == "A"])
    std_err <- sd(virus_counts) / sqrt(length(virus_counts))
    
    # Adjust for clustering effect
    std_err <- std_err * sqrt(1 + clustering_factor)
    
    # Calculate test statistic and p-value
    t_stat <- mean_diff / std_err
    df <- length(virus_counts) - n_strata  # Adjust df for stratification
    p_values[v] <- 2 * pt(-abs(t_stat), df)
    
    # For demonstration of stratum-specific effects:
    # Also assess each stratum individually
    for (s in 1:n_strata) {
      stratum_indices <- combined_metadata$stratum_id == s
      stratum_a <- virus_counts[stratum_indices & combined_metadata$group == "A"]
      stratum_b <- virus_counts[stratum_indices & combined_metadata$group == "B"]
      
      if (length(stratum_a) > 0 && length(stratum_b) > 0) {
        stratum_p <- tryCatch({
          t.test(stratum_a, stratum_b)$p.value
        }, error = function(e) 1)
        
        # Check if stratum_p is NA or NULL before comparison
        if (!is.null(stratum_p) && !is.na(stratum_p) && stratum_p < alpha) {
          stratum_sig[[s]] <- c(stratum_sig[[s]], v)
        }
      }
    }
  }
  
  # Apply multiple testing correction
  adjusted_p <- stats::p.adjust(p_values, method = "BH")
  
  # Determine significant taxa
  sig_taxa <- which(adjusted_p < alpha)
  
  list(
    sig_taxa = sig_taxa,
    stratum_sig = stratum_sig,
    p_values = p_values,
    adjusted_p = adjusted_p
  )
}

# Helper function for survey-adjusted analysis
analyze_survey_data <- function(strata_data, strata_weights, design_effect, alpha) {
  # In real implementation, would use survey package
  # Here we'll simulate the approach
  
  n_strata <- length(strata_data)
  n_viruses <- nrow(strata_data[[1]]$counts)
  
  # Initialize p-values
  p_values <- numeric(n_viruses)
  stratum_sig <- list()
  for (s in 1:n_strata) stratum_sig[[s]] <- numeric(0)
  
  # Analyze each virus across all strata
  for (v in 1:n_viruses) {
    # Initialize weighted estimates
    weighted_diff <- 0
    weighted_var <- 0
    
    # Analyze each stratum
    for (s in 1:n_strata) {
      stratum <- strata_data[[s]]
      counts <- stratum$counts
      metadata <- stratum$metadata
      
      # Extract data for this virus in this stratum
      group_a <- counts[v, metadata$group == "A"]
      group_b <- counts[v, metadata$group == "B"]
      
      # Calculate stratum-specific difference
      mean_a <- mean(group_a)
      mean_b <- mean(group_b)
      diff <- mean_b - mean_a
      
      # Estimate variance of the difference
      var_a <- var(group_a) / length(group_a)
      var_b <- var(group_b) / length(group_b)
      var_diff <- var_a + var_b
      
      # Apply design effect to variance
      var_diff <- var_diff * design_effect
      
      # Apply stratum weight
      weighted_diff <- weighted_diff + strata_weights[s] * diff
      weighted_var <- weighted_var + strata_weights[s]^2 * var_diff
      
      # Also check for significance within this stratum
      if (length(group_a) > 0 && length(group_b) > 0) {
        t_stat <- diff / sqrt(var_diff)
        df <- length(group_a) + length(group_b) - 2
        p_stratum <- 2 * pt(-abs(t_stat), df)
        
        if (p_stratum < alpha) {
          stratum_sig[[s]] <- c(stratum_sig[[s]], v)
        }
      }
    }
    
    # Calculate overall test statistic
    t_stat <- weighted_diff / sqrt(weighted_var)
    
    # Calculate p-value (using normal approximation for simplicity)
    p_values[v] <- 2 * (1 - pnorm(abs(t_stat)))
  }
  
  # Apply multiple testing correction
  adjusted_p <- stats::p.adjust(p_values, method = "BH")
  
  # Determine significant taxa
  sig_taxa <- which(adjusted_p < alpha)
  
  list(
    sig_taxa = sig_taxa,
    stratum_sig = stratum_sig,
    p_values = p_values,
    adjusted_p = adjusted_p
  )
}