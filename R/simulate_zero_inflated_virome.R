#' Simulate Zero-Inflated Virome Data
#'
#' Generates simulated virome data using a zero-inflated negative binomial (ZINB) model,
#' which is particularly suitable for viral metagenomic data with excess zeros beyond
#' what would be expected from a standard negative binomial distribution. The ZINB model
#' combines a point mass at zero with a negative binomial distribution, allowing for 
#' modeling of structural zeros (true absence) and sampling zeros (detection failures).
#'
#' @param n_samples Number of samples to simulate
#' @param n_viruses Number of viral taxa to simulate
#' @param structural_zeros Proportion of structural zeros (true absence) in the data (default: 0.7)
#' @param sampling_zeros Proportion of sampling zeros (detection failures) among non-structural zeros (default: 0.2)
#' @param dispersion Dispersion parameter for negative binomial distribution (default: 1.5)
#' @param effect_size Effect size for group differences (default: 2.0)
#' @param zero_inflation_difference Whether the zero-inflation probability differs between groups (default: TRUE)
#' @param groups Vector of group labels for samples (default: equal groups)
#' @param variable_zi_rates Logical; whether to use varying zero inflation rates across viral taxa (default: FALSE)
#' @param zi_alpha Shape parameter for beta distribution to generate variable structural zero rates (default: 2)
#' @param zi_beta Shape parameter for beta distribution to generate variable structural zero rates (default: 3)
#'
#' @return A list containing simulated count matrix and sample metadata, with additional
#'         information about the zero-inflation structure
#' @export
#'
#' @examples
#' # Basic zero-inflated simulation
#' sim_data <- simulate_zero_inflated_virome(n_samples = 30, n_viruses = 100)
#' 
#' # Simulation with varying zero-inflation rates across viral taxa
#' sim_data_var <- simulate_zero_inflated_virome(
#'   n_samples = 30, 
#'   n_viruses = 100,
#'   variable_zi_rates = TRUE,
#'   zi_alpha = 2,
#'   zi_beta = 5  # More right-skewed distribution of zi rates
#' )
simulate_zero_inflated_virome <- function(n_samples, 
                                         n_viruses, 
                                         structural_zeros = 0.7, 
                                         sampling_zeros = 0.2,
                                         dispersion = 1.5, 
                                         effect_size = 2.0,
                                         zero_inflation_difference = TRUE,
                                         groups = NULL,
                                         variable_zi_rates = FALSE,
                                         zi_alpha = 2,
                                         zi_beta = 3) {
  # Create group labels if not provided
  if (is.null(groups)) {
    # Default: Half samples in group A, half in group B
    if (n_samples %% 2 != 0) {
      # Handle odd number silently by ensuring even split
      n_group_a <- floor(n_samples/2)
      n_group_b <- n_samples - n_group_a
      groups <- c(rep("A", n_group_a), rep("B", n_group_b))
    } else {
      groups <- rep(c("A", "B"), each = n_samples/2)
    }
  } else {
    # Check if provided groups match n_samples
    if (length(groups) != n_samples) {
      stop("Length of groups must match n_samples")
    }
  }
  
  # Create sample metadata
  metadata <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    group = groups
  )
  
  # Generate abundance distribution following power law (common in virome data)
  # Use exponential distribution for abundance means
  abundance_means <- rexp(n_viruses, rate = 0.1)
  abundance_means <- abundance_means / sum(abundance_means) * 1e6  # Scale to reasonable counts
  
  # Make a subset of taxa differentially abundant between groups
  n_diff <- round(n_viruses * 0.15)  # 15% of taxa are differentially abundant
  diff_idx <- sample(n_viruses, n_diff)
  
  # Create matrix to store counts and structural zero indicators
  counts <- matrix(0, nrow = n_viruses, ncol = n_samples)
  structural_zero_matrix <- matrix(FALSE, nrow = n_viruses, ncol = n_samples)
  sampling_zero_matrix <- matrix(FALSE, nrow = n_viruses, ncol = n_samples)
  
  colnames(counts) <- metadata$sample_id
  rownames(counts) <- paste0("virus_", 1:n_viruses)
  
  # Determine which group indices are group A and B
  group_a_indices <- which(groups == "A")
  group_b_indices <- which(groups == "B")
  
  # Generate zero-inflation probabilities
  # For each virus, determine structural zero probability (can vary by group)
  structural_zero_probs <- matrix(0, nrow = n_viruses, ncol = 2)  # Col 1 = Group A, Col 2 = Group B
  
  # Track individual virus zero-inflation rates for output
  virus_specific_zi_rates <- numeric(n_viruses)
  
  # If variable_zi_rates is TRUE, generate virus-specific zero-inflation rates
  if (variable_zi_rates) {
    # Use beta distribution to generate varying rates - different shapes create different patterns
    # Higher alpha skews right (more viruses with lower ZI), higher beta skews left (more viruses with higher ZI)
    virus_zi_rates <- rbeta(n_viruses, shape1 = zi_alpha, shape2 = zi_beta)
    
    # Scale the rates to be centered around the specified structural_zeros parameter
    # This ensures the average matches the input parameter while allowing variability
    scaling_factor <- structural_zeros / mean(virus_zi_rates)
    virus_zi_rates <- virus_zi_rates * scaling_factor
    
    # Bound rates to reasonable values
    virus_zi_rates <- pmin(0.95, pmax(0.1, virus_zi_rates))
    
    # Save for later output
    virus_specific_zi_rates <- virus_zi_rates
  }
  
  for (i in 1:n_viruses) {
    # Set base probability based on either fixed or variable rates
    if (variable_zi_rates) {
      # Use pre-generated virus-specific rate
      base_prob <- virus_zi_rates[i]
    } else {
      # Use the traditional abundance-based approach
      base_prob <- structural_zeros * (1 - abundance_means[i] / max(abundance_means))
      base_prob <- min(0.95, max(0.3, base_prob))  # Keep within reasonable range
      virus_specific_zi_rates[i] <- base_prob
    }
    
    # Group A probability
    structural_zero_probs[i, 1] <- base_prob
    
    # Group B probability (may differ if zero_inflation_difference is TRUE)
    if (zero_inflation_difference && i %in% diff_idx) {
      # For differentially abundant taxa, reduce structural zeros in group B
      # This creates a situation where the taxa are more prevalent in group B
      structural_zero_probs[i, 2] <- max(0.1, base_prob - 0.3)
    } else {
      structural_zero_probs[i, 2] <- base_prob
    }
  }
  
  # Generate data with zero-inflation and effect size differences
  for (i in 1:n_viruses) {
    # Determine base abundance for this taxon
    lambda <- abundance_means[i]
    
    # Process Group A samples
    for (j in group_a_indices) {
      # Determine if this is a structural zero
      if (runif(1) < structural_zero_probs[i, 1]) {
        counts[i, j] <- 0
        structural_zero_matrix[i, j] <- TRUE
      } else {
        # Determine if this is a sampling zero
        if (runif(1) < sampling_zeros) {
          counts[i, j] <- 0
          sampling_zero_matrix[i, j] <- TRUE
        } else {
          # Generate non-zero count from negative binomial
          counts[i, j] <- rnbinom(1, mu = lambda, size = dispersion)
        }
      }
    }
    
    # Process Group B samples
    for (j in group_b_indices) {
      # Determine if this is a structural zero
      if (runif(1) < structural_zero_probs[i, 2]) {
        counts[i, j] <- 0
        structural_zero_matrix[i, j] <- TRUE
      } else {
        # Determine if this is a sampling zero
        if (runif(1) < sampling_zeros) {
          counts[i, j] <- 0
          sampling_zero_matrix[i, j] <- TRUE
        } else {
          # Generate non-zero count from negative binomial
          # Apply effect size for differentially abundant taxa
          if (i %in% diff_idx) {
            counts[i, j] <- rnbinom(1, mu = lambda * effect_size, size = dispersion)
          } else {
            counts[i, j] <- rnbinom(1, mu = lambda, size = dispersion)
          }
        }
      }
    }
  }
  
  # Add additional virome-specific characteristics
  
  # 1. Some extremely abundant phages in a few samples (common in virome data)
  n_abundant_phages <- rpois(1, lambda = 2) + 1  # 1-3 abundant phages
  if (n_abundant_phages > 0 && n_viruses >= 5) {
    abundant_idx <- sample(1:n_viruses, n_abundant_phages)
    abundant_samples <- sample(1:n_samples, n_abundant_phages)
    
    for (i in seq_along(abundant_idx)) {
      # Make extremely abundant
      counts[abundant_idx[i], abundant_samples[i]] <- rnbinom(1, mu = mean(abundance_means) * 50, size = 1)
      # Clear zero indicators for these positions
      structural_zero_matrix[abundant_idx[i], abundant_samples[i]] <- FALSE
      sampling_zero_matrix[abundant_idx[i], abundant_samples[i]] <- FALSE
    }
  }
  
  # Calculate final sparsity statistics
  total_zeros <- sum(counts == 0)
  structural_zeros_count <- sum(structural_zero_matrix)
  sampling_zeros_count <- sum(sampling_zero_matrix)
  observed_sparsity <- total_zeros / (n_viruses * n_samples)
  
  # Return simulated data with zero-inflation details
  list(
    counts = counts,
    metadata = metadata,
    diff_taxa = diff_idx,  # indices of differentially abundant taxa
    structural_zeros = structural_zero_matrix,
    sampling_zeros = sampling_zero_matrix,
    parameters = list(
      n_samples = n_samples,
      n_viruses = n_viruses,
      structural_zeros = structural_zeros,
      sampling_zeros = sampling_zeros,
      dispersion = dispersion,
      effect_size = effect_size,
      zero_inflation_difference = zero_inflation_difference,
      variable_zi_rates = variable_zi_rates,
      zi_alpha = zi_alpha,
      zi_beta = zi_beta
    ),
    sparsity_summary = list(
      observed_sparsity = observed_sparsity,
      structural_zeros_proportion = structural_zeros_count / (n_viruses * n_samples),
      sampling_zeros_proportion = sampling_zeros_count / (n_viruses * n_samples),
      nonzero_proportion = 1 - observed_sparsity
    ),
    virus_specific_zi = list(
      rates = virus_specific_zi_rates,
      mean = mean(virus_specific_zi_rates),
      median = median(virus_specific_zi_rates),
      min = min(virus_specific_zi_rates),
      max = max(virus_specific_zi_rates),
      sd = sd(virus_specific_zi_rates),
      structural_zero_probs = structural_zero_probs
    )
  )
}