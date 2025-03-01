#' Simulate Virome Data
#'
#' Generates simulated virome data with characteristics typical of viral metagenomic studies,
#' including high sparsity, low abundance, and high variability.
#'
#' @param n_samples Number of samples to simulate
#' @param n_viruses Number of viral taxa to simulate
#' @param sparsity Proportion of zeros in the data (default: 0.8)
#' @param dispersion Dispersion parameter for negative binomial distribution (default: 2)
#' @param effect_size Effect size for group differences (default: 1.5)
#' @param groups Vector of group labels for samples (default: equal groups)
#'
#' @return A list containing simulated count matrix and sample metadata
#' @export
#'
#' @examples
#' sim_data <- simulate_virome_data(n_samples = 20, n_viruses = 100)
simulate_virome_data <- function(n_samples, n_viruses, sparsity = 0.8, 
                                dispersion = 2, effect_size = 1.5, 
                                groups = NULL) {
  # Create group labels if not provided
  if (is.null(groups)) {
    # Default: Half samples in group A, half in group B
    if (n_samples %% 2 != 0) {
      warning("Odd number of samples, adjusting for balanced groups")
      n_samples <- n_samples + 1
    }
    groups <- rep(c("A", "B"), each = n_samples/2)
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
  # Use exponential distribution for abundance means (creates power law-like pattern)
  abundance_means <- rexp(n_viruses, rate = 0.1)
  abundance_means <- abundance_means / sum(abundance_means) * 1e6  # Scale to reasonable counts
  
  # Make a subset of taxa differentially abundant between groups
  n_diff <- round(n_viruses * 0.1)  # 10% of taxa are differentially abundant
  diff_idx <- sample(n_viruses, n_diff)
  
  # Create matrix to store counts
  counts <- matrix(0, nrow = n_viruses, ncol = n_samples)
  colnames(counts) <- metadata$sample_id
  rownames(counts) <- paste0("virus_", 1:n_viruses)
  
  # Generate data with sparsity and effect size differences
  for (i in 1:n_viruses) {
    # Determine base abundance for this taxon
    lambda <- abundance_means[i]
    
    # Generate sparsity pattern
    # Higher abundance taxa generally have lower sparsity
    # Add safety check to prevent NaN or extreme values
    zero_prob <- min(0.95, sparsity * exp(-lambda/(mean(abundance_means) + 1e-10)))
    zero_prob <- max(0.05, zero_prob)  # Always keep some data points
    is_present <- runif(n_samples) > zero_prob
    
    # Generate counts for present observations using negative binomial
    # (common distribution for microbiome/virome count data)
    counts_nonzero <- numeric(n_samples)
    counts_nonzero[is_present] <- rnbinom(
      sum(is_present), 
      mu = lambda,
      size = dispersion  # dispersion parameter (lower = more overdispersion)
    )
    
    # Apply effect size to differentially abundant taxa
    if (i %in% diff_idx) {
      # Apply effect size to group B (scales abundance by effect_size)
      is_group_b <- metadata$group == "B"
      counts_nonzero[is_present & is_group_b] <- rnbinom(
        sum(is_present & is_group_b),
        mu = lambda * effect_size,
        size = dispersion
      )
    }
    
    counts[i, ] <- counts_nonzero
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
    }
  }
  
  # 2. Add streak of zeros (barcode bleeding effect, common in sequencing)
  streak_zeros <- sample(1:n_viruses, round(n_viruses * 0.05))  # 5% of taxa affected
  streak_samples <- sample(1:n_samples, round(n_samples * 0.2))  # in 20% of samples
  counts[streak_zeros, streak_samples] <- 0
  
  # Return simulated data
  list(
    counts = counts,
    metadata = metadata,
    diff_taxa = diff_idx,  # indices of differentially abundant taxa
    parameters = list(
      n_samples = n_samples,
      n_viruses = n_viruses,
      sparsity = sparsity,
      dispersion = dispersion,
      effect_size = effect_size
    )
  )
}

# Null coalescing operator helper
`%||%` <- function(x, y) if (is.null(x)) y else x