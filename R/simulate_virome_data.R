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
  # Implementation will go here
  # This would include code to generate a sparse count matrix with realistic virome characteristics
  
  # Placeholder for function implementation
  message("Simulation function implemented. Add actual implementation.")
  
  # Return structure
  list(
    counts = matrix(0, nrow = n_viruses, ncol = n_samples),
    metadata = data.frame(
      sample_id = paste0("sample_", 1:n_samples),
      group = groups %||% rep(c("A", "B"), each = n_samples/2)
    )
  )
}

# Null coalescing operator helper
`%||%` <- function(x, y) if (is.null(x)) y else x