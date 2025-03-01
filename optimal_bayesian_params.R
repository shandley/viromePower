#!/usr/bin/env Rscript

# This script tests different parameter combinations to find settings
# that consistently produce meaningful power values for the README example

# Install the package locally
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::document()
devtools::install(build_vignettes = FALSE)

library(viromePower)

# Fixed seed for reproducibility
set.seed(123)

# Function to test a parameter combination
test_params <- function(n_samples, effect_size, n_viruses, sparsity, prior_strength, n_sim = 5) {
  cat("\nTesting parameters:\n")
  cat("  n_samples =", n_samples, "\n")
  cat("  effect_size =", effect_size, "\n")
  cat("  n_viruses =", n_viruses, "\n")
  cat("  sparsity =", sparsity, "\n")
  cat("  prior_strength =", prior_strength, "\n")
  
  result <- calc_bayesian_power(
    n_samples = n_samples,
    effect_size = effect_size,
    n_viruses = n_viruses,
    sparsity = sparsity,
    prior_strength = prior_strength,
    n_sim = n_sim,
    suppress_warnings = TRUE
  )
  
  cat("RESULTS:\n")
  cat("  Power:", round(result$power * 100, 1), "%\n")
  cat("  Expected discoveries:", round(result$expected_discoveries, 1), "\n")
  cat("  False discovery proportion:", round(result$false_discovery_proportion * 100, 1), "%\n")
  
  return(result)
}

# Test different combinations
results <- list()

# Test 1: Baseline from README
results[["baseline"]] <- test_params(
  n_samples = 20,
  effect_size = 3.0,
  n_viruses = 50,
  sparsity = 0.6,
  prior_strength = 2.5
)

# Test 2: Higher effect size
results[["high_effect"]] <- test_params(
  n_samples = 20,
  effect_size = 4.0,
  n_viruses = 50,
  sparsity = 0.6,
  prior_strength = 2.5
)

# Test 3: More samples
results[["more_samples"]] <- test_params(
  n_samples = 30, 
  effect_size = 3.0,
  n_viruses = 50,
  sparsity = 0.6,
  prior_strength = 2.5
)

# Test 4: Less sparsity
results[["less_sparse"]] <- test_params(
  n_samples = 20,
  effect_size = 3.0,
  n_viruses = 50,
  sparsity = 0.5,
  prior_strength = 2.5
)

# Test 5: Stronger prior
results[["stronger_prior"]] <- test_params(
  n_samples = 20,
  effect_size = 3.0,
  n_viruses = 50,
  sparsity = 0.6,
  prior_strength = 4.0
)

# Test 6: Optimal combination
results[["optimal"]] <- test_params(
  n_samples = 25,
  effect_size = 3.5,
  n_viruses = 40,
  sparsity = 0.5,
  prior_strength = 3.0
)

# Print summary of all results
cat("\n\n=== RESULTS SUMMARY ===\n")
for (name in names(results)) {
  cat(sprintf("%-15s: Power = %5.1f%%, Expected discoveries = %5.1f, FDP = %5.1f%%\n",
              name,
              results[[name]]$power * 100,
              results[[name]]$expected_discoveries,
              results[[name]]$false_discovery_proportion * 100))
}

cat("\nRecommended parameters for README example:\n")
cat("  n_samples = 25\n")
cat("  effect_size = 3.5\n")
cat("  n_viruses = 40\n")
cat("  sparsity = 0.5\n")
cat("  prior_strength = 3.0\n")
cat("  suppress_warnings = TRUE\n")