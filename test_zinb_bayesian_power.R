#!/usr/bin/env Rscript

# This test script demonstrates the new zero-inflated Bayesian power analysis functionality

library(viromePower)

# Calculate zero-inflated Bayesian power for a virome study
# This example uses settings that are realistic for highly sparse virome data
zinb_power <- calc_zinb_bayesian_power(
  n_samples = 25,             # 25 samples per group
  effect_size = 2.5,          # 2.5-fold difference between groups
  n_viruses = 150,            # 150 viral taxa
  structural_zeros = 0.75,    # 75% structural zeros (true absence)
  sampling_zeros = 0.15,      # 15% sampling zeros (detection failures)
  dispersion = 1.5,           # Moderate dispersion
  zero_inflation_difference = TRUE,  # Allow differences in zero-inflation between groups
  prior_strength = 2.0,       # Moderate prior strength
  n_sim = 5                   # 5 simulation iterations for quick testing
)

# Print the zero-inflated Bayesian power estimate
print(paste("ZINB Bayesian power:", round(zinb_power$power * 100, 1), "%"))

# Print the expected number of discoveries
print(paste("Expected discoveries:", round(zinb_power$expected_discoveries, 1)))

# Print the false discovery proportion
print(paste("False discovery proportion:", round(zinb_power$false_discovery_proportion * 100, 1), "%"))

# Print zero-inflation summary
print("\nZero-inflation summary:")
print(paste("  Structural zeros proportion:", round(zinb_power$zero_inflation_summary$structural_zeros_proportion * 100, 1), "%"))
print(paste("  Sampling zeros proportion:", round(zinb_power$zero_inflation_summary$sampling_zeros_proportion * 100, 1), "%"))
print(paste("  Total observed sparsity:", round(zinb_power$zero_inflation_summary$observed_sparsity * 100, 1), "%"))
print(paste("  Detection rate:", round(zinb_power$zero_inflation_summary$detection_rate * 100, 1), "%"))
print(paste("  Differential detection between groups:", round(zinb_power$zero_inflation_summary$differential_detection * 100, 1), "%"))

# Print Bayes factor summary
print("\nBayes factor summary:")
print(paste("  Mean:", round(zinb_power$bayes_factor_summary$mean, 2)))
print(paste("  Median:", round(zinb_power$bayes_factor_summary$median, 2)))
print(paste("  25th percentile:", round(zinb_power$bayes_factor_summary$q25, 2)))
print(paste("  75th percentile:", round(zinb_power$bayes_factor_summary$q75, 2)))

# Compare with standard Bayesian power analysis
print("\nComparing with standard Bayesian power analysis (ignoring zero-inflation):")
standard_power <- calc_bayesian_power(
  n_samples = 25,             # Same parameters as above
  effect_size = 2.5,
  n_viruses = 150,
  sparsity = 0.8,             # Approximating the combined sparsity
  dispersion = 1.5,
  prior_strength = 2.0,
  n_sim = 5                   # Same number of simulations
)

print(paste("Standard Bayesian power:", round(standard_power$power * 100, 1), "%"))
print(paste("ZINB Bayesian power:", round(zinb_power$power * 100, 1), "%"))
print(paste("Power difference:", round((zinb_power$power - standard_power$power) * 100, 1), "percentage points"))

print("\nTest completed successfully.")