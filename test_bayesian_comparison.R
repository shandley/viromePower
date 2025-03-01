#!/usr/bin/env Rscript

# This test script compares traditional and Bayesian power analysis approaches

library(viromePower)

# Shared parameters for both analyses
n_samples <- 15
effect_size <- 3.0
n_viruses <- 40
sparsity <- 0.7
n_sim <- 10  # Small number for testing purposes

# Run traditional power analysis
traditional_power <- calc_virome_power(
  n_samples = n_samples,
  effect_size = effect_size,
  n_viruses = n_viruses,
  sparsity = sparsity,
  method = "t.test",
  n_sim = n_sim
)

# Run Bayesian power analysis
bayesian_power <- calc_bayesian_power(
  n_samples = n_samples,
  effect_size = effect_size,
  n_viruses = n_viruses,
  sparsity = sparsity,
  prior_strength = 2.0,  # Moderate prior strength
  n_sim = n_sim
)

# Print comparison of results
print("POWER ANALYSIS COMPARISON")
print("========================")
print("Shared parameters:")
print(paste("  Samples per group:", n_samples))
print(paste("  Effect size:", effect_size))
print(paste("  Number of viruses:", n_viruses))
print(paste("  Sparsity:", sparsity))
print("")

print("Traditional Power Analysis:")
print(paste("  Power:", round(traditional_power$power * 100, 1), "%"))
print(paste("  False Discovery Rate:", round(traditional_power$fdr * 100, 1), "%"))
print(paste("  Expected Discoveries:", round(traditional_power$avg_detected, 1)))
print("")

print("Bayesian Power Analysis:")
print(paste("  Power:", round(bayesian_power$power * 100, 1), "%"))
print(paste("  False Discovery Proportion:", round(bayesian_power$false_discovery_proportion * 100, 1), "%"))
print(paste("  Expected Discoveries:", round(bayesian_power$expected_discoveries, 1)))
print(paste("  Average Posterior Probability:", round(mean(bayesian_power$posterior_prob_summary$mean), 2)))
print("")

# Generate both reports
traditional_report <- generate_power_report(
  power_results = traditional_power,
  output_file = "traditional_power_report.html",
  title = "Traditional Power Analysis for Virome Study"
)

bayesian_report <- generate_bayesian_power_report(
  bayesian_power_results = bayesian_power,
  output_file = "bayesian_power_report.html",
  title = "Bayesian Power Analysis for Virome Study"
)

print(paste("Traditional report generated at:", traditional_report))
print(paste("Bayesian report generated at:", bayesian_report))
print("Test completed successfully.")