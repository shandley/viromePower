#!/usr/bin/env Rscript

# This test script demonstrates the new Bayesian power analysis functionality

library(viromePower)

# Calculate Bayesian power for a virome study
# This example uses settings that should result in high power
power_result <- calc_bayesian_power(
  n_samples = 20,         # 20 samples per group
  effect_size = 3.0,      # 3-fold difference between groups
  n_viruses = 50,         # 50 viral taxa
  sparsity = 0.6,         # 60% zeros in the data
  prior_strength = 2.5,   # Moderate prior strength
  n_sim = 10              # 10 simulation iterations for quick testing
)

# Print the Bayesian power estimate
print(paste("Bayesian power:", round(power_result$power * 100, 1), "%"))

# Print the expected number of discoveries
print(paste("Expected discoveries:", round(power_result$expected_discoveries, 1)))

# Print the false discovery proportion
print(paste("False discovery proportion:", round(power_result$false_discovery_proportion * 100, 1), "%"))

# Print Bayes factor summary
print("Bayes factor summary:")
print(paste("  Mean:", round(power_result$bayes_factor_summary$mean, 2)))
print(paste("  Median:", round(power_result$bayes_factor_summary$median, 2)))

# Print posterior probability summary
print("Posterior probability summary:")
print(paste("  Mean:", round(power_result$posterior_prob_summary$mean, 2)))
print(paste("  Median:", round(power_result$posterior_prob_summary$median, 2)))

# Generate a Bayesian power report
report_path <- generate_bayesian_power_report(
  bayesian_power_results = power_result,
  output_file = "bayesian_power_report.html",
  title = "Bayesian Power Analysis for Virome Study",
  include_code = TRUE  # Include reproducible code in the report
)

print(paste("Report generated at:", report_path))
print("Test completed successfully.")