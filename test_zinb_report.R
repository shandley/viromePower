#!/usr/bin/env Rscript

#############################################################
# Test Script: Zero-Inflated Virome Power Analysis Workflow #
#############################################################

# Load the viromePower package
library(viromePower)

# Set seed for reproducibility
set.seed(123)

# Step 1: Simulate zero-inflated virome data
# Using parameters typical for a virome study
cat("Simulating zero-inflated virome data...\n")
sim_data <- simulate_zero_inflated_virome(
  n_samples = 20,             # 20 samples per group (40 total)
  n_viruses = 200,            # 200 viral taxa
  structural_zeros = 0.75,    # 75% structural zeros (true absence)
  sampling_zeros = 0.15,      # 15% sampling zeros (detection failures)
  dispersion = 1.5,           # Dispersion parameter for negative binomial
  effect_size = 2.0,          # 2-fold change for differential taxa
  zero_inflation_difference = TRUE  # Allow zero-inflation to differ between groups
)

# Print basic information about the simulated data
cat("\nSimulated virome dataset:\n")
cat("  Number of samples:", ncol(sim_data$counts), "\n")
cat("  Number of viruses:", nrow(sim_data$counts), "\n")
cat("  Observed sparsity:", round(sim_data$sparsity_summary$observed_sparsity * 100, 1), "%\n")
cat("  Number of differentially abundant taxa:", length(sim_data$diff_taxa), "\n\n")

# Step 2: Calculate zero-inflated Bayesian power
# Creating a mock power result for testing (usually you would use calc_zinb_bayesian_power)
# This is a simplified version for testing without requiring the full computation
zinb_power <- list(
  power = 0.72,  # 72% power
  parameters = list(
    n_samples = 20,
    n_viruses = 200,
    effect_size = 2.0,
    structural_zeros = 0.75,
    sampling_zeros = 0.15,
    dispersion = 1.5,
    posterior_prob_threshold = 0.95
  ),
  zero_inflation_summary = list(
    observed_sparsity = sim_data$sparsity_summary$observed_sparsity,
    detection_rate = 0.30,  # 30% detection rate
    differential_detection = 0.15  # 15% differential detection between groups
  ),
  false_discovery_proportion = 0.05,  # 5% FDR
  sim_summary = list(
    true_positives = 25,  # 25 true discoveries
    false_positives = 2   # 2 false discoveries
  ),
  expected_discoveries = 27  # Total discoveries
)

# Print power results
cat("\nZero-inflated Bayesian power analysis results:\n")
cat("  Statistical power:", round(zinb_power$power * 100, 1), "%\n")
cat("  Expected discoveries:", round(zinb_power$expected_discoveries), "\n")
cat("  False discovery proportion:", round(zinb_power$false_discovery_proportion * 100, 1), "%\n\n")

# Step 3: Generate HTML power analysis report
cat("Generating zero-inflated power analysis report...\n")
report_path <- generate_zero_inflated_power_report(
  zinb_power_results = zinb_power,
  output_file = "zinb_power_report.html",
  title = "Zero-Inflated Virome Power Analysis Report",
  include_code = FALSE,
  target_power = 0.8  # Target 80% power
)

cat("\nWorkflow completed successfully!\n")
cat("Power analysis report saved to:", report_path, "\n")