#!/usr/bin/env Rscript

# This test script demonstrates the fixed stratified power report generation

library(viromePower)

# Calculate stratified power with settings for high power
stratified_power <- calc_stratified_power(
  strata_sizes = c(25, 30),         # Larger sample sizes
  effect_sizes = c(4.0, 4.5),       # Larger effect sizes
  strata_weights = c(0.4, 0.6),     # Importance weights for strata
  n_viruses = 15,                   # Fewer viruses (less multiple testing burden)
  clustering_factor = 0.02,         # Lower clustering factor
  sparsity = 0.5,                   # Less sparsity
  dispersion = 1.2,                 # Lower dispersion (less variability)
  stratification_vars = "geography",
  method = "mixed_effects",
  n_sim = 10                        # Reduced simulations for faster testing
)

# Print the overall power
print(paste("Overall power:", round(stratified_power$overall_power * 100, 1), "%"))

# Print the stratum-specific power
print("Power by stratum:")
for (s in 1:length(stratified_power$stratum_specific_power)) {
  power_s <- stratified_power$stratum_specific_power[[s]]
  print(paste("  Stratum", s, ":", round(power_s * 100, 1), "%"))
}

# Generate the report
report_path <- generate_stratified_power_report(
  stratified_power_results = stratified_power,
  output_file = "stratified_power_report_fixed.html",
  title = "Fixed Stratified Virome Study Power Analysis",
  include_code = TRUE  # Include reproducible code in the report
)

print(paste("Report generated at:", report_path))
print("Test completed successfully.")