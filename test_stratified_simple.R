# Test script for generating a simplified stratified power analysis

# Source the necessary function
source("R/calc_stratified_power.R")

# Calculate stratified power with parameters that provide meaningful power
stratified_power <- calc_stratified_power(
  strata_sizes = c(15, 20),         # Samples per group in each stratum
  effect_sizes = c(3.0, 3.5),       # Effect sizes by stratum
  strata_weights = c(0.4, 0.6),     # Importance weights for strata
  n_viruses = 20,                   # Number of viral taxa
  clustering_factor = 0.05,         # Intra-class correlation
  sparsity = 0.6,                   # Proportion of zeros
  dispersion = 1.5,                 # Dispersion parameter
  stratification_vars = "geography", # Variable used for stratification
  method = "mixed_effects",         # Statistical method for analysis
  n_sim = 50                        # Reduce for quicker testing
)

# Print results in a nicely formatted way
cat("\n=== Stratified Power Analysis Report ===\n\n")
cat(sprintf("Overall Power: %.1f%%\n", stratified_power$overall_power * 100))

# Print stratum-specific power
cat("\nPower by Stratum:\n")
for (s in 1:length(stratified_power$stratum_specific_power)) {
  power_s <- stratified_power$stratum_specific_power[[s]]
  cat(sprintf("  Stratum %d (n=%d, effect=%.1f): %.1f%%\n", 
              s, stratified_power$parameters$strata_sizes[s], 
              stratified_power$parameters$effect_sizes[s], 
              power_s * 100))
}

# Print effective sample size
cat("\nSample Size Information:\n")
cat(sprintf("  Total samples: %d\n", sum(stratified_power$parameters$strata_sizes) * 2))
cat(sprintf("  Effective sample size: %.1f\n", stratified_power$effective_sample_size))
cat(sprintf("  Design effect: %.2f\n", stratified_power$design_effect_actual))

# Print false discovery rate information
cat("\nFalse Discovery Rate Analysis:\n")
cat(sprintf("  False discovery rate: %.1f%%\n", stratified_power$fdr * 100))
cat(sprintf("  Average detections per study: %.1f\n", stratified_power$avg_detected))
cat(sprintf("  True positives: %.2f\n", mean(stratified_power$sim_summary$true_positives)))
cat(sprintf("  False positives: %.2f\n", mean(stratified_power$sim_summary$false_positives)))

# Print recommendations
cat("\nRecommendations:\n")
if (stratified_power$overall_power < 0.8) {
  cat("  - Power is below the conventional 80% threshold\n")
  cat("  - Consider increasing sample sizes or focusing on strata with higher effect sizes\n")
} else {
  cat("  - Power exceeds the conventional 80% threshold\n")
  cat("  - The design has adequate power to detect the specified effect sizes\n")
}

if (stratified_power$design_effect_actual > 1.5) {
  cat("  - The stratification substantially reduces the effective sample size\n")
  cat("  - Consider reducing the clustering or allocating more samples to key strata\n")
} else if (stratified_power$design_effect_actual > 1.2) {
  cat("  - The stratification moderately reduces the effective sample size\n")
  cat("  - The benefits of stratification should be weighed against this reduction\n")
} else {
  cat("  - The stratification has minimal impact on the effective sample size\n")
  cat("  - This indicates an efficient stratified design\n")
}

cat("\n========================================\n")