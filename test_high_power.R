# Test script for high-power stratified analysis example

# Source the necessary function
source("R/calc_stratified_power.R")

# Calculate stratified power with parameters optimized for high power
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
  n_sim = 50                        
)

# Print results in a nicely formatted way
cat("\n=== High-Power Stratified Analysis Report ===\n\n")
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

cat("\n==============================================\n")