# Load the necessary source files
source("/Users/scott/Handley Lab Dropbox/Scott Handley/restructure/code/R/viromePower/R/calc_stratified_power.R")

# Try different parameter combinations to find one with better power
test_higher_power_example <- function() {
  # Use larger effect sizes and smaller number of viruses to increase power
  stratified_power <- calc_stratified_power(
    strata_sizes = c(15, 20),         # Larger sample sizes
    effect_sizes = c(3.0, 3.5),       # Larger effect sizes
    strata_weights = c(0.4, 0.6),     # Importance weights
    n_viruses = 20,                   # Fewer viruses (less multiple testing burden)
    clustering_factor = 0.05,         # Lower clustering (less variance)
    sparsity = 0.6,                   # Less sparsity (more non-zero data)
    dispersion = 1.5,                 # Lower dispersion (less variance)
    stratification_vars = "geography",
    method = "mixed_effects",
    n_sim = 100                       # More simulations for better estimate
  )
  
  # View overall power
  print(paste("Overall power:", round(stratified_power$overall_power * 100, 1), "%"))
  
  # View stratum-specific power
  print("Power by stratum:")
  for (s in 1:length(stratified_power$stratum_specific_power)) {
    power_s <- stratified_power$stratum_specific_power[[s]]
    print(paste("  Stratum", s, ":", round(power_s * 100, 1), "%"))
  }
  
  # View effective sample size
  print(paste("Effective sample size:", round(stratified_power$effective_sample_size, 1), 
              "(actual:", sum(stratified_power$parameters$strata_sizes) * 2, "samples)"))
  
  return(stratified_power)
}

# Run the example
result <- test_higher_power_example()