# Load the viromePower package
source("/Users/scott/Handley Lab Dropbox/Scott Handley/restructure/code/R/viromePower/R/calc_stratified_power.R")

# Create a wrapper for the README example
test_readme_example <- function() {
  # Calculate power for a stratified sampling design
  stratified_power <- calc_stratified_power(
    strata_sizes = c(10, 15),       # Samples per group in each stratum
    effect_sizes = c(1.8, 2.2),     # Different effect sizes by stratum
    strata_weights = c(0.4, 0.6),   # Importance weights for strata
    n_viruses = 50,                 # Number of viral taxa
    clustering_factor = 0.1,        # Intra-class correlation
    stratification_vars = "geography", # Variable used for stratification
    method = "mixed_effects"        # Statistical method for analysis
  )
  
  # View overall power
  print(paste("Overall power:", round(stratified_power$overall_power * 100, 1), "%"))
  
  # View stratum-specific power
  print("Power by stratum:")
  for (s in 1:length(stratified_power$stratum_specific_power)) {
    power_s <- stratified_power$stratum_specific_power[[s]]
    print(paste("  Stratum", s, ":", round(power_s * 100, 1), "%"))
  }
  
  # View effective sample size (accounting for design effect)
  print(paste("Effective sample size:", round(stratified_power$effective_sample_size, 1), 
              "(actual:", sum(stratified_power$parameters$strata_sizes) * 2, "samples)"))
  
  return(stratified_power)
}

# Run the example
result <- test_readme_example()