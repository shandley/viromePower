# Load the viromePower package
source("/Users/scott/Handley Lab Dropbox/Scott Handley/restructure/code/R/viromePower/R/calc_stratified_power.R")
source("/Users/scott/Handley Lab Dropbox/Scott Handley/restructure/code/R/viromePower/R/simulate_virome_data.R")
source("/Users/scott/Handley Lab Dropbox/Scott Handley/restructure/code/R/viromePower/R/calc_virome_power.R")

# Test the function
result <- calc_stratified_power(
  strata_sizes = c(10, 15),       # Samples per group in each stratum
  effect_sizes = c(1.8, 2.2),     # Different effect sizes by stratum
  strata_weights = c(0.4, 0.6),   # Importance weights for strata
  n_viruses = 50,                 # Number of viral taxa
  clustering_factor = 0.1,        # Intra-class correlation
  stratification_vars = "geography", # Variable used for stratification
  method = "mixed_effects"        # Statistical method for analysis
)

# Print the result
print(result)