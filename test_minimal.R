# Minimal test script for real implementations
# Set seed for reproducibility
set.seed(123)

# Source all R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in r_files) {
  cat("Sourcing", file, "\n")
  source(file)
}

# Test 1: Simulate data with tiny parameters
cat("\nTest 1: Simulate tiny dataset\n")
sim_tiny <- simulate_virome_data(n_samples = 6, n_viruses = 10)
cat("- Data dimensions:", dim(sim_tiny$counts)[1], "viruses x", dim(sim_tiny$counts)[2], "samples\n")
cat("- Sparsity (% zeros):", mean(sim_tiny$counts == 0) * 100, "%\n")

# Test 2: Calculate power with tiny parameters and very few simulations
cat("\nTest 2: Calculate power with tiny parameters\n")
power_tiny <- calc_virome_power(
  n_samples = 3,
  effect_size = 2.0,
  n_viruses = 10,
  n_sim = 2  # Very few simulations for quick testing
)
cat("- Power:", power_tiny$power, "\n")