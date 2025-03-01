# Test script for real implementations
# Set seed for reproducibility
set.seed(123)

# Source all R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in r_files) {
  cat("Sourcing", file, "\n")
  source(file)
}

# Test 1: Simulate data with small parameters
cat("\nTest 1: Simulate small dataset\n")
sim_small <- simulate_virome_data(n_samples = 10, n_viruses = 20)
cat("- Data dimensions:", dim(sim_small$counts)[1], "viruses x", dim(sim_small$counts)[2], "samples\n")
cat("- Number of true diff. taxa:", length(sim_small$diff_taxa), "\n")
cat("- Sparsity (% zeros):", mean(sim_small$counts == 0) * 100, "%\n")

# Test 2: Calculate power with small parameters
cat("\nTest 2: Calculate power with small parameters\n")
power_small <- calc_virome_power(
  n_samples = 5,
  effect_size = 2.0,
  n_viruses = 20,
  n_sim = 5
)
cat("- Power:", power_small$power, "\n")
cat("- FDR:", power_small$fdr, "\n")

# Test 3: Test all statistical methods
cat("\nTest 3: Compare statistical methods\n")
methods <- c("wilcoxon", "t.test", "deseq")
for (m in methods) {
  power_result <- calc_virome_power(
    n_samples = 5,
    effect_size = 2.0,
    n_viruses = 20,
    method = m,
    n_sim = 3
  )
  cat("- Method:", m, "| Power:", power_result$power, "\n")
}

# Test 4: Estimate sample size
cat("\nTest 4: Estimate sample size\n")
sample_size <- estimate_sample_size(
  power = 0.8,
  effect_size = 1.5,
  n_viruses = 20
)
cat("- Estimated sample size:", sample_size$sample_size, "\n")
cat("- Achieved power:", sample_size$achieved_power, "\n")