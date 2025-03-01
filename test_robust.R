# Test robustness improvements
# Set seed for reproducibility
set.seed(123)

# Source all R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in r_files) {
  cat("Sourcing", file, "\n")
  source(file)
}

# Test 1: Simulation with small parameters
cat("\nTest 1: Simulate dataset with robustness improvements\n")
sim_data <- simulate_virome_data(n_samples = 10, n_viruses = 20)
cat("- Data dimensions:", dim(sim_data$counts)[1], "viruses x", dim(sim_data$counts)[2], "samples\n")
cat("- Number of true diff. taxa:", length(sim_data$diff_taxa), "\n")
cat("- Sparsity (% zeros):", mean(sim_data$counts == 0) * 100, "%\n")

# Test 2: Calculate power with small parameters (should handle edge cases better)
cat("\nTest 2: Calculate power with improved robustness\n")
power_result <- calc_virome_power(
  n_samples = 5,
  effect_size = 2.0,
  n_viruses = 20,
  n_sim = 5
)
cat("- Power:", power_result$power, "\n")
cat("- FDR:", power_result$fdr, "\n")

# Test 3: Estimate sample size with bounded upper limit
cat("\nTest 3: Sample size estimation with reasonable bounds\n")
sample_size <- estimate_sample_size(
  power = 0.8,
  effect_size = 1.5,
  n_viruses = 20
)
cat("- Estimated sample size:", sample_size$sample_size, "\n")
cat("- Achieved power:", sample_size$achieved_power, "\n")
cat("- Min/max sample sizes tested:", min(sample_size$power_by_n$sample_size), 
    "/", max(sample_size$power_by_n$sample_size), "\n")

# Test 4: Power curve with smoothing
cat("\nTest 4: Power curve with smoothing\n")
power_curve <- plot_power_curve(
  effect_size = 1.5,
  n_viruses = 20,
  sample_sizes = c(3, 5, 7, 10, 15, 20)
)
cat("- Power curve sample sizes:", paste(attr(power_curve, "power_data")$sample_size, collapse=", "), "\n")
cat("- Power values:", paste(round(attr(power_curve, "power_data")$power, 2), collapse=", "), "\n")
cat("- Sample size for 80% power:", attr(power_curve, "sample_size_80"), "\n")