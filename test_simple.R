# Simple test script for viromePower package

# Source all R files directly instead of loading as a package
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in r_files) {
  cat("Sourcing", file, "\n")
  source(file)
}

# Test simulate_virome_data
cat("\n1. Testing simulate_virome_data...\n")
sim_data <- simulate_virome_data(n_samples = 20, n_viruses = 100)
cat("  - Function executed successfully\n")
cat("  - Data dimensions:", dim(sim_data$counts)[1], "viruses x", 
    dim(sim_data$counts)[2], "samples\n")
cat("  - Metadata rows:", nrow(sim_data$metadata), "\n")

# Test calc_virome_power
cat("\n2. Testing calc_virome_power...\n")
power_result <- calc_virome_power(
  n_samples = 10,
  effect_size = 1.5,
  n_viruses = 100
)
cat("  - Function executed successfully\n")
cat("  - Calculated power:", power_result$power, "\n")

# Test estimate_sample_size
cat("\n3. Testing estimate_sample_size...\n")
sample_size <- estimate_sample_size(
  power = 0.8,
  effect_size = 1.5,
  n_viruses = 100
)
cat("  - Function executed successfully\n")
cat("  - Estimated sample size:", sample_size$sample_size, "\n")

# Test generate_power_report
cat("\n4. Testing generate_power_report...\n")
report <- generate_power_report(
  power_results = power_result,
  sample_size_results = sample_size,
  power_curve = NULL,  # We'll skip the power curve for this test
  output_file = "test_report.html"
)

# Check if the report was created
if (file.exists("test_report.html")) {
  cat("  - Report file created: test_report.html\n")
} else {
  cat("  - Report file was not created\n")
}

cat("\nAll tests completed.\n")