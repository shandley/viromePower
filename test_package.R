# Script to test the viromePower package
# Run this script from the package directory with:
# Rscript test_package.R

# Load required packages
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
library(devtools)

# Load the package in development mode
load_all(".")
cat("Package loaded in development mode.\n")

# Set seed for reproducibility
set.seed(123)
cat("Testing functions...\n")

# Test simulate_virome_data
cat("\n1. Testing simulate_virome_data...\n")
tryCatch({
  sim_data <- simulate_virome_data(n_samples = 20, n_viruses = 100)
  cat("  - Function executed successfully\n")
  cat("  - Data dimensions:", dim(sim_data$counts)[1], "viruses x", 
      dim(sim_data$counts)[2], "samples\n")
  cat("  - Metadata rows:", nrow(sim_data$metadata), "\n")
}, error = function(e) {
  cat("  - ERROR:", conditionMessage(e), "\n")
})

# Test calc_virome_power
cat("\n2. Testing calc_virome_power...\n")
tryCatch({
  power_result <- calc_virome_power(
    n_samples = 10,
    effect_size = 1.5,
    n_viruses = 100
  )
  cat("  - Function executed successfully\n")
  cat("  - Calculated power:", power_result$power, "\n")
  cat("  - Parameters:", 
      "n_samples =", power_result$parameters$n_samples, 
      "effect_size =", power_result$parameters$effect_size,
      "n_viruses =", power_result$parameters$n_viruses, "\n")
}, error = function(e) {
  cat("  - ERROR:", conditionMessage(e), "\n")
})

# Test estimate_sample_size
cat("\n3. Testing estimate_sample_size...\n")
tryCatch({
  sample_size <- estimate_sample_size(
    power = 0.8,
    effect_size = 1.5,
    n_viruses = 100
  )
  cat("  - Function executed successfully\n")
  cat("  - Estimated sample size:", sample_size$sample_size, "\n")
  cat("  - Parameters:", 
      "power =", sample_size$parameters$power, 
      "effect_size =", sample_size$parameters$effect_size,
      "n_viruses =", sample_size$parameters$n_viruses, "\n")
}, error = function(e) {
  cat("  - ERROR:", conditionMessage(e), "\n")
})

# Test plot_power_curve
cat("\n4. Testing plot_power_curve...\n")
tryCatch({
  power_curve <- plot_power_curve(
    effect_size = 1.5,
    n_viruses = 100,
    sample_sizes = seq(5, 30, by = 5)
  )
  cat("  - Function executed successfully\n")
  cat("  - Plot object created (not displayed in console)\n")
  # Save plot to file if ggplot2 is available
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    ggplot2::ggsave("test_power_curve.png", power_curve, width = 8, height = 6)
    cat("  - Plot saved to test_power_curve.png\n")
  }
}, error = function(e) {
  cat("  - ERROR:", conditionMessage(e), "\n")
})

# Test generate_power_report
cat("\n5. Testing generate_power_report...\n")
tryCatch({
  report <- generate_power_report(
    power_results = power_result,
    sample_size_results = sample_size,
    power_curve = power_curve,
    output_file = "test_report.html"
  )
  cat("  - Function executed successfully\n")
  if (file.exists("test_report.html")) {
    cat("  - Report file created: test_report.html\n")
  } else {
    cat("  - Report file was not created\n")
  }
}, error = function(e) {
  cat("  - ERROR:", conditionMessage(e), "\n")
})

# Run tests
cat("\n6. Running automated tests...\n")
tryCatch({
  test_results <- test()
  cat("  - Tests completed\n")
}, error = function(e) {
  cat("  - ERROR:", conditionMessage(e), "\n")
})

cat("\nAll tests completed. Check the output above for any errors.\n")
cat("For a full package check, run: devtools::check()\n")