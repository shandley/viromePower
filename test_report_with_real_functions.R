# Test script for HTML report generation with real functions
# Set seed for reproducibility
set.seed(456)

# Source all R files directly
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in r_files) {
  cat("Sourcing", file, "\n")
  source(file)
}

# Check and install required packages for HTML report
needed_packages <- c("ggplot2", "base64enc")
for (pkg in needed_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("Installing package:", pkg, "\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Generate test data with small parameters for quick execution
power_result <- calc_virome_power(
  n_samples = 5, 
  effect_size = 2.0,
  n_viruses = 20,
  n_sim = 5  # Small number for quick test
)

# Generate a simple power curve
sample_sizes <- c(3, 5, 7, 10)
power_curve <- plot_power_curve(
  effect_size = 2.0,
  n_viruses = 20,
  sample_sizes = sample_sizes
)

# Estimate sample size
sample_size_result <- list(
  sample_size = 8,
  achieved_power = 0.82,
  parameters = list(
    power = 0.8,
    effect_size = 2.0,
    n_viruses = 20,
    alpha = 0.05,
    sparsity = 0.8,
    dispersion = 2,
    method = "wilcoxon"
  )
)

# Generate HTML report
report_file <- "virome_power_report_real.html"
cat("\nGenerating HTML report:", report_file, "\n")
report <- generate_power_report(
  power_results = power_result,
  sample_size_results = sample_size_result,
  power_curve = power_curve,
  output_file = report_file,
  title = "Virome Study Power Analysis Report"
)

# Check if the report was created
if (file.exists(report_file)) {
  cat("SUCCESS: HTML report created successfully!\n")
  cat("Open", report_file, "in your web browser to view the report.\n")
} else {
  cat("ERROR: Failed to create the HTML report.\n")
}

cat("\nTest completed.\n")