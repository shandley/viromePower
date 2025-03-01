# Test script for HTML report generation
# This script will test the HTML report generation without requiring the full package installation

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

# Generate test data
power_result <- calc_virome_power(
  n_samples = 15,
  effect_size = 1.5,
  n_viruses = 100
)

sample_size <- estimate_sample_size(
  power = 0.8,
  effect_size = 1.5,
  n_viruses = 100
)

# Generate HTML report
report_file <- "virome_power_report.html"
cat("\nGenerating HTML report:", report_file, "\n")
report <- generate_power_report(
  power_results = power_result,
  sample_size_results = sample_size,
  power_curve = NULL,  # Let the function generate the curve
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