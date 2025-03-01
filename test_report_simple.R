# Simplified test script for report generation
# Source the report generation function only
source("R/generate_power_report.R")

# Generate mock data for the report
power_result <- list(
  power = 0.78,
  fdr = 0.15,
  avg_detected = 12,
  parameters = list(
    n_samples = 10,
    effect_size = 1.5,
    n_viruses = 100,
    alpha = 0.05,
    sparsity = 0.8,
    method = "wilcoxon"
  )
)

sample_size_result <- list(
  sample_size = 15,
  achieved_power = 0.82,
  parameters = list(
    power = 0.8,
    effect_size = 1.5,
    n_viruses = 100,
    alpha = 0.05,
    sparsity = 0.8,
    method = "wilcoxon"
  )
)

# Generate HTML report
report_file <- "virome_power_report_simple.html"
cat("\nGenerating HTML report:", report_file, "\n")
report <- generate_power_report(
  power_results = power_result,
  sample_size_results = sample_size_result,
  power_curve = NULL,  # Skip the power curve
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