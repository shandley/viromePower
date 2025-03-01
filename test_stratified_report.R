# Test script for generating a stratified power report

# Source the necessary functions
source("R/calc_stratified_power.R")
source("R/generate_stratified_power_report.R")

# Calculate stratified power with parameters that provide meaningful power
stratified_power <- calc_stratified_power(
  strata_sizes = c(15, 20),         # Samples per group in each stratum
  effect_sizes = c(3.0, 3.5),       # Effect sizes by stratum
  strata_weights = c(0.4, 0.6),     # Importance weights for strata
  n_viruses = 20,                   # Number of viral taxa
  clustering_factor = 0.05,         # Intra-class correlation
  sparsity = 0.6,                   # Proportion of zeros
  dispersion = 1.5,                 # Dispersion parameter
  stratification_vars = "geography", # Variable used for stratification
  method = "mixed_effects",         # Statistical method for analysis
  n_sim = 50                        # Reduce for quicker testing
)

# Print basic results
cat("Overall power:", round(stratified_power$overall_power * 100, 1), "%\n")
cat("Power by stratum:\n")
for (s in 1:length(stratified_power$stratum_specific_power)) {
  power_s <- stratified_power$stratum_specific_power[[s]]
  cat("  Stratum", s, ":", round(power_s * 100, 1), "%\n")
}
cat("Effective sample size:", round(stratified_power$effective_sample_size, 1), 
    "(actual:", sum(stratified_power$parameters$strata_sizes) * 2, "samples)\n")

# Generate a report
cat("\nGenerating stratified power analysis report...\n")
output_file <- "stratified_power_report.html"

# Check if we have the required packages
required_packages <- c("rmarkdown", "ggplot2", "dplyr", "tidyr", "knitr", "kableExtra")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Missing packages required for report generation:", paste(missing_packages, collapse = ", "), "\n")
  cat("Please install them with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))\n")
} else {
  # Generate the report
  tryCatch({
    report_path <- generate_stratified_power_report(
      stratified_power_results = stratified_power,
      output_file = output_file,
      title = "Stratified Virome Power Analysis Report",
      include_code = TRUE
    )
    cat("Report generated successfully at:", report_path, "\n")
  }, error = function(e) {
    cat("Error generating report:", e$message, "\n")
  })
}