# Test script for the fixed stratified report generation function
# This will check if the function works but skips the report generation if packages are missing

# Source the necessary functions
source("R/calc_stratified_power.R")
source("R/generate_stratified_power_report.R")

# Calculate stratified power
cat("Calculating stratified power...\n")
stratified_power <- calc_stratified_power(
  strata_sizes = c(25, 30),     # Samples per group in each stratum
  effect_sizes = c(4.0, 4.5),   # Effect sizes by stratum
  strata_weights = c(0.4, 0.6), # Importance weights for strata
  n_viruses = 15,              # Number of viral taxa
  clustering_factor = 0.02,    # Intra-class correlation
  sparsity = 0.5,              # Proportion of zeros
  dispersion = 1.2,            # Lower dispersion
  stratification_vars = "geography",
  method = "mixed_effects",
  n_sim = 10                   # Small number for quick testing
)

# Check results
cat("\nPower analysis results:\n")
cat(sprintf("Overall power: %.1f%%\n", stratified_power$overall_power * 100))

# Check for required packages for visualization
required_packages <- c("rmarkdown", "ggplot2", "dplyr", "tidyr", "knitr", "kableExtra")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("\nMissing packages required for report generation:", paste(missing_packages, collapse = ", "), "\n")
  cat("If you want to generate HTML reports, install them with:\n")
  cat(sprintf('install.packages(c("%s"))\n', paste(missing_packages, collapse = '", "')))
  
  # Skip the report generation
  cat("\nSkipping HTML report generation due to missing packages.\n")
} else {
  # Try to generate the report
  tryCatch({
    cat("\nGenerating HTML report...\n")
    output_file <- "stratified_power_report.html"
    
    report_path <- generate_stratified_power_report(
      stratified_power_results = stratified_power,
      output_file = output_file,
      title = "Stratified Virome Study Power Analysis",
      include_code = TRUE
    )
    
    cat("\nReport successfully generated at:", report_path, "\n")
  }, error = function(e) {
    cat("\nError generating report:", e$message, "\n")
  })
}