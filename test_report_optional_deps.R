# Test script for report generation with optional dependencies
# This will show that the report can be generated even when kableExtra is not available

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

# Check if required packages are available
has_kableExtra <- requireNamespace("kableExtra", quietly = TRUE)
has_rmarkdown <- requireNamespace("rmarkdown", quietly = TRUE)
has_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)
has_knitr <- requireNamespace("knitr", quietly = TRUE)

cat("\nPackage availability:\n")
cat("kableExtra:", if(has_kableExtra) "Available" else "Not available", "\n")
cat("rmarkdown:", if(has_rmarkdown) "Available" else "Not available", "\n")
cat("ggplot2:", if(has_ggplot2) "Available" else "Not available", "\n")
cat("knitr:", if(has_knitr) "Available" else "Not available", "\n")

# Try to generate the report (will work with or without kableExtra)
if (has_rmarkdown && has_ggplot2 && has_knitr) {
  cat("\nGenerating stratified power report...\n")
  output_file <- "stratified_power_report_basic.html"
  
  tryCatch({
    report_path <- generate_stratified_power_report(
      stratified_power_results = stratified_power,
      output_file = output_file,
      title = "Stratified Virome Study Power Analysis (with optional dependencies)",
      include_code = TRUE
    )
    
    cat("\nReport successfully generated at:", report_path, "\n")
    cat("The report was generated", if(has_kableExtra) "WITH" else "WITHOUT", "kableExtra for enhanced formatting.\n")
    
  }, warning = function(w) {
    cat("Warning during report generation:", conditionMessage(w), "\n")
  }, error = function(e) {
    cat("\nError generating report:", conditionMessage(e), "\n")
  })
} else {
  cat("\nCannot generate report - missing essential packages (rmarkdown, ggplot2, or knitr).\n")
}