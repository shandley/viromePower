# Test script for generating a stratified power report with the fixed function

# Source the necessary functions
source("R/calc_stratified_power.R")
source("R/generate_stratified_power_report.R")

# Set up error handling
tryCatch({
  # Calculate stratified power with parameters that provide high power
  cat("Calculating stratified power analysis...\n")
  stratified_power <- calc_stratified_power(
    strata_sizes = c(25, 30),         # Larger sample sizes
    effect_sizes = c(4.0, 4.5),       # Larger effect sizes
    strata_weights = c(0.4, 0.6),     # Importance weights for strata
    n_viruses = 15,                   # Fewer viruses (less multiple testing burden)
    clustering_factor = 0.02,         # Lower clustering factor
    sparsity = 0.5,                   # Less sparsity
    dispersion = 1.2,                 # Lower dispersion (less variability)
    stratification_vars = "geography", 
    method = "mixed_effects",         
    n_sim = 10                        # Small number for quick testing
  )
  
  # Print basic results
  cat("\nPower analysis complete!\n")
  cat(sprintf("Overall Power: %.1f%%\n", stratified_power$overall_power * 100))
  
  # Check for required packages
  required_packages <- c("rmarkdown", "ggplot2", "dplyr", "tidyr", "knitr", "kableExtra")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    cat("\nMissing packages required for report generation:", paste(missing_packages, collapse = ", "), "\n")
    cat("Please install them with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))\n")
    cat("Report generation will be skipped.\n")
  } else {
    # Generate the report
    cat("\nGenerating stratified power report...\n")
    output_file <- "stratified_power_report.html"
    
    report_path <- generate_stratified_power_report(
      stratified_power_results = stratified_power,
      output_file = output_file,
      title = "Stratified Virome Study Power Analysis",
      include_code = TRUE
    )
    
    cat("\nReport generated successfully at:", report_path, "\n")
  }
}, error = function(e) {
  cat("\nError in report generation:", e$message, "\n")
  cat("Stack trace:\n")
  print(traceback())
})