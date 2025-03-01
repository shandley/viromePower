# Final test to check for messages
# Capture output to check for "Add actual implementation" messages
output <- capture.output({
  # Source all R files
  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  for (file in r_files) {
    source(file)
  }
  
  # Run minimal versions of each function
  sim_data <- simulate_virome_data(n_samples = 6, n_viruses = 10)
  
  power_result <- calc_virome_power(
    n_samples = 3,
    effect_size = 2.0,
    n_viruses = 10,
    n_sim = 2
  )
  
  sample_size <- list(
    sample_size = 5,
    achieved_power = 0.7,
    parameters = list(
      power = 0.8,
      effect_size = 2.0,
      n_viruses = 10,
      alpha = 0.05,
      sparsity = 0.8,
      method = "wilcoxon"
    )
  )
  
  # Create a temporary file for the report
  report_file <- tempfile(fileext = ".html")
  report <- generate_power_report(
    power_results = power_result,
    sample_size_results = sample_size,
    power_curve = NULL,
    output_file = report_file
  )
})

# Check output for implementation messages
implementation_msgs <- grep("Add actual implementation", output)
if (length(implementation_msgs) > 0) {
  cat("Found implementation messages in output:\n")
  for (i in implementation_msgs) {
    cat("  ", output[i], "\n")
  }
} else {
  cat("No implementation messages found in output. All good!\n")
}