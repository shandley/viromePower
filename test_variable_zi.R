# Test script for variable zero-inflation rates in virome power analysis
library(devtools)
load_all() # Load the package for testing without installing

# Set seed for reproducibility
set.seed(12345)

# Generate data with uniform zero-inflation rates (baseline)
sim_uniform <- simulate_zero_inflated_virome(
  n_samples = 40,
  n_viruses = 150,
  structural_zeros = 0.7,
  sampling_zeros = 0.2,
  variable_zi_rates = FALSE
)

# Generate data with variable zero-inflation rates using beta distribution
# Beta(2,5) creates right-skewed distribution (more viruses with lower ZI rates)
sim_variable_right <- simulate_zero_inflated_virome(
  n_samples = 40,
  n_viruses = 150,
  structural_zeros = 0.7,
  sampling_zeros = 0.2,
  variable_zi_rates = TRUE,
  zi_alpha = 2,
  zi_beta = 5
)

# Beta(5,2) creates left-skewed distribution (more viruses with higher ZI rates)
sim_variable_left <- simulate_zero_inflated_virome(
  n_samples = 40,
  n_viruses = 150,
  structural_zeros = 0.7,
  sampling_zeros = 0.2,
  variable_zi_rates = TRUE,
  zi_alpha = 5,
  zi_beta = 2
)

# Beta(2,2) creates symmetric distribution centered at 0.5
sim_variable_symmetric <- simulate_zero_inflated_virome(
  n_samples = 40,
  n_viruses = 150,
  structural_zeros = 0.7,
  sampling_zeros = 0.2,
  variable_zi_rates = TRUE,
  zi_alpha = 2,
  zi_beta = 2
)

# Visualize zero-inflation rate distributions
if (require(ggplot2) && require(gridExtra)) {
  # Plot distributions for each simulation
  p1 <- plot_virus_specific_zi_rates(sim_uniform) + 
    ggtitle("Uniform Zero-Inflation")
  
  p2 <- plot_virus_specific_zi_rates(sim_variable_right) + 
    ggtitle("Right-Skewed Zero-Inflation (Beta(2,5))")
  
  p3 <- plot_virus_specific_zi_rates(sim_variable_left) + 
    ggtitle("Left-Skewed Zero-Inflation (Beta(5,2))")
  
  p4 <- plot_virus_specific_zi_rates(sim_variable_symmetric) + 
    ggtitle("Symmetric Zero-Inflation (Beta(2,2))")
  
  # Arrange in a grid
  grid.arrange(p1, p2, p3, p4, ncol = 2)
  
  # Plot relationship between abundance and zero-inflation rates
  p5 <- plot_zi_by_abundance(sim_variable_right)
  p6 <- plot_zi_by_abundance(sim_variable_left)
  
  grid.arrange(p5, p6, ncol = 2)
  
  # Compare group-specific patterns
  p7 <- compare_group_zi_patterns(sim_variable_right, plot_type = "boxplot")
  p8 <- compare_group_zi_patterns(sim_variable_right, plot_type = "scatter")
  
  grid.arrange(p7, p8, ncol = 2)
}

# Run power analysis with variable zero-inflation rates
power_uniform <- calc_zinb_bayesian_power(
  n_samples = 30,
  effect_size = 2.0,
  n_viruses = 100,
  structural_zeros = 0.7,
  sampling_zeros = 0.2,
  n_sim = 10 # Reduced for quicker testing
)

power_variable <- calc_zinb_bayesian_power(
  n_samples = 30,
  effect_size = 2.0,
  n_viruses = 100,
  structural_zeros = 0.7,
  sampling_zeros = 0.2,
  n_sim = 10, # Reduced for quicker testing
  variable_zi_rates = TRUE,
  zi_alpha = 2,
  zi_beta = 5
)

# Compare power results
cat("\nPower with uniform zero-inflation:", round(power_uniform$power * 100, 1), "%\n")
cat("Power with variable zero-inflation:", round(power_variable$power * 100, 1), "%\n")
cat("Difference:", round((power_variable$power - power_uniform$power) * 100, 1), "%\n")

# Generate report
if (interactive()) {
  report_variable <- generate_zero_inflated_power_report(
    zinb_power_results = power_variable,
    output_file = "variable_zi_power_report.html",
    title = "Variable Zero-Inflation Virome Power Analysis"
  )
  
  # Open the report in browser
  browseURL(report_variable)
}