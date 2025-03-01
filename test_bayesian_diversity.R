#!/usr/bin/env Rscript

# This test script demonstrates the new Bayesian alpha diversity power analysis functionality

library(viromePower)

# Calculate Bayesian power for a virome study with alpha diversity analysis
# This example uses settings that should result in high power
power_result <- calc_bayesian_power(
  n_samples = 20,         # 20 samples per group
  effect_size = 3.0,      # 3-fold difference between groups
  n_viruses = 50,         # 50 viral taxa
  sparsity = 0.6,         # 60% zeros in the data
  prior_strength = 2.5,   # Moderate prior strength
  n_sim = 10,             # 10 simulation iterations for quick testing
  diversity_analysis = TRUE,  # Enable diversity analysis
  diversity_metrics = c("shannon", "simpson", "richness", "chao1", "fisher")  # Metrics to analyze
)

# Print the Bayesian power estimate for differential abundance
print(paste("Bayesian differential abundance power:", round(power_result$power * 100, 1), "%"))

# Print diversity results if available
if (!is.null(power_result$diversity_results)) {
  print("\nBayesian Alpha Diversity Results:")
  
  # Loop through each metric
  for (metric in names(power_result$diversity_results$summary)) {
    results <- power_result$diversity_results$summary[[metric]]
    
    print(paste("\n", toupper(metric), "Diversity:"))
    print(paste("  Power:", round(results$power * 100, 1), "%"))
    print(paste("  Bayes Factor:", round(results$bayes_factor, 2)))
    print(paste("  Effect Probability:", round(results$effect_probability, 3)))
    
    # Print posterior summaries
    print("  Posterior Summary:")
    print(paste("    Group A Mean:", round(results$posterior_summary$group_a$mean, 3)))
    print(paste("    Group B Mean:", round(results$posterior_summary$group_b$mean, 3)))
    print(paste("    Mean Difference:", round(results$posterior_summary$diff$mean, 3), 
                "(", round(results$posterior_summary$diff$lower, 3), "to", 
                round(results$posterior_summary$diff$upper, 3), ")"))
  }
}

print("\nTest completed successfully.")