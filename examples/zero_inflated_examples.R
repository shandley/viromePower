#!/usr/bin/env Rscript

#####################################################################
# Example Script: Using Zero-Inflated Model Visualization Functions #
#####################################################################

# Load required libraries
library(viromePower)
library(ggplot2)
# gridExtra is recommended for combined plots but optional
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  message("The 'gridExtra' package is recommended for optimal visualization. Install with: install.packages('gridExtra')")
}

# Set seed for reproducibility
set.seed(42)

#------------------------------------------------------
# 1. Simulate zero-inflated virome data
#------------------------------------------------------
# This function simulates data with both structural zeros (true absence)
# and sampling zeros (detection failures)
sim_data <- simulate_zero_inflated_virome(
  n_samples = 30,             # 30 samples (15 per group)
  n_viruses = 100,            # 100 viral taxa
  structural_zeros = 0.7,     # 70% structural zeros (true absence)
  sampling_zeros = 0.2,       # 20% sampling zeros among non-structural zeros
  dispersion = 1.5,           # Dispersion parameter for negative binomial
  effect_size = 2.0,          # 2-fold difference between groups for differential taxa
  zero_inflation_difference = TRUE  # Allow zero-inflation to differ between groups
)

# Display basic information about the simulated data
cat("Simulated virome dataset:\n")
cat("  Number of samples:", sim_data$parameters$n_samples, "\n")
cat("  Number of viruses:", sim_data$parameters$n_viruses, "\n")
cat("  Observed sparsity:", round(sim_data$sparsity_summary$observed_sparsity * 100, 1), "%\n")
cat("  Number of differentially abundant taxa:", length(sim_data$diff_taxa), "\n\n")

#------------------------------------------------------
# 2. Create mock model fits
#------------------------------------------------------
# a) Define standard negative binomial model parameters
standard_model_params <- list(
  mu = 10,    # Mean count
  size = 0.5  # Dispersion parameter
)

# b) Define zero-inflated negative binomial model parameters
zinb_model_params <- list(
  mu = 10,       # Mean count
  size = 0.5,    # Dispersion parameter
  zi_prob = 0.7  # Zero-inflation probability
)

# c) Create a mock zero-inflation parameter matrix (simulates a fitted model)
n_samples <- nrow(sim_data$counts)
n_viruses <- ncol(sim_data$counts)

# Mock ZINB model fit with per-sample and per-taxon parameters
zinb_model_fit <- list(
  mu = matrix(rgamma(n_samples * n_viruses, 5, 0.5), nrow = n_samples),
  size = 0.5,
  zi_prob = matrix(rbeta(n_samples * n_viruses, 2, 5), nrow = n_samples),
  pi = matrix(rbeta(n_samples * n_viruses, 2, 5), nrow = n_samples),  # Zero-inflation probabilities
  taxa_names = colnames(sim_data$counts)
)

#------------------------------------------------------
# 3. Calculate mock power results for different sample sizes
#------------------------------------------------------
# Define sample sizes to evaluate
sample_sizes <- seq(10, 50, by = 10)

# Create mock power results for standard and zero-inflated models
# In practice, you would use calc_bayesian_power and calc_zinb_bayesian_power 
# but they may not exist yet in your implementation
std_power <- c(0.25, 0.45, 0.60, 0.72, 0.82)
zinb_power <- c(0.40, 0.65, 0.78, 0.85, 0.92)

# Standard model results
std_power_results <- list(
  n_samples = sample_sizes,
  power = std_power,
  effect_size = rep(2.0, length(sample_sizes))
)

# Zero-inflated model results
zinb_power_results <- list(
  n_samples = sample_sizes,
  power = zinb_power,
  effect_size = rep(2.0, length(sample_sizes))
)

#------------------------------------------------------
# 4. Create visualizations using the four functions
#------------------------------------------------------

# a) Plot observed vs expected zeros
# This plot compares observed zero counts with expected zero counts
# under standard and zero-inflated negative binomial models
cat("Creating plot 1: Observed vs Expected Zeros...\n")
zeros_plot <- plot_observed_vs_expected_zeros(
  obs_data = t(sim_data$counts),  # Transpose to get samples as rows
  standard_model_params = standard_model_params,
  zinb_model_params = zinb_model_params,
  group_factor = sim_data$metadata$group,  # Group factor for coloring points
  title = "Observed vs Expected Zeros in Virome Data"
)

# b) Plot zero-inflation distribution
# This plot shows the distribution of zero-inflation parameters
# across samples or taxa
cat("Creating plot 2: Zero-Inflation Distribution...\n")
# By sample
zi_sample_plot <- plot_zero_inflation_distribution(
  zinb_model_fit = zinb_model_fit,
  by_taxa = FALSE,  # View distribution by sample
  group_factor = sim_data$metadata$group,
  highlight_threshold = 0.5  # Highlight this threshold
)

# By taxa
zi_taxa_plot <- plot_zero_inflation_distribution(
  zinb_model_fit = zinb_model_fit,
  by_taxa = TRUE,  # View distribution by taxa
  highlight_threshold = 0.7  # Highlight this threshold
)

# c) Compare power curves
# This plot compares statistical power between standard and
# zero-inflated negative binomial models
cat("Creating plot 3: Power Curve Comparison...\n")
power_curve_plot <- compare_power_curves(
  standard_power_results = std_power_results,
  zinb_power_results = zinb_power_results,
  x_variable = "n_samples",  # Compare across sample sizes
  highlight_thresh = 0.8     # Highlight 80% power threshold
)

# d) Plot ZINB diagnostics
# This creates diagnostic plots for assessing the fit of
# zero-inflated negative binomial models
cat("Creating plot 4: ZINB Diagnostics...\n")
diagnostics <- plot_zinb_diagnostics(
  obs_data = t(sim_data$counts),  # Transpose to get samples as rows
  zinb_model_fit = zinb_model_fit,
  n_taxa_to_plot = 4,  # Number of taxa to show in detailed plots
  seed = 123           # For reproducible taxa selection
)

#------------------------------------------------------
# 5. Display the visualizations
#------------------------------------------------------
cat("\nPlot objects created successfully!\n")
cat("In an interactive R session, you can view the plots using:\n")
cat("  - zeros_plot            # Observed vs expected zeros plot\n")
cat("  - zi_sample_plot        # Zero-inflation distribution by sample\n")
cat("  - zi_taxa_plot          # Zero-inflation distribution by taxa\n")
cat("  - power_curve_plot      # Power curve comparison\n")
cat("  - diagnostics$zero_proportion_plot  # Zero proportion diagnostic plot\n")
cat("  - diagnostics$qq_plot              # QQ plot for model fit\n")
cat("  - diagnostics$taxa_fit_plots       # Detailed taxa fit plots\n\n")

cat("You can save plots to files using ggsave(), for example:\n")
cat("  ggsave('zeros_plot.png', zeros_plot, width=10, height=6)\n\n")

#------------------------------------------------------
# 6. Print summary information
#------------------------------------------------------
cat("Zero-Inflation Summary Statistics:\n")
cat("  Structural zeros:", round(sim_data$sparsity_summary$structural_zeros_proportion * 100, 1), "%\n")
cat("  Sampling zeros:", round(sim_data$sparsity_summary$sampling_zeros_proportion * 100, 1), "%\n")
cat("  Total sparsity:", round(sim_data$sparsity_summary$observed_sparsity * 100, 1), "%\n\n")

cat("Power Analysis Comparison:\n")
cat("  Standard model power (n=50):", round(std_power_results$power[length(std_power_results$power)] * 100, 1), "%\n")
cat("  Zero-inflated model power (n=50):", round(zinb_power_results$power[length(zinb_power_results$power)] * 100, 1), "%\n")