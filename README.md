# viromePower

## Overview
An R package for statistical power analysis of virome studies. This package helps researchers design robust studies by estimating required sample sizes and statistical power for detecting differences in viral abundances between groups, accounting for the unique characteristics of virome data.

## Features
- Simulate realistic virome data with appropriate sparsity and dispersion
- Calculate statistical power for various sample sizes and effect sizes
- Estimate required sample size for desired statistical power
- Generate interactive HTML reports with publication-quality visualizations
- Account for multiple testing correction and compositional nature of virome data
- Support for Bayesian and frequentist analysis approaches
- Stratified sampling designs for complex study structures
- Variable zero-inflation modeling for taxon-specific structural zero rates

## Robust Statistical Analysis
viromePower incorporates specially designed features to handle unique challenges in virome data:

- **Extreme sparsity handling**: Accounts for high proportion of zeros in virome datasets
- **Variable zero-inflation modeling**: Models taxon-specific structural zero rates using Beta distributions
- **Sample size boundaries**: Prevents unrealistic sample size estimates for small effect sizes
- **Statistical test safeguards**: Automatically adjusts for low sample sizes and high variability
- **Error handling and recovery**: Gracefully handles edge cases without crashing
- **Multiple statistical methods**: Supports Wilcoxon, t-test, DESeq-like, and Bayesian approaches
- **Smoothed power curves**: Reduces noise in visualization for clearer interpretation
- **Base64-embedded visualizations**: Self-contained HTML reports that work without external files

## Installation

```r
# Install from GitHub
devtools::install_github("username/viromePower")
```

## Quick Start

```r
library(viromePower)

# 1. Simulate realistic virome data 
# (handles sparsity, overdispersion, and rare taxa)
sim_data <- simulate_virome_data(
  n_samples = 30,           # Total number of samples
  n_viruses = 200,          # Number of viral taxa
  sparsity = 0.8,           # Proportion of zeros
  effect_size = 1.5         # Fold change between groups
)

# 2. Estimate sample size with robust bounds
# (avoids unrealistic extrapolation)
sample_size_result <- estimate_sample_size(
  power = 0.8,              # Target power
  effect_size = 3.0,        # Large effect size for clearer visualization
  n_viruses = 50,           # Moderate number of viral taxa
  sparsity = 0.5,           # Moderate sparsity level
  method = "t.test"         # Parametric test for normally-distributed data
)
# Extract the sample size from the result
sample_size_value <- sample_size_result$sample_size
print(paste("Recommended sample size:", sample_size_value, "samples per group"))

# 3. Calculate power for a specific sample size
# (handles edge cases and small sample sizes)
power_result <- calc_virome_power(
  n_samples = 15,           # Samples per group
  effect_size = 3.0,        # Same effect size as sample size estimation
  n_viruses = 50,           # Same number of viral taxa
  method = "t.test",        # Same statistical test
  sparsity = 0.5,           # Same sparsity level
  n_sim = 100               # Number of simulations
)

# Print power and false discovery rate 
power_percent <- round(power_result$power * 100, 1)
print(paste("Power:", power_percent, "%"))

# Print the false discovery rate
fdr_percent <- round(power_result$fdr * 100, 1)
print(paste("False discovery rate:", fdr_percent, "%"))

# 4. Plot power curve with smoothing
# (ensures reliable visualization)
power_curve <- plot_power_curve(
  effect_size = 3.0,            # Same effect size as previous examples
  n_viruses = 50,               # Same number of viral taxa
  sparsity = 0.5,               # Same sparsity level
  method = "t.test",            # Same statistical test
  sample_sizes = c(5, 10, 15, 20, 25, 30)  # Range of sample sizes to evaluate
)

# Access the sample size for 80% power
sample_size_80 <- attr(power_curve, "sample_size_80")
print(paste("Sample size for 80% power:", round(sample_size_80, 1)))

# 5. Generate comprehensive HTML report
# (embedded visualizations and interactive elements)
report <- generate_power_report(
  power_results = power_result,
  sample_size_results = sample_size_result,
  power_curve = power_curve,
  output_file = "virome_power_report.html",
  title = "Virome Study Power Analysis Report - High Effect Size Example"
)
```

The package handles the unique challenges of virome data:
- High sparsity and rare taxa detection
- Extreme abundance distributions
- Multiple testing correction
- Small sample sizes
- Edge cases with robust error handling
- Complex study designs with stratified sampling

## Stratified Sampling Example

For studies with structured populations (e.g., different geographical regions or host characteristics), the stratified power analysis provides more accurate estimates:

```r
# Calculate power for a stratified sampling design
stratified_power <- calc_stratified_power(
  strata_sizes = c(15, 20),         # Samples per group in each stratum
  effect_sizes = c(3.0, 3.5),       # Effect sizes by stratum
  strata_weights = c(0.4, 0.6),     # Importance weights for strata
  n_viruses = 20,                   # Number of viral taxa
  clustering_factor = 0.05,         # Intra-class correlation
  sparsity = 0.6,                   # Proportion of zeros
  dispersion = 1.5,                 # Dispersion parameter
  stratification_vars = "geography", # Variable used for stratification
  method = "mixed_effects"          # Statistical method for analysis
)

# View overall power
print(paste("Overall power:", round(stratified_power$overall_power * 100, 1), "%"))
# Example output: "Overall power: 47.5 %"

# View stratum-specific power
print("Power by stratum:")
for (s in 1:length(stratified_power$stratum_specific_power)) {
  power_s <- stratified_power$stratum_specific_power[[s]]
  print(paste("  Stratum", s, ":", round(power_s * 100, 1), "%"))
}
# Example output:
# "Power by stratum:"
# "  Stratum 1 : 2.5 %"
# "  Stratum 2 : 6.0 %"

# View effective sample size (accounting for design effect)
print(paste("Effective sample size:", round(stratified_power$effective_sample_size, 1), 
            "(actual:", sum(stratified_power$parameters$strata_sizes) * 2, "samples)"))
# Example output: "Effective sample size: 19.2 (actual: 70 samples)"

# High-power stratified analysis example
stratified_power <- calc_stratified_power(
  strata_sizes = c(25, 30),         # Larger sample sizes
  effect_sizes = c(4.0, 4.5),       # Larger effect sizes
  strata_weights = c(0.4, 0.6),     # Importance weights for strata
  n_viruses = 15,                   # Fewer viruses (less multiple testing burden)
  clustering_factor = 0.02,         # Lower clustering factor
  sparsity = 0.5,                   # Less sparsity
  dispersion = 1.2,                 # Lower dispersion (less variability)
  stratification_vars = "geography",
  method = "mixed_effects"
)

# Example output:
# Overall Power: 97.0%
# Power by Stratum:
#   Stratum 1 (n=25, effect=4.0): 39.0%
#   Stratum 2 (n=30, effect=4.5): 58.0%

# Generate comprehensive HTML report with visualizations
# Note: For enhanced formatting, install optional package: install.packages("kableExtra")
report_path <- generate_stratified_power_report(
  stratified_power_results = stratified_power,
  output_file = "stratified_power_report.html",
  title = "Stratified Virome Study Power Analysis",
  include_code = TRUE  # Include reproducible code in the report
)
```

The generated report includes:
- Interactive visualizations comparing power across strata
- Design effect and effective sample size analysis
- False discovery rate assessment
- Customized recommendations based on the analysis results
- Reproducible code for all calculations and visualizations

![Stratified Power Report Example](https://github.com/username/viromePower/raw/main/man/figures/stratified_report_example.png)

## Viral Diversity Power Analysis

For studies focused on community-level metrics rather than individual taxa, viromePower offers comprehensive diversity power analysis:

```r
# Calculate power for alpha diversity (Shannon index)
alpha_power <- calc_viral_diversity_power(
  n_samples = 20,             # 20 samples per group
  effect_size = 1.5,          # Moderate-large effect size (Cohen's d)
  n_viruses = 150,            # Number of viral taxa
  diversity_measure = "shannon", # Shannon diversity index
  sparsity = 0.75,            # Typical virome sparsity level
  dispersion = 1.8,           # Moderate dispersion parameter
  n_sim = 50                  # Number of simulation iterations
)

# Print power estimate for alpha diversity
print(paste("Alpha diversity power:", round(alpha_power$power * 100, 1), "%"))
# Example output: "Alpha diversity power: 82.0%"

# Calculate power for beta diversity (Bray-Curtis)
beta_power <- calc_viral_diversity_power(
  n_samples = 18,             # 18 samples per group
  effect_size = 0.25,         # Moderate effect size for beta diversity (R-squared scale)
  n_viruses = 200,            # Number of viral taxa
  diversity_measure = "bray", # Bray-Curtis dissimilarity
  sparsity = 0.75,            # Typical virome sparsity level
  dispersion = 2.0,           # Realistic viral dispersion
  n_sim = 50                  # Number of simulation iterations
)

# Print power estimate for beta diversity
print(paste("Beta diversity power:", round(beta_power$power * 100, 1), "%"))
# Example output: "Beta diversity power: 74.0%"

# Generate plot showing how power changes with sample size for alpha diversity
shannon_curve <- plot_diversity_power_curve(
  param_range = seq(10, 30, by = 5),  # Sample sizes to test
  param_type = "n_samples",           # Varying sample size
  effect_size = 1.5,                  # Fixed effect size (same as above)
  n_viruses = 150,                    # Number of viral taxa
  diversity_measure = "shannon",      # Shannon diversity
  sparsity = 0.75,                    # Typical virome sparsity
  dispersion = 1.8,                   # Consistent with above example
  n_sim = 30                          # Simulations per point
)

# Generate plot showing how power changes with effect size for beta diversity
bray_curve <- plot_diversity_power_curve(
  param_range = seq(0.1, 0.4, by = 0.05),    # Effect sizes to test (realistic range)
  param_type = "effect_size",                # Varying effect size
  n_samples = 18,                            # Fixed sample size (same as above)
  n_viruses = 200,                           # Number of viral taxa
  diversity_measure = "bray",                # Bray-Curtis dissimilarity
  sparsity = 0.75,                           # Typical virome sparsity
  dispersion = 2.0,                          # Consistent with above example
  n_sim = 30                                 # Simulations per point
)

# Generate comprehensive alpha diversity power report
alpha_report <- generate_diversity_power_report(
  n_samples = 20, 
  effect_size = 1.5,
  n_viruses = 150,
  diversity_measure = "shannon",
  sparsity = 0.75,
  dispersion = 1.8,
  output_file = "shannon_diversity_power_report.html"
)

# Generate comprehensive beta diversity power report
beta_report <- generate_diversity_power_report(
  n_samples = 18,
  effect_size = 0.25,
  n_viruses = 200,
  diversity_measure = "bray",
  sparsity = 0.75,
  dispersion = 2.0,
  output_file = "beta_diversity_power_report.html"
)
```

Viral diversity analysis offers several benefits:
- Captures community-level differences that may be missed by taxon-by-taxon analysis
- Reduces statistical burden (single test vs. multiple tests)
- Provides holistic view of virome structure changes
- Addresses compositional data challenges through dissimilarity measures
- Works well even with previously unknown or poorly characterized viruses

The diversity power reports include:
- Interactive visualizations of alpha or beta diversity patterns
- Power curves for both sample size and effect size
- Sample size recommendations for achieving target power
- Simulated data visualizations (boxplots for alpha, NMDS/PCoA for beta)
- Detailed interpretations and practical recommendations

![Diversity Power Report Example](https://github.com/username/viromePower/raw/main/man/figures/diversity_report_example.png)

## Bayesian Power Analysis

For researchers preferring Bayesian approaches, viromePower offers Bayesian power analysis that handles sparse data through prior incorporation and yields posterior probabilities instead of p-values:

```r
# Calculate Bayesian power for a virome study with stable settings
bayesian_power <- calc_bayesian_power(
  n_samples = 30,           # 30 samples per group
  effect_size = 4.0,        # 4-fold difference between groups
  n_viruses = 30,           # 30 viral taxa
  sparsity = 0.4,           # 40% zeros in the data  
  prior_strength = 5.0,     # Strong prior for stability
  n_sim = 10,               # 10 simulation iterations for example
  suppress_warnings = TRUE  # Suppress numerical calculation warnings
)

# Print Bayesian power estimate
print(paste("Bayesian power:", round(bayesian_power$power * 100, 1), "%"))
# Example output: "Bayesian power: 75.5%"

# Print expected discoveries
print(paste("Expected discoveries:", round(bayesian_power$expected_discoveries, 1)))
# Example output: "Expected discoveries: 18.2"

# Print false discovery proportion
print(paste("False discovery proportion:", round(bayesian_power$false_discovery_proportion * 100, 1), "%"))
# Example output: "False discovery proportion: 12.8%"

# Generate Bayesian power report
report_path <- generate_bayesian_power_report(
  bayesian_power_results = bayesian_power,
  output_file = "bayesian_power_report.html",
  title = "Bayesian Power Analysis for Virome Study",
  include_code = TRUE  # Include reproducible code in the report
)
```

Bayesian analysis offers several advantages for virome studies:
- Better handling of sparse data through prior information
- More intuitive credible intervals instead of confidence intervals
- Direct probability statements about effect existence
- No need for multiple testing correction
- Graceful performance with small sample sizes when using informative priors

The Bayesian report includes:
- Posterior probability distributions for differential abundance
- Bayes factor interpretations for evidence strength
- Comparison of frequentist and Bayesian approaches
- Precision-recall metrics for detection performance
- Customized recommendations based on the Bayesian analysis results

![Bayesian Power Report Example](https://github.com/username/viromePower/raw/main/man/figures/bayesian_report_example.png)

## Zero-Inflated Bayesian Power Analysis

Virome data is characterized by extreme sparsity that can't be adequately modeled by standard approaches. viromePower implements a zero-inflated Bayesian model that distinguishes between structural zeros (true absence) and sampling zeros (detection failures):

```r
# Generate zero-inflated virome data that realistically models both
# true absence (structural zeros) and detection failures (sampling zeros)
zinb_data <- simulate_zero_inflated_virome(
  n_samples = 30,             # 30 samples per group
  n_viruses = 200,            # 200 viral taxa
  structural_zeros = 0.7,     # 70% structural zeros (true absence)
  sampling_zeros = 0.2,       # 20% sampling zeros (detection failures)
  dispersion = 1.5,           # Dispersion parameter for negative binomial
  effect_size = 2.5,          # 2.5-fold difference between groups
  zero_inflation_difference = TRUE  # Allow prevalence differences between groups
)

# Calculate power using the zero-inflated Bayesian model
zinb_power <- calc_zinb_bayesian_power(
  n_samples = 25,             # 25 samples per group
  effect_size = 2.5,          # 2.5-fold difference between groups
  n_viruses = 150,            # 150 viral taxa
  structural_zeros = 0.7,     # 70% structural zeros (true absence)
  sampling_zeros = 0.2,       # 20% sampling zeros (detection failures)
  dispersion = 1.5,           # Dispersion parameter for NB distribution
  zero_inflation_difference = TRUE,  # Allow prevalence differences between groups
  prior_strength = 2.0,       # Moderate prior strength
  n_sim = 10                  # 10 simulation iterations for quick example
)

# Print zero-inflated Bayesian power estimate
print(paste("ZINB Bayesian power:", round(zinb_power$power * 100, 1), "%"))
# Example output: "ZINB Bayesian power: 68.5%"

# Compare with standard Bayesian approach (combined sparsity model)
standard_power <- calc_bayesian_power(
  n_samples = 25,             # Same sample size
  effect_size = 2.5,          # Same effect size
  n_viruses = 150,            # Same number of taxa
  sparsity = 0.8,             # Combined sparsity (structural + sampling)
  dispersion = 1.5,           # Same dispersion
  prior_strength = 2.0,       # Same prior strength
  n_sim = 10                  # Same simulation count
)

# Print comparison
print(paste("Standard Bayesian power:", round(standard_power$power * 100, 1), "%"))
print(paste("ZINB Bayesian power:", round(zinb_power$power * 100, 1), "%"))

# Examine zero-inflation statistics
print(paste("Structural zeros:", round(zinb_power$zero_inflation_summary$structural_zeros_proportion * 100, 1), "%"))
print(paste("Sampling zeros:", round(zinb_power$zero_inflation_summary$sampling_zeros_proportion * 100, 1), "%"))
print(paste("Total sparsity:", round(zinb_power$zero_inflation_summary$observed_sparsity * 100, 1), "%"))

# Generate comprehensive report for zero-inflated model analysis
# with specialized visualizations and diagnostics
zinb_report <- generate_zero_inflated_power_report(
  zinb_power_results = zinb_power,
  output_file = "zero_inflated_power_report.html",
  title = "Zero-Inflated Virome Power Analysis",
  target_power = 0.8  # Target 80% power for sample size recommendation
)

# Report includes:
# - Zero-inflation pattern visualization (structural vs sampling zeros)
# - Model comparison between standard and zero-inflated approaches
# - Power curves accounting for both zero types
# - Sample size recommendations that prevent underpowered designs
# - Discovery rate analysis with expected true/false positives
# - Decision tools for appropriate model selection
```

### Advanced Zero-Inflation Modeling with Variable Zero Rates Across Taxa

viromePower now supports variable zero-inflation rates across viral taxa, reflecting the biological reality that different viruses have different prevalence patterns:

```r
# Generate data with variable zero-inflation rates using Beta distribution
# Beta(2,5) creates right-skewed distribution (more viruses with lower ZI rates)
sim_variable <- simulate_zero_inflated_virome(
  n_samples = 40,              # 40 samples per group
  n_viruses = 150,             # 150 viral taxa
  structural_zeros = 0.7,      # Mean structural zero rate
  sampling_zeros = 0.2,        # Sampling zero rate
  variable_zi_rates = TRUE,    # Enable variable zero-inflation
  zi_alpha = 2,                # Beta distribution shape parameter
  zi_beta = 5,                 # Beta distribution shape parameter
  effect_size = 2.0            # 2-fold difference between groups
)

# Visualize the distribution of zero-inflation rates across viral taxa
if (require(ggplot2)) {
  # Plot the distribution of structural zero rates
  zi_dist_plot <- plot_virus_specific_zi_rates(sim_variable)
  print(zi_dist_plot)
  
  # Examine relationship between abundance and zero-inflation
  zi_abundance_plot <- plot_zi_by_abundance(sim_variable)
  print(zi_abundance_plot)
  
  # Compare zero-inflation patterns between groups
  zi_group_plot <- compare_group_zi_patterns(sim_variable, plot_type = "scatter")
  print(zi_group_plot)
}

# Calculate power using variable zero-inflation model
variable_zi_power <- calc_zinb_bayesian_power(
  n_samples = 30,              # 30 samples per group
  effect_size = 2.0,           # 2-fold difference
  n_viruses = 150,             # 150 viral taxa
  structural_zeros = 0.7,      # Mean structural zero rate
  sampling_zeros = 0.2,        # Sampling zero rate
  variable_zi_rates = TRUE,    # Enable variable zero-inflation
  zi_alpha = 2,                # Beta distribution shape parameter
  zi_beta = 5,                 # Beta distribution shape parameter
  n_sim = 10                   # 10 simulations for example
)

# Compare with uniform zero-inflation model
uniform_zi_power <- calc_zinb_bayesian_power(
  n_samples = 30,
  effect_size = 2.0,
  n_viruses = 150,
  structural_zeros = 0.7,
  sampling_zeros = 0.2,
  variable_zi_rates = FALSE,   # Use uniform zero-inflation
  n_sim = 10
)

# Print power comparison
print(paste("Variable ZI power:", round(variable_zi_power$power * 100, 1), "%"))
print(paste("Uniform ZI power:", round(uniform_zi_power$power * 100, 1), "%"))

# Generate report with variable zero-inflation diagnostics
var_zi_report <- generate_zero_inflated_power_report(
  zinb_power_results = variable_zi_power,
  output_file = "variable_zi_power_report.html",
  title = "Variable Zero-Inflation Power Analysis"
)
```

## Advanced Features for Variable Zero-Inflation Analysis

viromePower includes several advanced tools for working with variable zero-inflation in real virome studies:

```r
# FEATURE 1: ESTIMATE ZERO-INFLATION PARAMETERS FROM REAL DATA
# ---------------------------------------------------------
# Let's say you have a real virome count matrix (rows=viruses, cols=samples)
# We can estimate the best Beta distribution parameters to model the zero-inflation pattern

zi_params <- estimate_zi_parameters_from_data(
  count_matrix = my_virome_counts,         # Your virome count matrix
  presence_threshold = 5,                  # Minimum count to consider present
  method = "mle",                          # Maximum likelihood estimation
  plot_fit = TRUE                          # Create visualization of the fit
)

# Estimated parameters can then be used in simulations and power analysis
print(paste("Estimated alpha:", round(zi_params$alpha, 2)))
print(paste("Estimated beta:", round(zi_params$beta, 2)))
print(paste("Mean zero-inflation rate:", round(zi_params$mean_zi, 2)))

# FEATURE 2: COMPARE POWER ACROSS DIFFERENT ZERO-INFLATION PATTERNS
# ----------------------------------------------------------------
# Evaluate how different zero-inflation patterns affect statistical power

power_comparison <- plot_power_by_zi_pattern(
  sample_sizes = c(20, 30, 40, 50),        # Sample sizes to evaluate
  effect_sizes = c(1.5, 2.0, 2.5),         # Effect sizes to evaluate
  zi_patterns = list(                       # Define ZI patterns to compare
    right_skewed = c(alpha = 2, beta = 5),  # More taxa with lower ZI 
    left_skewed = c(alpha = 5, beta = 2),   # More taxa with higher ZI
    uniform = c(alpha = 1, beta = 1),       # Flat distribution
    fixed = c(alpha = NA, beta = NA)        # Traditional fixed rate
  ),
  fixed_zi_rate = 0.7,                     # Mean ZI rate for all patterns
  n_sim = 5                                # Use higher value (20+) for actual analysis
)

# Show comparison plot
print(power_comparison$plots$combined)

# FEATURE 3: COMPARE FIXED VS VARIABLE ZERO-INFLATION MODELS
# --------------------------------------------------------
# Determine which approach better fits your data

model_comparison <- compare_zi_models(
  count_matrix = my_virome_counts,         # Your virome count matrix
  variable_zi_method = "beta",             # Method for modeling variable ZI
  n_bootstrap = 100,                       # Bootstrap samples for confidence
  visualize = TRUE                         # Create comparison plots
)

# See which model fits better
print(paste("Better model:", model_comparison$best_model))
print(paste("AIC difference:", round(model_comparison$aic_values$difference, 2)))
print(paste("P-value:", format(model_comparison$lrt_p_value, digits = 4)))

# Visualization showing model comparison
print(model_comparison$plots$model_comparison)

# FEATURE 4: OPTIMIZE SAMPLING STRATEGY
# -----------------------------------
# Find optimal sample size and sequencing depth for a given ZI pattern

sampling_strategy <- optimize_sampling_strategy(
  target_power = 0.85,                     # Target statistical power
  effect_size = 2.0,                       # Expected effect size
  zi_alpha = zi_params$alpha,              # Use estimated alpha
  zi_beta = zi_params$beta,                # Use estimated beta
  total_sequencing_effort = 100,           # Total sequencing capacity (million reads)
  test_variable_zi = TRUE                  # Compare fixed vs variable ZI approaches
)

# Print recommendation
cat(sampling_strategy$recommendation)

# Plot showing power by sample size and depth
print(sampling_strategy$plots$contour_variable)

# FEATURE 5: EXPORT SIMULATED DATA FOR EXTERNAL ANALYSIS
# ----------------------------------------------------
# Export data to common formats for use with external tools

sim_data <- simulate_zero_inflated_virome(
  n_samples = 40,
  n_viruses = 200,
  variable_zi_rates = TRUE,
  zi_alpha = zi_params$alpha,
  zi_beta = zi_params$beta
)

export_results <- export_zinb_simulation(
  sim_data = sim_data,
  formats = c("csv", "biom", "phyloseq"),  # Formats to export
  output_dir = "virome_simulation",        # Directory for output files
  base_filename = "variable_zi_data",      # Base name for files
  include_zi_parameters = TRUE,            # Include ZI info in metadata
  normalize = TRUE,                        # Export normalized counts
  normalization_method = "CPM"             # Counts per million normalization
)

# Summary of exported files
print(paste("Exported", length(export_results$exported_files), "files"))
print(paste("Formats:", paste(export_results$summary$formats_exported, collapse = ", ")))
```

Zero-inflated Bayesian analysis provides several critical advantages for virome studies:
- Distinguishes between true absence and failed detection
- Accounts for differential viral prevalence between groups
- Models taxon-specific zero-inflation rates for better biological realism
- More accurate power estimates in extremely sparse datasets
- Identifies detection limits vs. biological effects
- Improves sensitivity for rare viral taxa discovery
- Provides detection-adjusted effect size estimates

The variable zero-inflation model offers these additional benefits:
- Each viral taxon has its own structural zero rate
- Reflects the biological reality that some viruses are genuinely rarer
- Captures abundance-prevalence relationships
- More accurate modeling of high taxonomic diversity
- Better estimation of true viral community structure

Using zero-inflated models is particularly important when:
1. Virome data has >80% zeros
2. Sequencing depth varies significantly between samples
3. Viral communities have many low-abundance members
4. Groups may differ in both abundance AND prevalence
5. Different viral taxa have distinct prevalence patterns

## Interactive Zero-Inflation Explorer Dashboard

viromePower includes an interactive Shiny dashboard for exploring zero-inflation patterns without writing code:

```r
# Launch the interactive Zero-Inflation Explorer dashboard
run_zi_explorer()
```

This dashboard provides easy-to-use interfaces for:
- Simulating and visualizing different zero-inflation patterns
- Estimating parameters from real virome data
- Comparing statistical power across different zero-inflation scenarios
- Determining optimal study designs for various zero-inflation patterns
- Statistically comparing fixed vs. variable zero-inflation models

![Zero-Inflation Explorer Dashboard](https://github.com/username/viromePower/raw/main/man/figures/zi_explorer_dashboard.png)

The dashboard is particularly useful for:
- Teaching concepts of zero-inflation in virome data
- Exploratory analysis of real datasets
- Quick visualization of different modeling approaches
- Interactive power analysis with real-time visualization
- Communicating results to collaborators and stakeholders

## Citation

If you use viromePower in your research, please cite:

[Citation information to be added]

## License

This project is licensed under the MIT License - see the LICENSE file for details.