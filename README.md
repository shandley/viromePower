# viromePower

## Overview
An R package for statistical power analysis of virome studies. This package helps researchers design robust studies by estimating required sample sizes and statistical power for detecting differences in viral abundances between groups, accounting for the unique characteristics of virome data.

## Features
- Simulate realistic virome data with appropriate sparsity and dispersion
- Calculate statistical power for various sample sizes and effect sizes
- Estimate required sample size for desired statistical power
- Generate interactive HTML reports with publication-quality visualizations
- Account for multiple testing correction and compositional nature of virome data

## Robust Statistical Analysis
viromePower incorporates specially designed features to handle unique challenges in virome data:

- **Extreme sparsity handling**: Accounts for high proportion of zeros in virome datasets
- **Sample size boundaries**: Prevents unrealistic sample size estimates for small effect sizes
- **Statistical test safeguards**: Automatically adjusts for low sample sizes and high variability
- **Error handling and recovery**: Gracefully handles edge cases without crashing
- **Multiple statistical methods**: Supports Wilcoxon, t-test, and DESeq-like approaches
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
```

## Citation

If you use viromePower in your research, please cite:

[Citation information to be added]

## License

This project is licensed under the MIT License - see the LICENSE file for details.
