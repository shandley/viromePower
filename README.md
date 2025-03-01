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
  sample_size_results = sample_size,
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

## Citation

If you use viromePower in your research, please cite:

[Citation information to be added]

## License

This project is licensed under the MIT License - see the LICENSE file for details.
