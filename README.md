# viromePower

## Overview
An R package for statistical power analysis of virome studies. This package helps researchers design robust studies by estimating required sample sizes and statistical power for detecting differences in viral abundances between groups, accounting for the unique characteristics of virome data.

## Features
- Simulate realistic virome data with appropriate sparsity and dispersion
- Calculate statistical power for various sample sizes and effect sizes
- Estimate required sample size for desired statistical power
- Generate interactive HTML reports with publication-quality visualizations
- Account for multiple testing correction and compositional nature of virome data

## Installation

```r
# Install from GitHub
devtools::install_github("username/viromePower")
```

## Quick Start

```r
library(viromePower)

# Estimate sample size needed for 80% power
sample_size <- estimate_sample_size(
  power = 0.8,
  effect_size = 1.5,
  n_viruses = 100
)

# Calculate power for a specific sample size
power_result <- calc_virome_power(
  n_samples = 15,
  effect_size = 1.5,
  n_viruses = 100
)

# Plot power curve
power_curve <- plot_power_curve(
  effect_size = 1.5,
  n_viruses = 100,
  sample_sizes = seq(5, 30, by = 5)
)

# Generate comprehensive HTML report
report <- generate_power_report(
  power_results = power_result,
  sample_size_results = sample_size,
  power_curve = power_curve
)
```

## Citation

If you use viromePower in your research, please cite:

[Citation information to be added]

## License

This project is licensed under the MIT License - see the LICENSE file for details.
