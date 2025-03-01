# Testing Guide for viromePower

This document provides a comprehensive guide to testing the viromePower package during development.

## 1. Setting Up the Testing Environment

First, make sure you have the necessary development packages:

```r
install.packages(c("devtools", "testthat", "roxygen2", "knitr", "rmarkdown"))
```

## 2. Loading the Package in Development Mode

When testing during development, load the package directly from the source:

```r
library(devtools)
setwd("/path/to/viromePower")  # Set to your package directory
load_all()  # Loads the package in development mode
```

## 3. Running Automated Tests

### Run All Tests

To run the entire test suite:

```r
devtools::test()
```

### Run a Specific Test File

To run only specific test files:

```r
devtools::test_file("tests/testthat/test-simulate_virome_data.R")
devtools::test_file("tests/testthat/test-calc_virome_power.R")
```

## 4. Testing Individual Functions

Test each main function to verify it works as expected:

### Simulate Virome Data

```r
# Basic simulation
sim_data <- simulate_virome_data(
  n_samples = 20, 
  n_viruses = 100
)

# Check structure
str(sim_data)

# Verify dimensions
dim(sim_data$counts)
nrow(sim_data$metadata)

# Try with custom parameters
sim_data2 <- simulate_virome_data(
  n_samples = 30,
  n_viruses = 50,
  sparsity = 0.9,
  dispersion = 3,
  effect_size = 2.0
)
```

### Calculate Power

```r
# Basic power calculation
power_result <- calc_virome_power(
  n_samples = 10,
  effect_size = 1.5,
  n_viruses = 100
)

# View results
power_result

# Try different methods
power_wilcoxon <- calc_virome_power(
  n_samples = 10,
  effect_size = 1.5,
  n_viruses = 100,
  method = "wilcoxon"
)

power_ttest <- calc_virome_power(
  n_samples = 10,
  effect_size = 1.5,
  n_viruses = 100,
  method = "t.test"
)

# Compare results
power_wilcoxon$power
power_ttest$power
```

### Estimate Sample Size

```r
# Basic sample size estimation
sample_size <- estimate_sample_size(
  power = 0.8,
  effect_size = 1.5,
  n_viruses = 100
)

# View results
sample_size

# Try different parameters
sample_size_high_power <- estimate_sample_size(
  power = 0.9,
  effect_size = 1.5,
  n_viruses = 100
)

sample_size_low_effect <- estimate_sample_size(
  power = 0.8,
  effect_size = 1.2,
  n_viruses = 100
)

# Compare results
sample_size$sample_size
sample_size_high_power$sample_size
sample_size_low_effect$sample_size
```

### Plot Power Curve

```r
# Basic power curve
power_curve <- plot_power_curve(
  effect_size = 1.5,
  n_viruses = 100,
  sample_sizes = seq(5, 30, by = 5)
)

# Display plot
print(power_curve)
```

### Generate HTML Report

```r
# Generate a report with the results from previous steps
report <- generate_power_report(
  power_results = power_result,
  sample_size_results = sample_size,
  power_curve = power_curve,
  output_file = "test_report.html"
)

# Check if the file was created
file.exists("test_report.html")
```

## 5. Package Checks

Run full package checks to ensure everything is working properly:

```r
devtools::check()
```

This will:
- Check documentation consistency
- Run all tests
- Verify package structure
- Check for CRAN submission issues

## 6. Documentation Testing

Check that all documentation is properly rendered:

```r
?simulate_virome_data
?calc_virome_power
?estimate_sample_size
?plot_power_curve
?generate_power_report
```

## 7. Building a Vignette

Build and view the vignette:

```r
devtools::build_vignettes()
browseVignettes("viromePower")
```

## 8. Installation Testing

Test installing the package locally:

```r
devtools::install()
```

Then restart R and try using it as a user would:

```r
library(viromePower)
# Run examples from the README
```

## 9. Troubleshooting

If you encounter issues:

1. Check error messages carefully
2. Use `browser()` for interactive debugging
3. Add more `message()` calls in functions to trace execution
4. Consider using `options(error = recover)` to debug errors 

## 10. Performance Testing

For larger simulations, test performance:

```r
system.time({
  power_result <- calc_virome_power(
    n_samples = 20,
    effect_size = 1.5,
    n_viruses = 1000,
    n_sim = 500
  )
})
```