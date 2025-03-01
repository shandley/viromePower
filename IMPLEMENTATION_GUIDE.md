# viromePower Implementation Guide

This document provides guidelines for implementing the statistical methods in the viromePower package. Currently, the package contains function templates with placeholders where the actual statistical calculations should go.

## Statistical Foundations for Virome Power Analysis

When implementing the full statistical methods, consider these virome-specific characteristics:

1. **High Sparsity**: Virome data typically contains many zeros (often >80%)
2. **Low Abundance**: Viral counts are generally lower than bacterial counts
3. **High Diversity**: Virome samples can contain thousands of viral taxa
4. **High Variability**: Viral abundances can vary widely between samples
5. **Compositional Data**: Sequence count data is compositional by nature
6. **Multiple Testing**: Many tests are performed simultaneously

## Implementation Roadmap

### 1. `simulate_virome_data.R`

This function should simulate realistic virome count data:

```r
simulate_virome_data <- function(n_samples, n_viruses, sparsity = 0.8, 
                                dispersion = 2, effect_size = 1.5, 
                                groups = NULL) {
  # Implementation:
  
  # 1. Generate mean abundances for each virus (often follows a power law)
  abundance_means <- rexp(n_viruses, rate = 0.1)
  abundance_means <- abundance_means / sum(abundance_means) * 1e6  # Scale to reasonable counts
  
  # 2. Apply sparsity (introduce zeros)
  sparsity_probs <- rbeta(n_viruses, 1, 4) * sparsity  # Different sparsity per virus
  
  # 3. Generate count matrix using negative binomial distribution
  # with appropriate dispersion parameter
  counts <- matrix(0, nrow = n_viruses, ncol = n_samples)
  
  # 4. Split samples into groups and apply effect size difference
  groups <- groups %||% rep(c("A", "B"), each = n_samples/2)
  
  # 5. Generate counts with appropriate distribution and effect sizes
  # between groups
  
  # Return data
  list(
    counts = counts,
    metadata = data.frame(
      sample_id = paste0("sample_", 1:n_samples),
      group = groups
    )
  )
}
```

### 2. `calc_virome_power.R`

This function should calculate power based on simulations:

```r
calc_virome_power <- function(n_samples, effect_size, n_viruses,
                             alpha = 0.05, sparsity = 0.8, dispersion = 2,
                             method = "wilcoxon", n_sim = 100) {
  # Implementation:
  
  # 1. Run multiple simulations
  power_results <- numeric(n_sim)
  
  for (i in 1:n_sim) {
    # 2. Simulate data with the given parameters
    sim_data <- simulate_virome_data(
      n_samples = n_samples * 2,  # Total samples (both groups)
      n_viruses = n_viruses,
      sparsity = sparsity,
      dispersion = dispersion,
      effect_size = effect_size
    )
    
    # 3. Apply statistical tests based on method parameter
    # (wilcoxon, t.test, deseq, etc.)
    
    # 4. Apply multiple testing correction
    
    # 5. Calculate power as proportion of true positives detected
  }
  
  # Return results
  list(
    power = mean(power_results),
    parameters = list(
      n_samples = n_samples,
      effect_size = effect_size,
      n_viruses = n_viruses,
      alpha = alpha,
      sparsity = sparsity,
      method = method
    )
  )
}
```

### 3. `estimate_sample_size.R`

This function should estimate required sample size:

```r
estimate_sample_size <- function(power = 0.8, effect_size, n_viruses,
                               alpha = 0.05, sparsity = 0.8, dispersion = 2,
                               method = "wilcoxon") {
  # Implementation:
  
  # 1. Define a range of sample sizes to try
  sample_sizes <- seq(5, 50, by = 5)
  
  # 2. Calculate power for each sample size
  powers <- numeric(length(sample_sizes))
  
  for (i in seq_along(sample_sizes)) {
    # 3. Calculate power for this sample size
    result <- calc_virome_power(
      n_samples = sample_sizes[i],
      effect_size = effect_size,
      n_viruses = n_viruses,
      alpha = alpha,
      sparsity = sparsity,
      dispersion = dispersion,
      method = method,
      n_sim = 50  # Fewer simulations for efficiency
    )
    
    powers[i] <- result$power
  }
  
  # 4. Find the minimum sample size that achieves the desired power
  # (could use interpolation for better precision)
  required_n <- min(sample_sizes[powers >= power], 
                   Inf)
  
  # Return results
  list(
    sample_size = required_n,
    parameters = list(
      power = power,
      effect_size = effect_size,
      n_viruses = n_viruses,
      alpha = alpha,
      sparsity = sparsity,
      method = method
    )
  )
}
```

### 4. `plot_power_curve.R` 

This function should create a power curve visualization:

```r
plot_power_curve <- function(effect_size, n_viruses, 
                           sample_sizes = seq(5, 30, by = 5),
                           alpha = 0.05, sparsity = 0.8, dispersion = 2,
                           method = "wilcoxon") {
  # Implementation:
  
  # 1. Calculate power for each sample size
  powers <- numeric(length(sample_sizes))
  
  for (i in seq_along(sample_sizes)) {
    # 2. Calculate power for this sample size
    result <- calc_virome_power(
      n_samples = sample_sizes[i],
      effect_size = effect_size,
      n_viruses = n_viruses,
      alpha = alpha,
      sparsity = sparsity,
      dispersion = dispersion,
      method = method,
      n_sim = 50  # Fewer simulations for efficiency
    )
    
    powers[i] <- result$power
  }
  
  # 3. Create data frame for plotting
  plot_data <- data.frame(
    sample_size = sample_sizes,
    power = powers
  )
  
  # 4. Create plot
  ggplot2::ggplot(plot_data, ggplot2::aes(x = sample_size, y = power)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    ggplot2::labs(
      x = "Sample Size per Group",
      y = "Statistical Power",
      title = "Power Curve for Virome Study",
      subtitle = paste("Effect size:", effect_size, 
                      "| Number of viruses:", n_viruses)
    ) +
    ggplot2::theme_minimal()
}
```

### 5. `generate_power_report.R`

This function should create an HTML report:

```r
generate_power_report <- function(power_results, sample_size_results, power_curve,
                              output_file = "virome_power_report.html",
                              title = "Virome Study Power Analysis Report",
                              include_code = FALSE) {
  # Implementation:
  
  # 1. Create a temporary directory for the report
  temp_dir <- tempdir()
  temp_file <- file.path(temp_dir, "report.Rmd")
  
  # 2. Get the template
  template <- system.file("templates", "power_report_template.Rmd", 
                         package = "viromePower")
  
  # 3. Read the template
  rmd_content <- readLines(template)
  
  # 4. Replace placeholders with actual values
  rmd_content <- gsub("{{title}}", title, rmd_content)
  rmd_content <- gsub("{{include_code}}", tolower(include_code), rmd_content)
  rmd_content <- gsub("{{n_viruses}}", power_results$parameters$n_viruses, rmd_content)
  rmd_content <- gsub("{{effect_size}}", power_results$parameters$effect_size, rmd_content)
  # Add more replacements for other parameters
  
  # 5. Write the modified template to a temporary file
  writeLines(rmd_content, temp_file)
  
  # 6. Render the report
  rmarkdown::render(
    input = temp_file,
    output_file = basename(output_file),
    output_dir = dirname(output_file),
    quiet = TRUE
  )
  
  # Return output file path
  invisible(output_file)
}
```

## Statistical Methods for Virome Analysis

When implementing the full statistical methods, consider these approaches:

1. **Zero-Inflation Handling**:
   - Zero-inflated negative binomial models
   - Zero-imputation methods based on taxonomic relationships
   - Separate modeling of presence/absence and abundance

2. **Compositional Data Analysis**:
   - Centered log-ratio (CLR) transformations
   - Additive log-ratio (ALR) transformations
   - Isometric log-ratio (ILR) transformations

3. **Multiple Testing Correction**:
   - Benjamini-Hochberg procedure
   - Benjamini-Yekutieli procedure
   - Permutation-based FDR control

4. **Statistical Tests**:
   - Wilcoxon rank-sum test (robust to non-normality)
   - DESeq2-like negative binomial models
   - ALDEx2-like Bayesian approaches

## Recommended Resources

For implementing these statistical methods, consider these resources:

1. McMurdie, P.J. & Holmes, S. (2014). Waste not, want not: why rarefying microbiome data is inadmissible. PLoS Computational Biology, 10(4), e1003531.

2. Love, M.I., Huber, W., & Anders, S. (2014). Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biology, 15(12), 550.

3. Fernandes, A.D., et al. (2013). Unifying the analysis of high-throughput sequencing datasets: characterizing RNA-seq, 16S rRNA gene sequencing and selective growth experiments by compositional data analysis. Microbiome, 2(1), 15.

4. Thorsen, J., et al. (2016). Large-scale benchmarking reveals false discoveries and count transformation sensitivity in 16S rRNA gene amplicon data analysis methods used in microbiome studies. Microbiome, 4(1), 62.