#' Generate Viral Diversity Power Analysis Report
#'
#' Generates a comprehensive HTML report for viral diversity power analysis,
#' including sample size recommendations, power curves, and visualizations.
#'
#' @param n_samples Number of samples per group to analyze
#' @param effect_size Expected effect size between groups
#' @param n_viruses Number of viral taxa in the dataset
#' @param diversity_measure Type of diversity to analyze: "shannon", "simpson", "richness", "evenness", "chao1", "ace", "inv_simpson", "fisher_alpha", "bray", "jaccard", or "unifrac" (default: "shannon")
#' @param alpha Significance level (default: 0.05)
#' @param sparsity Proportion of zeros in the data (default: 0.8)
#' @param dispersion Dispersion parameter for viral abundance (default: 2)
#' @param n_sim Number of simulations for power calculation (default: 100)
#' @param output_file Path to save the HTML report (default: "diversity_power_report.html")
#' @param show_individual_sims Logical indicating whether to show individual simulation points (default: FALSE)
#'
#' @return Path to the generated HTML report
#' @export
#'
#' @examples
#' # Generate a power report for Shannon diversity
#' report_file <- generate_diversity_power_report(
#'   n_samples = 15, 
#'   effect_size = 1.2,
#'   n_viruses = 200,
#'   diversity_measure = "shannon"
#' )
#'
#' # Generate a power report for Inverse Simpson index
#' inv_simpson_report <- generate_diversity_power_report(
#'   n_samples = 15,
#'   effect_size = 1.0,
#'   n_viruses = 200,
#'   diversity_measure = "inv_simpson",
#'   output_file = "inv_simpson_power_report.html"
#' )
#'
#' # Generate a power report for ACE richness estimator
#' ace_report <- generate_diversity_power_report(
#'   n_samples = 18,
#'   effect_size = 1.5,
#'   n_viruses = 250,
#'   diversity_measure = "ace",
#'   output_file = "ace_power_report.html"
#' )
#'
#' # Generate a power report for Bray-Curtis dissimilarity
#' beta_report <- generate_diversity_power_report(
#'   n_samples = 20,
#'   effect_size = 0.15,
#'   n_viruses = 300,
#'   diversity_measure = "bray",
#'   output_file = "beta_diversity_power_report.html"
#' )
generate_diversity_power_report <- function(n_samples, effect_size, n_viruses,
                                          diversity_measure = "shannon",
                                          alpha = 0.05, sparsity = 0.8, dispersion = 2,
                                          n_sim = 100, output_file = "diversity_power_report.html",
                                          show_individual_sims = FALSE) {
  
  # Validate input parameters
  if (n_samples < 3) {
    stop("n_samples must be at least 3 per group")
  }
  if (effect_size <= 0) {
    stop("effect_size must be positive")
  }
  if (is.null(output_file) || output_file == "") {
    stop("Output file path must be specified")
  }
  
  # Determine if this is alpha or beta diversity
  is_beta <- diversity_measure %in% c("bray", "jaccard", "unifrac")
  
  # Run power analysis for the specified parameters
  power_result <- calc_viral_diversity_power(
    n_samples = n_samples,
    effect_size = effect_size,
    n_viruses = n_viruses,
    alpha = alpha,
    sparsity = sparsity,
    dispersion = dispersion,
    diversity_measure = diversity_measure,
    n_sim = n_sim
  )
  
  # Pre-compute and format values for the report
  power_value <- sprintf("%.2f", power_result$power * 100)
  recommended_value <- FALSE
  if (power_result$power < 0.8) {
    # If power is less than 80%, generate a curve to find recommended sample size
    # Only do this for sample sizes from current n_samples up to 50 or double current, whichever is smaller
    max_samples <- min(50, n_samples * 2)
    sample_range <- seq(n_samples, max_samples, by = max(1, ceiling((max_samples - n_samples) / 10)))
    
    # Calculate power for each sample size
    sample_powers <- numeric(length(sample_range))
    for (i in seq_along(sample_range)) {
      temp_result <- calc_viral_diversity_power(
        n_samples = sample_range[i],
        effect_size = effect_size,
        n_viruses = n_viruses,
        alpha = alpha,
        sparsity = sparsity,
        dispersion = dispersion,
        diversity_measure = diversity_measure,
        n_sim = max(30, floor(n_sim / 2))  # Use fewer simulations for efficiency
      )
      sample_powers[i] <- temp_result$power
    }
    
    # Find the minimum sample size that achieves at least 80% power
    recommended_sample_size <- NA
    for (i in seq_along(sample_range)) {
      if (sample_powers[i] >= 0.8) {
        recommended_sample_size <- sample_range[i]
        break
      }
    }
    
    if (!is.na(recommended_sample_size)) {
      recommended_value <- recommended_sample_size
    }
  }
  
  # Format diversity measure name for display
  diversity_title <- switch(diversity_measure,
                           shannon = "Shannon Diversity",
                           simpson = "Simpson Diversity",
                           inv_simpson = "Inverse Simpson Diversity",
                           richness = "Species Richness",
                           evenness = "Pielou's Evenness",
                           chao1 = "Chao1 Richness Estimator",
                           ace = "ACE Richness Estimator",
                           fisher_alpha = "Fisher's Alpha",
                           bray = "Bray-Curtis Dissimilarity",
                           jaccard = "Jaccard Distance",
                           unifrac = "UniFrac Distance",
                           "Diversity")
  
  # Generate sample size curve
  sample_curve <- plot_diversity_power_curve(
    param_range = seq(max(3, floor(n_samples * 0.5)), ceiling(n_samples * 1.5), 
                     length.out = 6),
    param_type = "n_samples",
    effect_size = effect_size,
    n_viruses = n_viruses,
    diversity_measure = diversity_measure,
    alpha = alpha,
    sparsity = sparsity,
    dispersion = dispersion,
    n_sim = max(30, floor(n_sim / 2)),  # Use fewer simulations for efficiency
    show_individual_sims = show_individual_sims
  )
  
  # Generate effect size curve
  # For beta diversity, effect sizes are typically smaller
  if (is_beta) {
    effect_range <- seq(max(0.05, effect_size * 0.5), effect_size * 1.5, length.out = 6)
  } else {
    effect_range <- seq(max(0.3, effect_size * 0.5), effect_size * 1.5, length.out = 6)
  }
  
  effect_curve <- plot_diversity_power_curve(
    param_range = effect_range,
    param_type = "effect_size",
    n_samples = n_samples,
    n_viruses = n_viruses,
    diversity_measure = diversity_measure,
    alpha = alpha,
    sparsity = sparsity,
    dispersion = dispersion,
    n_sim = max(30, floor(n_sim / 2)),  # Use fewer simulations for efficiency
    show_individual_sims = show_individual_sims
  )
  
  # Save plots for the report
  temp_dir <- tempdir()
  sample_curve_file <- file.path(temp_dir, "sample_curve.png")
  effect_curve_file <- file.path(temp_dir, "effect_curve.png")
  
  ggplot2::ggsave(sample_curve_file, sample_curve, width = 8, height = 6, dpi = 100)
  ggplot2::ggsave(effect_curve_file, effect_curve, width = 8, height = 6, dpi = 100)
  
  # Generate example simulation to visualize the data
  set.seed(123)  # For reproducibility
  sim_data <- simulate_virome_data(
    n_samples = n_samples * 2,  # Both groups
    n_viruses = n_viruses,
    sparsity = sparsity,
    dispersion = dispersion,
    effect_size = effect_size
  )
  
  # Calculate diversity values for visualization
  counts <- sim_data$counts
  metadata <- sim_data$metadata
  counts_t <- t(counts)
  
  # Calculate diversity values based on the chosen measure
  diversity_values <- numeric(nrow(metadata))
  
  if (!is_beta) {
    # Alpha diversity calculation
    for (j in 1:nrow(metadata)) {
      sample_counts <- counts_t[j, ]
      nonzero_counts <- sample_counts[sample_counts > 0]
      
      if (length(nonzero_counts) == 0) {
        diversity_values[j] <- 0
        next
      }
      
      relative_abundance <- nonzero_counts / sum(nonzero_counts)
      
      if (diversity_measure == "shannon") {
        diversity_values[j] <- -sum(relative_abundance * log(relative_abundance))
      } else if (diversity_measure == "simpson") {
        diversity_values[j] <- 1 - sum(relative_abundance^2)
      } else if (diversity_measure == "richness") {
        diversity_values[j] <- length(nonzero_counts)
      } else if (diversity_measure == "evenness") {
        shannon <- -sum(relative_abundance * log(relative_abundance))
        richness <- length(nonzero_counts)
        diversity_values[j] <- if(richness > 1) shannon / log(richness) else 0
      }
    }
    
    # Create data frame for diversity visualization
    diversity_df <- data.frame(
      sample_id = metadata$sample_id,
      group = metadata$group,
      diversity = diversity_values
    )
    
    # Create alpha diversity boxplot
    diversity_plot <- ggplot2::ggplot(diversity_df, ggplot2::aes(x = group, y = diversity, fill = group)) +
      ggplot2::geom_boxplot(alpha = 0.7) +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.5) +
      ggplot2::labs(title = paste("Simulated", diversity_title, "Values by Group"),
                   x = "Group", y = diversity_title) +
      ggplot2::theme_minimal() +
      ggplot2::scale_fill_brewer(palette = "Set1")
    
  } else {
    # Beta diversity visualization (NMDS or PCoA plot)
    if (diversity_measure == "bray") {
      dist_matrix <- as.matrix(vegan::vegdist(counts_t, method = "bray"))
    } else if (diversity_measure == "jaccard") {
      binary_counts <- ifelse(counts_t > 0, 1, 0)
      dist_matrix <- as.matrix(vegan::vegdist(binary_counts, method = "jaccard"))
    } else if (diversity_measure == "unifrac") {
      binary_counts <- ifelse(counts_t > 0, 1, 0)
      dist_matrix <- as.matrix(vegan::vegdist(binary_counts, method = "jaccard"))
      
      for (j in 1:nrow(dist_matrix)) {
        for (k in 1:ncol(dist_matrix)) {
          if (j != k) {
            abundance_weight <- sum(abs(counts_t[j,] - counts_t[k,])) / sum(counts_t[j,] + counts_t[k,])
            dist_matrix[j,k] <- dist_matrix[j,k] * abundance_weight
          }
        }
      }
    }
    
    # Perform NMDS ordination
    nmds_result <- tryCatch({
      vegan::metaMDS(dist_matrix, k = 2, trymax = 50)
    }, error = function(e) {
      # If NMDS fails, use PCoA as fallback
      pcoa <- stats::cmdscale(dist_matrix, k = 2, eig = TRUE)
      list(points = pcoa$points, stress = NA)
    })
    
    # Create data frame for ordination plot
    if (!is.null(nmds_result$points)) {
      ord_df <- data.frame(
        NMDS1 = nmds_result$points[,1],
        NMDS2 = nmds_result$points[,2],
        group = metadata$group
      )
      
      # Create ordination plot
      diversity_plot <- ggplot2::ggplot(ord_df, ggplot2::aes(x = NMDS1, y = NMDS2, color = group)) +
        ggplot2::geom_point(size = 3, alpha = 0.7) +
        ggplot2::stat_ellipse(type = "norm", level = 0.95) +
        ggplot2::labs(
          title = paste("NMDS Ordination based on", diversity_title),
          subtitle = if(!is.na(nmds_result$stress)) paste("Stress:", round(nmds_result$stress, 3)) else "PCoA Ordination"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::scale_color_brewer(palette = "Set1")
    } else {
      # Fallback if ordination fails
      diversity_plot <- ggplot2::ggplot() +
        ggplot2::annotate("text", x = 0, y = 0, label = "Ordination failed") +
        ggplot2::theme_minimal()
    }
  }
  
  # Save diversity visualization
  diversity_plot_file <- file.path(temp_dir, "diversity_plot.png")
  ggplot2::ggsave(diversity_plot_file, diversity_plot, width = 8, height = 6, dpi = 100)
  
  # Prepare report content
  report_content <- paste0(
    "---\n",
    "title: \"", "Viral ", diversity_title, " Power Analysis Report", "\"\n",
    "date: \"", format(Sys.time(), "%Y-%m-%d"), "\"\n",
    "output: html_document\n",
    "---\n\n",
    
    "```{r setup, include=FALSE}\n",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)\n",
    "library(knitr)\n",
    "```\n\n",
    
    "## Study Parameters\n\n",
    
    "```{r parameters}\n",
    "parameters <- data.frame(\n",
    "  Parameter = c('Diversity Measure', 'Number of Samples per Group', 'Effect Size',\n",
    "                 'Number of Viral Taxa', 'Significance Level (α)', 'Sparsity', 'Dispersion', 'Simulations'),\n",
    "  Value = c('", diversity_title, "', ", n_samples, ", ", effect_size, ", ", 
                 n_viruses, ", ", alpha, ", ", sparsity, ", ", dispersion, ", ", n_sim, ")\n",
    ")\n",
    "knitr::kable(parameters, format = 'html')\n",
    "```\n\n",
    
    "## Power Analysis Results\n\n",
    
    "```{r power_result}\n",
    "power_pct <- ", power_value, "\n",
    "result_color <- if(power_pct >= 80) 'green' else if(power_pct >= 60) 'orange' else 'red'\n",
    "```\n\n",
    
    "<div style=\"background-color: #f8f8f8; padding: 15px; border-radius: 5px; margin-bottom: 20px;\">\n",
    "  <h3>Statistical Power</h3>\n",
    "  <p>The estimated statistical power for detecting ", if(is_beta) "community differences" else "diversity differences", " with the specified parameters is:</p>\n",
    "  <h2 style=\"color: `r result_color`\">`r power_pct`% Power</h2>\n",
    "  ", if(is.logical(recommended_value) && !recommended_value) "" else paste0("<p><strong>Recommendation:</strong> To achieve at least 80% power, you should use ", recommended_value, " samples per group.</p>"), "\n",
    "</div>\n\n",
    
    "## Power by Sample Size\n\n",
    
    "This curve shows how statistical power changes with different sample sizes, keeping all other parameters constant.\n\n",
    
    "```{r sample_curve, fig.cap='Power curve by sample size'}\n",
    "knitr::include_graphics('", sample_curve_file, "')\n",
    "```\n\n",
    
    "## Power by Effect Size\n\n",
    
    "This curve shows how statistical power changes with different effect sizes, keeping all other parameters constant.\n\n",
    
    "```{r effect_curve, fig.cap='Power curve by effect size'}\n",
    "knitr::include_graphics('", effect_curve_file, "')\n",
    "```\n\n",
    
    "## Simulated Data Visualization\n\n",
    
    "This plot shows simulated data for the specified parameters, illustrating the expected ", if(is_beta) "community differences" else "diversity differences", " between groups.\n\n",
    
    "```{r diversity_viz, fig.cap='Visualization of simulated data'}\n",
    "knitr::include_graphics('", diversity_plot_file, "')\n",
    "```\n\n",
    
    "## Interpretation\n\n",
    
    "```{r interpretation}\n",
    "power_interpretation <- if(power_pct >= 80) {\n",
    "  paste('Your study design has <span style=\"color:green;\">**adequate power**</span> (≥80%) to detect', \n",
    "        'the specified effect size. This means there is a high probability of detecting',\n",
    "        'true differences if they exist.')\n",
    "} else if(power_pct >= 60) {\n",
    "  paste('Your study design has <span style=\"color:orange;\">**marginal power**</span> to detect',\n",
    "        'the specified effect size. Consider increasing sample size if possible.')\n",
    "} else {\n",
    "  paste('Your study design has <span style=\"color:red;\">**low power**</span> to detect',\n",
    "        'the specified effect size. It is recommended to increase sample size',\n",
    "        'or focus on larger effect sizes.')\n",
    "}\n",
    "```\n\n",
    
    "### Power Interpretation\n\n",
    
    "`r power_interpretation`\n\n",
    
    "### Statistical Background\n\n",
    
    if(is_beta) {
      paste(
        "Beta diversity measures examine differences in viral community composition between samples or groups.",
        "The analysis shown uses", diversity_title, "to calculate distances between samples, and",
        "PERMANOVA (adonis) tests to determine if groups have statistically different community structures.\n\n",
        
        "Effect sizes for beta diversity are typically expressed as R² values from PERMANOVA,",
        "representing the proportion of variance explained by group differences.",
        "Values around 0.1-0.2 are common in microbiome studies, with values over 0.3 indicating strong effects.\n\n"
      )
    } else {
      paste(
        "Alpha diversity measures quantify the diversity within each sample.",
        "The", diversity_title, "index captures", 
        ifelse(diversity_measure == "shannon", 
               "both richness and evenness, giving more weight to rare taxa.",
               ifelse(diversity_measure == "simpson", 
                      "both richness and evenness, giving more weight to abundant taxa.",
                      ifelse(diversity_measure == "inv_simpson", 
                             "both richness and evenness, emphasizing community evenness and rare taxa more than Simpson.",
                             ifelse(diversity_measure == "richness", 
                                    "the number of different viral taxa present in a sample.",
                                    ifelse(diversity_measure == "evenness", 
                                           "how equally abundant the viral taxa are in a sample.",
                                           ifelse(diversity_measure == "chao1", 
                                                  "richness with a correction for unobserved rare species based on singleton and doubleton counts.",
                                                  ifelse(diversity_measure == "ace", 
                                                         "richness using abundance data to estimate rare species that may have been missed in sampling.",
                                                         ifelse(diversity_measure == "fisher_alpha", 
                                                                "the relationship between species and individuals using a logarithmic distribution model.",
                                                                "aspects of viral community structure within samples.")))))))),
        "\n\n",
        
        "Effect sizes for alpha diversity are typically expressed as standardized mean differences (Cohen's d).",
        "Values around 0.5 indicate moderate effects, while values over 0.8 are considered large effects.\n\n"
      )
    },
    
    "## Conclusions\n\n",
    
    if(power_value >= 80) {
      paste(
        "**The planned study design is likely to have sufficient statistical power.**\n\n",
        
        "With", n_samples, "samples per group, you can reliably detect differences in", diversity_title,
        "with an effect size of", effect_size, "or larger."
      )
    } else if(is.logical(recommended_value) && !recommended_value) {
      paste(
        "**The planned study design may have insufficient statistical power.**\n\n",
        
        "With", n_samples, "samples per group, you may struggle to reliably detect differences in", diversity_title,
        "with an effect size of", effect_size, ". Consider one of the following options:\n\n",
        
        "1. Increase sample size beyond the tested range if possible\n",
        "2. Focus on detecting larger effect sizes\n",
        "3. Use more sensitive analytical methods\n",
        "4. Consider longitudinal sampling to increase statistical power"
      )
    } else {
      paste(
        "**The planned study should be adjusted to achieve adequate statistical power.**\n\n",
        
        "Based on the simulations, **increasing to", recommended_value, "samples per group** would",
        "provide sufficient power (≥80%) to detect differences in", diversity_title,
        "with an effect size of", effect_size, "."
      )
    },
    
    "\n\n---\n\n",
    "*Report generated by the viromePower package*"
  )
  
  # Create temporary Rmd file
  report_rmd <- tempfile(fileext = ".Rmd")
  writeLines(report_content, report_rmd)
  
  # Render the report
  rmarkdown::render(
    input = report_rmd,
    output_file = basename(output_file),
    output_dir = dirname(output_file),
    quiet = TRUE
  )
  
  # Return the path to the generated report
  return(output_file)
}