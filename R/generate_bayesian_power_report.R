#' Generate Bayesian Power Analysis Report
#'
#' Creates an HTML report with visualizations and summaries of Bayesian power analysis results.
#' This report includes Bayesian-specific metrics like posterior probabilities and Bayes factors,
#' along with credible intervals and effect size distributions.
#'
#' @param bayesian_power_results List of results from calc_bayesian_power()
#' @param output_file File path for the output HTML report
#' @param title Title for the report (default: "Bayesian Virome Power Analysis Report")
#' @param include_code Boolean indicating whether to include R code in the report (default: FALSE)
#' @param custom_css Optional custom CSS for the report styling
#'
#' @return The path to the generated HTML report
#' @export
#'
#' @examples
#' # Calculate Bayesian power
#' power_result <- calc_bayesian_power(
#'   n_samples = 20,
#'   effect_size = 3.0,
#'   n_viruses = 50,
#'   sparsity = 0.6,
#'   prior_strength = 2.5,
#'   n_sim = 50
#' )
#'
#' # Generate report
#' report_path <- generate_bayesian_power_report(
#'   bayesian_power_results = power_result,
#'   output_file = "bayesian_power_report.html",
#'   title = "Bayesian Virome Power Analysis"
#' )
generate_bayesian_power_report <- function(bayesian_power_results,
                                         output_file = "bayesian_power_report.html",
                                         title = "Bayesian Virome Power Analysis Report",
                                         include_code = FALSE,
                                         custom_css = NULL) {
  
  # Validate input
  if (!is.list(bayesian_power_results) || 
      !all(c("power", "expected_discoveries", "parameters") %in% names(bayesian_power_results))) {
    stop("bayesian_power_results must be the output from calc_bayesian_power()")
  }
  
  # Check for essential packages
  essential_packages <- c("rmarkdown", "ggplot2", "knitr")
  missing_essential <- essential_packages[!sapply(essential_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_essential) > 0) {
    stop(paste0(
      "Missing essential packages for report generation: ", 
      paste(missing_essential, collapse = ", "), 
      ". Please install them with: install.packages(c('", 
      paste(missing_essential, collapse = "', '"), 
      "'))"
    ))
  }
  
  # Check for optional packages and set flags for conditional features
  has_kableExtra <- requireNamespace("kableExtra", quietly = TRUE)
  has_dplyr <- requireNamespace("dplyr", quietly = TRUE)
  has_tidyr <- requireNamespace("tidyr", quietly = TRUE)
  
  # Warn about optional packages
  optional_missing <- c()
  if (!has_kableExtra) optional_missing <- c(optional_missing, "kableExtra")
  if (!has_dplyr) optional_missing <- c(optional_missing, "dplyr")
  if (!has_tidyr) optional_missing <- c(optional_missing, "tidyr")
  
  if (length(optional_missing) > 0) {
    warning(paste0(
      "Some optional packages for enhanced report formatting are missing: ", 
      paste(optional_missing, collapse = ", "), 
      ". Reports will still be generated with basic formatting. ",
      "For enhanced formatting, install the missing packages with: ",
      "install.packages(c('", paste(optional_missing, collapse = "', '"), "'))"
    ))
  }
  
  # Create a temporary Rmd file with the report content
  rmd_file <- tempfile(fileext = ".Rmd")
  
  # Generate the R markdown content
  rmd_content <- generate_bayesian_report_rmd(
    bayesian_power_results = bayesian_power_results,
    title = title,
    include_code = include_code,
    custom_css = custom_css,
    output_file = output_file,
    has_kableExtra = has_kableExtra,
    has_dplyr = has_dplyr,
    has_tidyr = has_tidyr
  )
  
  # Write the content to the temporary file
  writeLines(rmd_content, rmd_file)
  
  # Render the report using rmarkdown
  tryCatch({
    # Pass the data directly to the rendering environment
    env <- new.env()
    env$bayesian_power_results <- bayesian_power_results
    
    rmarkdown::render(
      input = rmd_file,
      output_file = basename(output_file),
      output_dir = dirname(output_file),
      quiet = TRUE,
      envir = env
    )
  }, error = function(e) {
    stop(paste0("Error rendering report: ", e$message))
  })
  
  # Return the path to the generated report
  return(normalizePath(output_file))
}

#' Generate R Markdown content for Bayesian power report
#'
#' Helper function to create the R markdown content for the Bayesian power report
#'
#' @param bayesian_power_results List of results from calc_bayesian_power()
#' @param title Title for the report
#' @param include_code Boolean indicating whether to include R code
#' @param custom_css Optional custom CSS for styling
#' @param output_file File path for the output HTML report
#' @param has_kableExtra Boolean indicating whether kableExtra is available
#' @param has_dplyr Boolean indicating whether dplyr is available
#' @param has_tidyr Boolean indicating whether tidyr is available
#'
#' @return A character string with the R markdown content
#' @keywords internal
generate_bayesian_report_rmd <- function(bayesian_power_results,
                                       title,
                                       include_code,
                                       custom_css,
                                       output_file,
                                       has_kableExtra = TRUE,
                                       has_dplyr = TRUE,
                                       has_tidyr = TRUE) {
  
  # Keep local reference for building the report
  power_data <- bayesian_power_results
  
  # Create the header section
  header <- c(
    "---",
    paste0("title: \"", title, "\""),
    "date: \"`r format(Sys.time(), '%B %d, %Y')`\"",
    "output:",
    "  html_document:",
    "    theme: cosmo",
    "    toc: true",
    "    toc_float: true",
    "    highlight: tango",
    "    code_folding: hide",
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)",
    "library(ggplot2)",
    "library(knitr)",
    if (has_dplyr) "library(dplyr)" else "# dplyr not available",
    if (has_tidyr) "library(tidyr)" else "# tidyr not available",
    if (has_kableExtra) "library(kableExtra)" else "# kableExtra not available",
    "```",
    "",
    "```{css, echo=FALSE}",
    "/* Custom styling */",
    ".main-container {",
    "  max-width: 1200px;",
    "  margin-left: auto;",
    "  margin-right: auto;",
    "}",
    ".power-summary {",
    "  padding: 15px;",
    "  background-color: #f8f9fa;",
    "  border-radius: 5px;",
    "  margin-bottom: 20px;",
    "}",
    ".power-value {",
    "  font-size: 24px;",
    "  font-weight: bold;",
    "  color: #2c3e50;",
    "}",
    ".parameter-table th {",
    "  background-color: #eaecef;",
    "}",
    ".bayesian-explanation {",
    "  background-color: #f5f5f5;",
    "  border-left: 4px solid #3498db;",
    "  padding: 15px;",
    "  margin-bottom: 20px;",
    "}",
    if (!is.null(custom_css)) custom_css else "",
    "```",
    ""
  )
  
  # Get the data setup code
  data_setup <- c(
    "```{r data-setup, include=FALSE}",
    "# Make the input data available to all chunks",
    "power_data <- bayesian_power_results",
    "```",
    ""
  )
  
  # Introduction section
  intro <- c(
    "## Overview {.tabset}",
    "",
    "### Summary",
    "",
    "<div class=\"power-summary\">",
    paste0("<p>Bayesian statistical power: <span class=\"power-value\">", 
           round(power_data$power * 100, 1), "%</span></p>"),
    paste0("<p>Expected discoveries: <span class=\"power-value\">", 
           round(power_data$expected_discoveries, 1), 
           "</span> out of <span class=\"power-value\">", 
           power_data$parameters$n_viruses, "</span> viral taxa</p>"),
    paste0("<p>False discovery proportion: <span class=\"power-value\">", 
           round(power_data$false_discovery_proportion * 100, 1), "%</span></p>"),
    "</div>",
    "",
    "This report presents a comprehensive Bayesian power analysis for a virome study design. ",
    "The Bayesian approach incorporates prior knowledge and provides more nuanced decision-making ",
    "by quantifying the evidence for differential abundance rather than solely relying on p-values.",
    "",
    "<div class=\"bayesian-explanation\">",
    "<strong>Why Bayesian Power Analysis?</strong><br>",
    "Bayesian methods offer several advantages for virome studies:",
    "<ul>",
    "<li>Better handling of sparse data through incorporation of prior knowledge</li>",
    "<li>Direct probability statements about effect existence</li>",
    "<li>Credible intervals that have more intuitive interpretation than confidence intervals</li>",
    "<li>No need for multiple testing correction as in frequentist approaches</li>",
    "<li>Ability to incorporate prior knowledge from previous virome studies</li>",
    "</ul>",
    "</div>",
    "",
    "### Study Design",
    "",
    "```{r study-design, echo=FALSE}",
    "# Create a data frame for the study parameters",
    "params <- power_data$parameters",
    "param_df <- data.frame(",
    "  Parameter = c(",
    "    'Number of Samples per Group',",
    "    'Effect Size (Fold Change)',",
    "    'Number of Viral Taxa',",
    "    'Data Sparsity',",
    "    'Dispersion Parameter',",
    "    'Prior Strength',",
    "    'Credible Interval Width',",
    "    'Fold Change Threshold',",
    "    'Posterior Probability Threshold',",
    "    'Number of Simulations'",
    "  ),",
    "  Value = c(",
    "    params$n_samples,",
    "    params$effect_size,",
    "    params$n_viruses,",
    "    params$sparsity,",
    "    params$dispersion,",
    "    params$prior_strength,",
    "    params$credible_interval,",
    "    params$fold_change_threshold,",
    "    params$posterior_prob_threshold,",
    "    params$n_sim",
    "  )",
    ")",
    "",
    "# Display the parameter table",
    if (has_kableExtra) {
      c(
        "kable(param_df, col.names = c('Parameter', 'Value'),",
        "      caption = 'Study Design Parameters') %>%",
        "  kable_styling(bootstrap_options = c('striped', 'hover'), full_width = FALSE,",
        "                position = 'center') %>%",
        "  column_spec(1, bold = TRUE)"
      )
    } else {
      c(
        "kable(param_df, col.names = c('Parameter', 'Value'),",
        "      caption = 'Study Design Parameters')"
      )
    },
    "```",
    "",
    "This study design has an estimated Bayesian power of ", 
    paste0("**", round(power_data$power * 100, 1), "%**"),
    " for detecting viral taxa with a fold change of at least ", 
    power_data$parameters$fold_change_threshold, 
    " between groups, using a posterior probability threshold of ", 
    power_data$parameters$posterior_prob_threshold, 
    " and a ", paste0(power_data$parameters$credible_interval * 100, "%"), 
    " credible interval.",
    "",
    "## Bayesian Power Analysis Results",
    "",
    "### Posterior Probability Distribution",
    "",
    "```{r posterior-prob, echo=FALSE, fig.width=8, fig.height=5}",
    "# Extract posterior probabilities from all simulations",
    "posterior_probs <- unlist(lapply(power_data$simulation_results, function(x) x$posterior_probs))",
    "posterior_df <- data.frame(posterior_probability = posterior_probs)",
    "",
    "# Create histogram of posterior probabilities",
    "ggplot(posterior_df, aes(x = posterior_probability)) +",
    "  geom_histogram(bins = 30, fill = '#4292c6', color = 'black', alpha = 0.8) +",
    "  geom_vline(xintercept = power_data$parameters$posterior_prob_threshold, ",
    "             linetype = 'dashed', color = 'red', size = 1) +",
    "  annotate('text', x = power_data$parameters$posterior_prob_threshold - 0.05, ",
    "           y = max(hist(posterior_probs, breaks = 30, plot = FALSE)$counts) * 0.8, ",
    "           label = paste0('Threshold: ', power_data$parameters$posterior_prob_threshold), ",
    "           color = 'red', hjust = 1) +",
    "  labs(title = 'Distribution of Posterior Probabilities',",
    "       x = 'Posterior Probability of Effect',",
    "       y = 'Frequency') +",
    "  theme_minimal() +",
    "  theme(plot.title = element_text(hjust = 0.5, size = 16),",
    "        axis.title = element_text(size = 14),",
    "        axis.text = element_text(size = 12))",
    "```",
    "",
    "The histogram shows the distribution of posterior probabilities from the simulations. ",
    "The dashed red line indicates the threshold (", 
    power_data$parameters$posterior_prob_threshold, 
    ") used to determine statistical significance. ",
    "Higher posterior probabilities indicate stronger evidence for true effects.",
    "",
    "### Bayes Factor Distribution",
    "",
    "```{r bayes-factor, echo=FALSE, fig.width=8, fig.height=5}",
    "# Extract Bayes factors from all simulations",
    "bayes_factors <- unlist(lapply(power_data$simulation_results, function(x) x$bayes_factors))",
    "# Cap extremely large values for better visualization",
    "bayes_factors <- pmin(bayes_factors, 100)",
    "bayes_df <- data.frame(bayes_factor = bayes_factors)",
    "",
    "# Create histogram of Bayes factors",
    "ggplot(bayes_df, aes(x = bayes_factor)) +",
    "  geom_histogram(bins = 30, fill = '#7fbc41', color = 'black', alpha = 0.8) +",
    "  geom_vline(xintercept = 10, linetype = 'dashed', color = 'red', size = 1) +",
    "  geom_vline(xintercept = 3, linetype = 'dotted', color = 'orange', size = 1) +",
    "  annotate('text', x = 10, y = max(hist(bayes_factors, breaks = 30, plot = FALSE)$counts) * 0.9, ",
    "           label = 'Strong evidence', color = 'red', hjust = -0.1) +",
    "  annotate('text', x = 3, y = max(hist(bayes_factors, breaks = 30, plot = FALSE)$counts) * 0.8, ",
    "           label = 'Moderate evidence', color = 'orange', hjust = -0.1) +",
    "  labs(title = 'Distribution of Bayes Factors',",
    "       x = 'Bayes Factor (capped at 100)',",
    "       y = 'Frequency') +",
    "  theme_minimal() +",
    "  theme(plot.title = element_text(hjust = 0.5, size = 16),",
    "        axis.title = element_text(size = 14),",
    "        axis.text = element_text(size = 12))",
    "```",
    "",
    "Bayes factors quantify the relative evidence in favor of the alternative hypothesis ",
    "(that there is a difference between groups) compared to the null hypothesis (no difference). ",
    "Values above 3 indicate moderate evidence, while values above 10 indicate strong evidence ",
    "for a true effect.",
    "",
    "### Effect Size and Posterior Probability",
    "",
    "```{r effect-posterior, echo=FALSE, fig.width=8, fig.height=5}",
    "# Extract effect sizes and posterior probabilities from simulations",
    "# Use only a subset of simulations for clarity",
    "n_sims_to_use <- min(5, length(power_data$simulation_results))",
    "effect_sizes <- unlist(lapply(power_data$simulation_results[1:n_sims_to_use], function(x) x$effect_sizes))",
    "post_probs <- unlist(lapply(power_data$simulation_results[1:n_sims_to_use], function(x) x$posterior_probs))",
    "sim_ids <- rep(1:n_sims_to_use, each = length(power_data$simulation_results[[1]]$effect_sizes))",
    "",
    "scatter_df <- data.frame(",
    "  effect_size = effect_sizes,",
    "  posterior_probability = post_probs,",
    "  simulation = factor(sim_ids)",
    ")",
    "",
    "# Create scatter plot",
    "ggplot(scatter_df, aes(x = log2(effect_size), y = posterior_probability, color = simulation)) +",
    "  geom_point(alpha = 0.7, size = 3) +",
    "  geom_hline(yintercept = power_data$parameters$posterior_prob_threshold, ",
    "             linetype = 'dashed', color = 'black') +",
    "  geom_vline(xintercept = log2(power_data$parameters$fold_change_threshold), ",
    "             linetype = 'dashed', color = 'black') +",
    "  geom_vline(xintercept = -log2(power_data$parameters$fold_change_threshold), ",
    "             linetype = 'dashed', color = 'black') +",
    "  labs(title = 'Effect Size vs. Posterior Probability',",
    "       x = 'Log2 Fold Change',",
    "       y = 'Posterior Probability',",
    "       color = 'Simulation') +",
    "  theme_minimal() +",
    "  theme(plot.title = element_text(hjust = 0.5, size = 16),",
    "        axis.title = element_text(size = 14),",
    "        axis.text = element_text(size = 12),",
    "        legend.position = 'bottom')",
    "```",
    "",
    "This scatter plot shows the relationship between effect size (as log2 fold change) and ",
    "posterior probability. Points in the upper right or upper left quadrants represent viral taxa ",
    "with both high posterior probability and substantial effect size, which are more likely ",
    "to be detected as significant.",
    "",
    "## Statistical Performance",
    "",
    "### Detection Rates",
    "",
    "```{r detection-rates, echo=FALSE, fig.width=8, fig.height=5}",
    "# Create bar plot of detection metrics",
    "detection_df <- data.frame(",
    "  Metric = c('True Positives', 'False Positives', 'True Negatives', 'False Negatives'),",
    "  Count = c(",
    "    mean(sapply(power_data$simulation_results, function(x) length(x$true_positives))),",
    "    mean(sapply(power_data$simulation_results, function(x) length(x$false_positives))),",
    "    mean(sapply(power_data$simulation_results, function(x) length(x$true_negatives))),",
    "    mean(sapply(power_data$simulation_results, function(x) length(x$false_negatives)))",
    "  ),",
    "  Type = c('Correct', 'Incorrect', 'Correct', 'Incorrect')",
    ")",
    "",
    "# Create bar plot",
    "ggplot(detection_df, aes(x = Metric, y = Count, fill = Type)) +",
    "  geom_bar(stat = 'identity', color = 'black', alpha = 0.8) +",
    "  geom_text(aes(label = round(Count, 1)), vjust = -0.5, size = 4) +",
    "  labs(title = 'Average Detection Metrics',",
    "       x = '',",
    "       y = 'Average Count') +",
    "  scale_fill_manual(values = c('Correct' = '#4daf4a', 'Incorrect' = '#e41a1c')) +",
    "  theme_minimal() +",
    "  theme(plot.title = element_text(hjust = 0.5, size = 16),",
    "        axis.title = element_text(size = 14),",
    "        axis.text = element_text(size = 12),",
    "        legend.position = 'bottom')",
    "```",
    "",
    "This bar chart shows the average counts of true positives (correctly identified effects), ",
    "false positives (incorrectly identified effects), true negatives (correctly identified non-effects), ",
    "and false negatives (missed effects) across all simulations.",
    "",
    "### Precision and Recall",
    "",
    "```{r precision-recall, echo=FALSE}",
    "# Calculate precision, recall, and F1 score",
    "tp <- mean(sapply(power_data$simulation_results, function(x) length(x$true_positives)))",
    "fp <- mean(sapply(power_data$simulation_results, function(x) length(x$false_positives)))",
    "fn <- mean(sapply(power_data$simulation_results, function(x) length(x$false_negatives)))",
    "",
    "precision <- tp / (tp + fp)",
    "recall <- tp / (tp + fn)",
    "f1_score <- 2 * precision * recall / (precision + recall)",
    "",
    "# Create data frame for precision-recall metrics",
    "pr_df <- data.frame(",
    "  Metric = c('Precision', 'Recall', 'F1 Score'),",
    "  Value = c(precision, recall, f1_score)",
    ")",
    "",
    "# Display the precision-recall table",
    if (has_kableExtra) {
      c(
        "kable(pr_df, col.names = c('Metric', 'Value'),",
        "      digits = 3,",
        "      caption = 'Precision and Recall Metrics') %>%",
        "  kable_styling(bootstrap_options = c('striped', 'hover'), full_width = FALSE,",
        "                position = 'center') %>%",
        "  column_spec(1, bold = TRUE)"
      )
    } else {
      c(
        "kable(pr_df, col.names = c('Metric', 'Value'),",
        "      digits = 3,",
        "      caption = 'Precision and Recall Metrics')"
      )
    },
    "```",
    "",
    paste0("With a Bayesian power of ", round(power_data$power * 100, 1), "%, ",
           "this study design achieves a precision of ", round(precision * 100, 1), "% ",
           "(proportion of detected effects that are true) and a recall of ", round(recall * 100, 1), "% ",
           "(proportion of true effects that are detected). The F1 score, a harmonic mean of precision and recall, ",
           "is ", round(f1_score * 100, 1), "%."),
    "",
    "## Bayesian Method Comparison",
    "",
    "```{r bayesian-comparison, echo=FALSE}",
    "# Create data frame for comparison of Bayesian vs frequentist approaches",
    "comparison_df <- data.frame(",
    "  Feature = c(",
    "    'Interpretation',",
    "    'Prior Knowledge',",
    "    'Multiple Testing',",
    "    'Small Sample Handling',",
    "    'Uncertainty Quantification',",
    "    'Result Format'",
    "  ),",
    "  Bayesian = c(",
    "    'Direct probability of effect existence',",
    "    'Explicitly incorporated',",
    "    'No explicit correction needed',",
    "    'More stable with informative priors',",
    "    'Credible intervals',",
    "    'Posterior probabilities and Bayes factors'",
    "  ),",
    "  Frequentist = c(",
    "    'Evidence against null hypothesis',",
    "    'Not explicitly incorporated',",
    "    'Requires explicit correction',",
    "    'Can be unstable or biased',",
    "    'Confidence intervals',",
    "    'P-values and test statistics'",
    "  )",
    ")",
    "",
    "# Display the comparison table",
    if (has_kableExtra) {
      c(
        "kable(comparison_df, col.names = c('Feature', 'Bayesian Approach', 'Frequentist Approach'),",
        "      caption = 'Comparison of Bayesian and Frequentist Methods') %>%",
        "  kable_styling(bootstrap_options = c('striped', 'hover'), full_width = TRUE,",
        "                position = 'center') %>%",
        "  column_spec(1, bold = TRUE)"
      )
    } else {
      c(
        "kable(comparison_df, col.names = c('Feature', 'Bayesian Approach', 'Frequentist Approach'),",
        "      caption = 'Comparison of Bayesian and Frequentist Methods')"
      )
    },
    "```",
    "",
    "## Recommendations",
    "",
    "Based on the Bayesian power analysis results:",
    "",
    paste0("- The Bayesian power is ", round(power_data$power * 100, 1), "%, ",
           if (power_data$power < 0.8) {
             "which is below the conventionally desired level of 80%. Consider increasing sample size, focusing on stronger effects, or using more informative priors."
           } else {
             "which exceeds the conventionally desired level of 80%. The design has adequate power to detect the specified effect sizes."
           }),
    "",
    paste0("- With the current design, we expect to detect approximately ", 
           round(power_data$expected_discoveries, 1), " viral taxa ",
           "out of ", power_data$parameters$n_viruses, " total taxa."),
    "",
    paste0("- The false discovery proportion is estimated at ", 
           round(power_data$false_discovery_proportion * 100, 1), "%, ",
           if (power_data$false_discovery_proportion > 0.2) {
             "which is relatively high. Consider using a more stringent posterior probability threshold or increasing sample size."
           } else {
             "which is acceptably low, indicating good precision in the results."
           }),
    "",
    paste0("- The prior strength parameter is set to ", power_data$parameters$prior_strength, ". ",
           if (power_data$parameters$prior_strength < 1) {
             "This represents a weakly informative prior. Consider using a stronger prior if you have reliable prior knowledge from previous studies."
           } else if (power_data$parameters$prior_strength < 5) {
             "This represents a moderately informative prior, balancing prior knowledge with the current data."
           } else {
             "This represents a strongly informative prior. Ensure this level of prior knowledge is justified by previous studies."
           }),
    "",
    "## Methods",
    "",
    "This analysis was conducted using the viromePower package, which implements Bayesian power analysis for virome studies. The statistical power was calculated through Monte Carlo simulation, generating realistic virome data and analyzing it with Bayesian methods.",
    "",
    "The key aspects of the Bayesian approach include:",
    "",
    "- A Dirichlet-multinomial model for count data, which naturally handles overdispersion and sparsity",
    "- Prior strength parameter to incorporate existing knowledge about virome abundance patterns",
    "- Calculation of Bayes factors to quantify evidence for differential abundance",
    "- Posterior probability thresholding for decision-making",
    "- Credible intervals to quantify uncertainty in effect size estimates",
    "",
    "Unlike frequentist approaches, the Bayesian framework allows direct statements about the probability of a true effect, incorporates prior knowledge, and does not require explicit adjustment for multiple testing."
  )
  
  # Create the code section if requested
  if (include_code) {
    code_section <- c(
      "## Session Info",
      "",
      "```{r session-info, echo=FALSE}",
      "sessionInfo()",
      "```",
      "",
      "## Reproducibility Code",
      "",
      "```{r reproducible-code, eval=FALSE}",
      "# Code to reproduce this analysis",
      "library(viromePower)",
      "",
      "# Calculate Bayesian power",
      "bayesian_power <- calc_bayesian_power(",
      paste0("  n_samples = ", power_data$parameters$n_samples, ","),
      paste0("  effect_size = ", power_data$parameters$effect_size, ","),
      paste0("  n_viruses = ", power_data$parameters$n_viruses, ","),
      paste0("  sparsity = ", power_data$parameters$sparsity, ","),
      paste0("  dispersion = ", power_data$parameters$dispersion, ","),
      paste0("  prior_strength = ", power_data$parameters$prior_strength, ","),
      paste0("  credible_interval = ", power_data$parameters$credible_interval, ","),
      paste0("  fold_change_threshold = ", power_data$parameters$fold_change_threshold, ","),
      paste0("  posterior_prob_threshold = ", power_data$parameters$posterior_prob_threshold, ","),
      paste0("  n_sim = ", power_data$parameters$n_sim),
      ")",
      "",
      "# Generate report",
      "report <- generate_bayesian_power_report(",
      "  bayesian_power_results = bayesian_power,",
      paste0("  output_file = \"", basename(output_file), "\","),
      paste0("  title = \"", title, "\""),
      ")",
      "```"
    )
  } else {
    code_section <- c(
      "## Session Info",
      "",
      "```{r session-info, echo=FALSE}",
      "sessionInfo()",
      "```"
    )
  }
  
  # Combine all sections
  rmd_content <- c(header, data_setup, intro, code_section)
  
  # Return the complete R markdown content
  return(paste(rmd_content, collapse = "\n"))
}