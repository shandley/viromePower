#' Generate Stratified Power Analysis Report
#'
#' Creates an HTML report with visualizations and summaries of stratified power analysis results.
#' This report includes comparisons of power across strata, effective sample size calculations,
#' and design effect visualizations for complex sampling designs.
#'
#' @param stratified_power_results List of results from calc_stratified_power()
#' @param output_file File path for the output HTML report
#' @param title Title for the report (default: "Stratified Virome Power Analysis Report")
#' @param include_code Boolean indicating whether to include R code in the report (default: FALSE)
#' @param custom_css Optional custom CSS for the report styling
#'
#' @return The path to the generated HTML report
#' @export
#'
#' @examples
#' # Calculate stratified power
#' power_result <- calc_stratified_power(
#'   strata_sizes = c(15, 20),
#'   effect_sizes = c(3.0, 3.5),
#'   strata_weights = c(0.4, 0.6),
#'   n_viruses = 20,
#'   clustering_factor = 0.05,
#'   sparsity = 0.6,
#'   dispersion = 1.5,
#'   method = "mixed_effects"
#' )
#'
#' # Generate report
#' report_path <- generate_stratified_power_report(
#'   stratified_power_results = power_result,
#'   output_file = "stratified_power_report.html",
#'   title = "Stratified Virome Power Analysis"
#' )
generate_stratified_power_report <- function(stratified_power_results,
                                            output_file = "stratified_power_report.html",
                                            title = "Stratified Virome Power Analysis Report",
                                            include_code = FALSE,
                                            custom_css = NULL) {
  
  # Validate input
  if (!is.list(stratified_power_results) || 
      !all(c("overall_power", "stratum_specific_power", "parameters") %in% names(stratified_power_results))) {
    stop("stratified_power_results must be the output from calc_stratified_power()")
  }
  
  # Check for required packages
  required_packages <- c("rmarkdown", "ggplot2", "dplyr", "tidyr", "knitr", "kableExtra")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    stop(paste0(
      "Missing required packages for report generation: ", 
      paste(missing_packages, collapse = ", "), 
      ". Please install them with: install.packages(c('", 
      paste(missing_packages, collapse = "', '"), 
      "'))"
    ))
  }
  
  # Create a temporary Rmd file with the report content
  rmd_file <- tempfile(fileext = ".Rmd")
  
  # Generate the R markdown content
  rmd_content <- generate_stratified_report_rmd(
    stratified_power_results = stratified_power_results,
    title = title,
    include_code = include_code,
    custom_css = custom_css,
    output_file = output_file
  )
  
  # Write the content to the temporary file
  writeLines(rmd_content, rmd_file)
  
  # Render the report using rmarkdown
  tryCatch({
    rmarkdown::render(
      input = rmd_file,
      output_file = basename(output_file),
      output_dir = dirname(output_file),
      quiet = TRUE
    )
  }, error = function(e) {
    stop(paste0("Error rendering report: ", e$message))
  })
  
  # Return the path to the generated report
  return(normalizePath(output_file))
}

#' Generate R Markdown content for stratified power report
#'
#' Helper function to create the R markdown content for the stratified power report
#'
#' @param stratified_power_results List of results from calc_stratified_power()
#' @param title Title for the report
#' @param include_code Boolean indicating whether to include R code
#' @param custom_css Optional custom CSS for styling
#' @param output_file File path for the output HTML report
#'
#' @return A character string with the R markdown content
#' @keywords internal
generate_stratified_report_rmd <- function(stratified_power_results,
                                          title,
                                          include_code,
                                          custom_css,
                                          output_file) {
  
  # Extract data for easier use in the report
  power_data <- stratified_power_results
  
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
    "library(dplyr)",
    "library(tidyr)",
    "library(knitr)",
    "library(kableExtra)",
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
    if (!is.null(custom_css)) custom_css else "",
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
    paste0("<p>Overall statistical power: <span class=\"power-value\">", 
           round(power_data$overall_power * 100, 1), "%</span></p>"),
    paste0("<p>Effective sample size: <span class=\"power-value\">", 
           round(power_data$effective_sample_size, 1), 
           "</span> out of <span class=\"power-value\">", 
           sum(power_data$parameters$strata_sizes) * 2, "</span> total samples</p>"),
    paste0("<p>Design effect: <span class=\"power-value\">", 
           round(power_data$design_effect_actual, 2), "</span></p>"),
    "</div>",
    "",
    "This report presents a comprehensive analysis of statistical power for a stratified virome study design. ",
    "The stratified approach accounts for the complex structure of the population and provides more accurate ",
    "power estimates than a simple design.",
    "",
    "### Study Design",
    "",
    "```{r study-design, echo=FALSE}",
    "# Create a data frame for the stratum characteristics",
    "strata_df <- data.frame(",
    "  Stratum = seq_along(power_data$parameters$strata_sizes),",
    "  Sample_Size = power_data$parameters$strata_sizes,",
    "  Effect_Size = power_data$parameters$effect_sizes,",
    "  Weight = power_data$parameters$strata_weights,",
    "  Sparsity = if(length(power_data$parameters$sparsity) == 1) {",
    "    rep(power_data$parameters$sparsity, length(power_data$parameters$strata_sizes))",
    "  } else {",
    "    power_data$parameters$sparsity",
    "  },",
    "  Dispersion = if(length(power_data$parameters$dispersion) == 1) {",
    "    rep(power_data$parameters$dispersion, length(power_data$parameters$strata_sizes))",
    "  } else {",
    "    power_data$parameters$dispersion",
    "  },",
    "  Power = sapply(power_data$stratum_specific_power, function(x) x * 100)",
    ")",
    "",
    "# Display the stratum characteristics table",
    "kable(strata_df, col.names = c('Stratum', 'Sample Size (per group)', ",
    "                              'Effect Size', 'Weight', 'Sparsity', ",
    "                              'Dispersion', 'Power (%)'),",
    "      digits = c(0, 0, 2, 2, 2, 2, 1),",
    "      caption = 'Characteristics and Power by Stratum') %>%",
    "  kable_styling(bootstrap_options = 'striped', full_width = FALSE,",
    "                position = 'center', row_label_position = 'c') %>%",
    "  column_spec(7, background = '#e6f3ff')",
    "```",
    "",
    "The study design incorporates ", length(power_data$parameters$strata_sizes), 
    " strata with different characteristics as shown above. ",
    "The analysis was conducted using the '", power_data$parameters$method, "' method ",
    "with a significance level (α) of ", power_data$parameters$alpha, ".",
    "",
    "## Power Analysis Results",
    "",
    "### Power by Stratum",
    "",
    "```{r power-by-stratum, echo=FALSE, fig.width=8, fig.height=5}",
    "# Create a bar plot of power by stratum",
    "ggplot(strata_df, aes(x = factor(Stratum), y = Power, fill = factor(Stratum))) +",
    "  geom_bar(stat = 'identity', color = 'black', alpha = 0.8) +",
    "  geom_text(aes(label = paste0(round(Power, 1), '%')), vjust = -0.5, size = 4) +",
    "  geom_hline(yintercept = power_data$overall_power * 100, linetype = 'dashed', color = 'red', size = 1) +",
    "  annotate('text', x = length(power_data$parameters$strata_sizes) / 2 + 0.5, ",
    "           y = power_data$overall_power * 100 + 5, ",
    "           label = paste0('Overall: ', round(power_data$overall_power * 100, 1), '%'), ",
    "           color = 'red', size = 4.5) +",
    "  labs(title = 'Statistical Power by Stratum',",
    "       x = 'Stratum',",
    "       y = 'Power (%)') +",
    "  theme_minimal() +",
    "  scale_fill_brewer(palette = 'Blues') +",
    "  theme(legend.position = 'none',",
    "        plot.title = element_text(hjust = 0.5, size = 16),",
    "        axis.title = element_text(size = 14),",
    "        axis.text = element_text(size = 12))",
    "```",
    "",
    "The bar plot shows the statistical power for each stratum, with the overall power indicated by the red dashed line. ",
    "Power variations across strata reflect differences in sample size, effect size, and data characteristics.",
    "",
    "### Effect Size vs. Power",
    "",
    "```{r effect-size-power, echo=FALSE, fig.width=8, fig.height=5}",
    "# Create a scatter plot of effect size vs. power",
    "ggplot(strata_df, aes(x = Effect_Size, y = Power, size = Sample_Size, color = factor(Stratum))) +",
    "  geom_point(alpha = 0.8) +",
    "  geom_text(aes(label = paste0('Stratum ', Stratum)), hjust = -0.3, vjust = 1.5, size = 4) +",
    "  labs(title = 'Effect Size vs. Power by Stratum',",
    "       x = 'Effect Size',",
    "       y = 'Power (%)',",
    "       size = 'Sample Size',",
    "       color = 'Stratum') +",
    "  theme_minimal() +",
    "  scale_size_continuous(range = c(5, 15)) +",
    "  theme(plot.title = element_text(hjust = 0.5, size = 16),",
    "        axis.title = element_text(size = 14),",
    "        axis.text = element_text(size = 12),",
    "        legend.position = 'right')",
    "```",
    "",
    "This scatter plot illustrates the relationship between effect size and power across strata. ",
    "The size of each point corresponds to the sample size in that stratum. ",
    "As expected, strata with larger effect sizes and sample sizes generally show higher power.",
    "",
    "### Sample Size and Stratification Impact",
    "",
    "```{r design-effect, echo=FALSE, fig.width=8, fig.height=5}",
    "# Create a visualization of the design effect",
    "design_df <- data.frame(",
    "  Design = c('Simple Random Sampling', 'Stratified Design'),",
    "  Sample_Size = c(sum(power_data$parameters$strata_sizes) * 2, power_data$effective_sample_size),",
    "  Type = c('Actual', 'Effective')",
    ")",
    "",
    "ggplot(design_df, aes(x = Design, y = Sample_Size, fill = Type)) +",
    "  geom_bar(stat = 'identity', position = 'dodge', color = 'black', alpha = 0.8) +",
    "  geom_text(aes(label = round(Sample_Size, 1)), position = position_dodge(width = 0.9),",
    "            vjust = -0.5, size = 4) +",
    "  labs(title = 'Design Effect: Actual vs. Effective Sample Size',",
    "       x = '',",
    "       y = 'Sample Size') +",
    "  theme_minimal() +",
    "  scale_fill_brewer(palette = 'Set2') +",
    "  theme(plot.title = element_text(hjust = 0.5, size = 16),",
    "        axis.title = element_text(size = 14),",
    "        axis.text = element_text(size = 12))",
    "```",
    "",
    paste0("The design effect (", round(power_data$design_effect_actual, 2), ") quantifies how the stratified sampling design affects the effective sample size. ", 
           "A value greater than 1 indicates that the stratification and clustering reduce the effective sample size ", 
           "compared to a simple random sample of the same total size."),
    "",
    "## Study Parameters",
    "",
    "```{r parameters-table, echo=FALSE}",
    "# Create a data frame for the parameters",
    "params <- power_data$parameters",
    "param_df <- data.frame(",
    "  Parameter = c(",
    "    'Number of Viruses',",
    "    'Clustering Factor',",
    "    'Stratification Variable',",
    "    'Allocation Method',",
    "    'Alpha (Significance Level)',",
    "    'Method',",
    "    'Number of Simulations'",
    "  ),",
    "  Value = c(",
    "    params$n_viruses,",
    "    params$clustering_factor,",
    "    params$stratification_vars,",
    "    params$allocation_method,",
    "    params$alpha,",
    "    params$method,",
    "    params$n_sim",
    "  )",
    ")",
    "",
    "# Display the parameters table",
    "kable(param_df, col.names = c('Parameter', 'Value'),",
    "      caption = 'Analysis Parameters') %>%",
    "  kable_styling(bootstrap_options = c('striped', 'hover'), full_width = FALSE,",
    "                position = 'center') %>%",
    "  column_spec(1, bold = TRUE)",
    "```",
    "",
    "## False Discovery Rate Analysis",
    "",
    "```{r fdr-analysis, echo=FALSE, fig.width=8, fig.height=5}",
    "# Create a bar plot of true/false positives",
    "discovery_df <- data.frame(",
    "  Type = c('True Positives', 'False Positives'),",
    "  Count = c(mean(power_data$sim_summary$true_positives), ",
    "            mean(power_data$sim_summary$false_positives))",
    ")",
    "",
    "ggplot(discovery_df, aes(x = Type, y = Count, fill = Type)) +",
    "  geom_bar(stat = 'identity', color = 'black', alpha = 0.8) +",
    "  geom_text(aes(label = round(Count, 2)), vjust = -0.5, size = 4) +",
    "  labs(title = paste0('Discovery Analysis (FDR: ', round(power_data$fdr * 100, 1), '%)'),",
    "       x = '',",
    "       y = 'Average Count') +",
    "  theme_minimal() +",
    "  scale_fill_manual(values = c('True Positives' = '#4daf4a', 'False Positives' = '#e41a1c')) +",
    "  theme(plot.title = element_text(hjust = 0.5, size = 16),",
    "        axis.title = element_text(size = 14),",
    "        axis.text = element_text(size = 12),",
    "        legend.position = 'none')",
    "```",
    "",
    paste0("The multiple testing correction controls the False Discovery Rate (FDR) at approximately ", 
           round(power_data$fdr * 100, 1), "%. ", 
           "On average, the analysis detects ", round(power_data$avg_detected, 2), 
           " differentially abundant viral taxa, with an expected ", 
           round(mean(power_data$sim_summary$true_positives), 2), 
           " true positives and ", 
           round(mean(power_data$sim_summary$false_positives), 2), 
           " false positives per study."),
    "",
    "## Recommendations",
    "",
    "Based on the power analysis results for this stratified design:",
    "",
    paste0("- The overall power is ", round(power_data$overall_power * 100, 1), "%, ",
           if (power_data$overall_power < 0.8) {
             "which is below the conventionally desired level of 80%. Consider increasing sample sizes or focusing on strata with higher effect sizes."
           } else {
             "which exceeds the conventionally desired level of 80%. The design has adequate power to detect the specified effect sizes."
           }),
    "",
    "- Stratum-specific powers vary, with ",
    paste0(paste0("Stratum ", which.max(sapply(power_data$stratum_specific_power, function(x) x)), 
                 " showing the highest power (", 
                 round(max(sapply(power_data$stratum_specific_power, function(x) x)) * 100, 1), "%)"), 
           "."),
    "",
    paste0("- The design effect of ", round(power_data$design_effect_actual, 2), 
           " suggests that the stratification ", 
           if (power_data$design_effect_actual > 1.5) {
             "substantially reduces the effective sample size. Consider reducing the clustering or allocating more samples to key strata."
           } else if (power_data$design_effect_actual > 1.2) {
             "moderately reduces the effective sample size. The benefits of stratification should be weighed against this reduction."
           } else {
             "has minimal impact on the effective sample size, indicating an efficient stratified design."
           }),
    "",
    "## Methods",
    "",
    "This analysis was conducted using the viromePower package, which simulates realistic virome data with appropriate sparsity and dispersion characteristics. The statistical power was calculated through Monte Carlo simulation with multiple testing correction to control the false discovery rate.",
    "",
    paste0("For the '", power_data$parameters$method, "' method, the analysis accounts for the stratified design by ", 
           if (power_data$parameters$method == "mixed_effects") {
             "incorporating random effects for strata and modeling the hierarchical structure of the data."
           } else if (power_data$parameters$method == "stratified_test") {
             "conducting separate tests within each stratum and combining the results using stratum weights."
           } else {
             "adjusting the variance estimates to account for the design complexity."
           }),
    "",
    "The calculations incorporate:",
    "- Multiple testing correction using Benjamini-Hochberg procedure",
    "- Design effect adjustments for clustered sampling",
    "- Stratum-specific weights for importance sampling",
    "- Realistic simulation of virome count data with appropriate sparsity and dispersion"
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
      "# Calculate stratified power",
      "stratified_power <- calc_stratified_power(",
      paste0("  strata_sizes = c(", paste(power_data$parameters$strata_sizes, collapse = ", "), "),"),
      if (length(power_data$parameters$effect_sizes) > 1) {
        paste0("  effect_sizes = c(", paste(power_data$parameters$effect_sizes, collapse = ", "), "),")
      } else {
        paste0("  effect_sizes = ", power_data$parameters$effect_sizes, ",")
      },
      if (!is.null(power_data$parameters$strata_weights)) {
        paste0("  strata_weights = c(", paste(power_data$parameters$strata_weights, collapse = ", "), "),")
      },
      paste0("  n_viruses = ", power_data$parameters$n_viruses, ","),
      paste0("  clustering_factor = ", power_data$parameters$clustering_factor, ","),
      if (!is.null(power_data$parameters$stratification_vars)) {
        paste0("  stratification_vars = \"", power_data$parameters$stratification_vars, "\",")
      },
      if (!is.null(power_data$parameters$allocation_method)) {
        paste0("  allocation_method = \"", power_data$parameters$allocation_method, "\",")
      },
      paste0("  alpha = ", power_data$parameters$alpha, ","),
      if (length(power_data$parameters$sparsity) > 1) {
        paste0("  sparsity = c(", paste(power_data$parameters$sparsity, collapse = ", "), "),")
      } else {
        paste0("  sparsity = ", power_data$parameters$sparsity, ",")
      },
      if (length(power_data$parameters$dispersion) > 1) {
        paste0("  dispersion = c(", paste(power_data$parameters$dispersion, collapse = ", "), "),")
      } else {
        paste0("  dispersion = ", power_data$parameters$dispersion, ",")
      },
      paste0("  method = \"", power_data$parameters$method, "\","),
      paste0("  n_sim = ", power_data$parameters$n_sim),
      ")",
      "",
      "# Generate report",
      "report <- generate_stratified_power_report(",
      "  stratified_power_results = stratified_power,",
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
  rmd_content <- c(header, intro, code_section)
  
  # Return the complete R markdown content
  return(paste(rmd_content, collapse = "\n"))
}