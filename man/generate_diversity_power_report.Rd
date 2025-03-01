% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_diversity_power_report.R
\name{generate_diversity_power_report}
\alias{generate_diversity_power_report}
\title{Generate Viral Diversity Power Analysis Report}
\usage{
generate_diversity_power_report(
  n_samples,
  effect_size,
  n_viruses,
  diversity_measure = "shannon",
  alpha = 0.05,
  sparsity = 0.8,
  dispersion = 2,
  n_sim = 100,
  output_file = "diversity_power_report.html",
  show_individual_sims = FALSE
)
}
\arguments{
\item{n_samples}{Number of samples per group to analyze}

\item{effect_size}{Expected effect size between groups}

\item{n_viruses}{Number of viral taxa in the dataset}

\item{diversity_measure}{Type of diversity to analyze: "shannon", "simpson", "richness", "evenness", "chao1", "ace", "inv_simpson", "fisher_alpha", "goods_coverage", "berger_parker", "bray", "jaccard", "sorensen", "morisita_horn", "ruzicka", "cca" (default: "shannon")}

\item{alpha}{Significance level (default: 0.05)}

\item{sparsity}{Proportion of zeros in the data (default: 0.8)}

\item{dispersion}{Dispersion parameter for viral abundance (default: 2)}

\item{n_sim}{Number of simulations for power calculation (default: 100)}

\item{output_file}{Path to save the HTML report (default: "diversity_power_report.html")}

\item{show_individual_sims}{Logical indicating whether to show individual simulation points (default: FALSE)}
}
\value{
Path to the generated HTML report
}
\description{
Generates a comprehensive HTML report for viral diversity power analysis,
including sample size recommendations, power curves, and visualizations.
}
\examples{
# Generate a power report for Shannon diversity
report_file <- generate_diversity_power_report(
  n_samples = 15, 
  effect_size = 1.2,
  n_viruses = 200,
  diversity_measure = "shannon"
)

# Generate a power report for Inverse Simpson index
inv_simpson_report <- generate_diversity_power_report(
  n_samples = 15,
  effect_size = 1.0,
  n_viruses = 200,
  diversity_measure = "inv_simpson",
  output_file = "inv_simpson_power_report.html"
)

# Generate a power report for ACE richness estimator
ace_report <- generate_diversity_power_report(
  n_samples = 18,
  effect_size = 1.5,
  n_viruses = 250,
  diversity_measure = "ace",
  output_file = "ace_power_report.html"
)

# Generate a power report for Good's coverage
goods_report <- generate_diversity_power_report(
  n_samples = 20,
  effect_size = 1.2,
  n_viruses = 200,
  diversity_measure = "goods_coverage",
  output_file = "goods_coverage_report.html"
)

# Generate a power report for Berger-Parker dominance index
bp_report <- generate_diversity_power_report(
  n_samples = 15,
  effect_size = 1.0,
  n_viruses = 150,
  diversity_measure = "berger_parker",
  output_file = "berger_parker_report.html"
)

# Generate a power report for Bray-Curtis dissimilarity
beta_report <- generate_diversity_power_report(
  n_samples = 20,
  effect_size = 0.15,
  n_viruses = 300,
  diversity_measure = "bray",
  output_file = "beta_diversity_power_report.html"
)
}
