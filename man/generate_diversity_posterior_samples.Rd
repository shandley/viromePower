% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_bayesian_power.R
\name{generate_diversity_posterior_samples}
\alias{generate_diversity_posterior_samples}
\title{Generate Posterior Samples for Diversity Metrics}
\usage{
generate_diversity_posterior_samples(
  counts_by_group,
  metric,
  prior_strength = 2,
  n_samples = 1000
)
}
\arguments{
\item{counts_by_group}{List of count matrices split by group}

\item{metric}{Diversity metric to analyze}

\item{prior_strength}{Prior strength parameter}

\item{n_samples}{Number of posterior samples to generate}
}
\value{
List with posterior samples for each group and their difference
}
\description{
Internal function to generate MCMC samples for Bayesian analysis of diversity metrics
}
\keyword{internal}
