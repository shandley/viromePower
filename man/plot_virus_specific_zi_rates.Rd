\name{plot_virus_specific_zi_rates}
\alias{plot_virus_specific_zi_rates}
\title{Plot Virus-Specific Zero Inflation Rates}
\usage{
plot_virus_specific_zi_rates(
  sim_data,
  binwidth = 0.05,
  include_density = TRUE,
  color_palette = c("navy", "firebrick")
)
}
\arguments{
\item{sim_data}{A list output from simulate_zero_inflated_virome() function}

\item{binwidth}{Width of bins for histogram (default: 0.05)}

\item{include_density}{Whether to overlay density curve (default: TRUE)}

\item{color_palette}{Vector of colors for plot (default: c("navy", "firebrick"))}
}
\value{
A ggplot object showing the distribution of zero-inflation rates
}
\description{
Visualizes the distribution of zero-inflation rates across viral taxa, which is useful
for understanding the heterogeneity of structural zeros in virome data. Different viral
taxa often have different prevalence patterns, leading to varying rates of structural zeros.
}
\examples{
\dontrun{
# Generate simulated data with variable zero-inflation rates
sim_data <- simulate_zero_inflated_virome(
  n_samples = 30, 
  n_viruses = 100,
  variable_zi_rates = TRUE
)

# Plot distribution of zero-inflation rates
plot_virus_specific_zi_rates(sim_data)
}
}