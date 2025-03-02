\name{plot_zi_by_abundance}
\alias{plot_zi_by_abundance}
\title{Plot Zero Inflation by Viral Abundance}
\usage{
plot_zi_by_abundance(
  sim_data,
  abundance_scale = "log",
  add_regression = TRUE,
  color_palette = c("navy", "firebrick")
)
}
\arguments{
\item{sim_data}{A list output from simulate_zero_inflated_virome() function}

\item{abundance_scale}{Scale for abundance values ("log" or "linear", default: "log")}

\item{add_regression}{Whether to add regression line (default: TRUE)}

\item{color_palette}{Vector of colors for plot (default: c("navy", "firebrick"))}
}
\value{
A ggplot object showing relationship between abundance and zero-inflation rates
}
\description{
Visualizes the relationship between viral abundance and zero-inflation rates. This helps
understand whether rare viral taxa have higher rates of structural zeros, which is a common
pattern in virome data.
}
\examples{
\dontrun{
# Generate simulated data with variable zero-inflation rates
sim_data <- simulate_zero_inflated_virome(
  n_samples = 30, 
  n_viruses = 100,
  variable_zi_rates = TRUE
)

# Plot relationship between abundance and zero-inflation rates
plot_zi_by_abundance(sim_data)
}
}