\name{compare_group_zi_patterns}
\alias{compare_group_zi_patterns}
\title{Compare Group-Specific Zero Inflation Patterns}
\usage{
compare_group_zi_patterns(
  sim_data,
  plot_type = "boxplot",
  color_palette = c("steelblue", "firebrick")
)
}
\arguments{
\item{sim_data}{A list output from simulate_zero_inflated_virome() function}

\item{plot_type}{Type of comparison to plot: "boxplot" or "scatter" (default: "boxplot")}

\item{color_palette}{Vector of colors for plot (default: c("steelblue", "firebrick"))}
}
\value{
A ggplot object showing comparison of zero-inflation rates between groups
}
\description{
Visualizes and compares zero-inflation patterns between experimental groups, which is 
useful for understanding differential detection rates across groups.
}
\examples{
\dontrun{
# Generate simulated data with variable zero-inflation rates and group differences
sim_data <- simulate_zero_inflated_virome(
  n_samples = 30, 
  n_viruses = 100,
  variable_zi_rates = TRUE,
  zero_inflation_difference = TRUE
)

# Compare zero-inflation patterns between groups
compare_group_zi_patterns(sim_data)
}
}