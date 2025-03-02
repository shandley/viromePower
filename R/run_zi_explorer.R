#' Launch the Zero-Inflation Explorer Dashboard
#'
#' Opens an interactive Shiny dashboard for exploring zero-inflation patterns
#' in virome data and their impact on statistical power. This dashboard provides
#' an intuitive interface for simulating, analyzing, and visualizing variable
#' zero-inflation models without requiring programming expertise.
#'
#' @param browser Whether to open the dashboard in a web browser (default: TRUE)
#' @param ... Additional parameters to pass to shiny::runApp()
#'
#' @return Returns the Shiny app object invisibly
#'
#' @details
#' The Zero-Inflation Explorer dashboard includes six primary tabs:
#' \itemize{
#'   \item \strong{Overview:} Introduction to variable zero-inflation concepts
#'   \item \strong{Simulation:} Generate and visualize simulated virome data
#'   \item \strong{Parameter Estimation:} Estimate zero-inflation parameters from real data
#'   \item \strong{Power Analysis:} Compare statistical power across different ZI patterns
#'   \item \strong{Model Comparison:} Statistically compare fixed vs. variable ZI models
#'   \item \strong{Sampling Strategy:} Optimize sample size and sequencing depth
#' }
#'
#' This tool is particularly useful for:
#' \itemize{
#'   \item Researchers designing virome studies who need to determine sample sizes
#'   \item Analysts seeking to understand the impact of zero-inflation on their results
#'   \item Teachers and students learning about zero-inflation modeling
#'   \item Method developers comparing different approaches to handling excess zeros
#' }
#'
#' The dashboard requires several additional R packages beyond viromePower:
#' \itemize{
#'   \item shiny, shinydashboard: For the interactive interface
#'   \item ggplot2, plotly: For visualization
#'   \item DT: For interactive tables
#' }
#' If any required packages are missing, the function will prompt you to install them.
#'
#' @examples
#' \dontrun{
#' # Launch the dashboard with default settings
#' run_zi_explorer()
#'
#' # Launch with specific port and host settings
#' run_zi_explorer(port = 4321, host = "0.0.0.0")
#'
#' # Launch without opening a browser (e.g., for server deployment)
#' run_zi_explorer(browser = FALSE)
#' }
#'
#' @export
run_zi_explorer <- function(browser = TRUE, ...) {
  # Check for required packages
  required_pkgs <- c("shiny", "shinydashboard", "ggplot2", "plotly", "DT")
  missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
  
  if (length(missing_pkgs) > 0) {
    msg <- paste0(
      "The following packages are required but not installed: ",
      paste(missing_pkgs, collapse = ", "), ".\n",
      "Would you like to install them now? (y/n)"
    )
    
    answer <- readline(msg)
    
    if (tolower(answer) %in% c("y", "yes")) {
      utils::install.packages(missing_pkgs)
    } else {
      stop("Cannot launch Zero-Inflation Explorer without required packages.")
    }
  }
  
  # Get the app directory
  app_dir <- system.file("shiny_apps/zi_explorer", package = "viromePower")
  
  if (app_dir == "") {
    stop("Could not find the Zero-Inflation Explorer app. Please reinstall the viromePower package.")
  }
  
  # Launch the app
  shiny::runApp(app_dir, launch.browser = browser, ...)
}