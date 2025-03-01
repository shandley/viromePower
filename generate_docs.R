# Script to generate documentation for viromePower package

# Check for required packages
if (!requireNamespace("roxygen2", quietly = TRUE)) {
  cat("Installing roxygen2 package...\n")
  install.packages("roxygen2", repos = "https://cloud.r-project.org")
}

# Load roxygen2
library(roxygen2)

# Set working directory to package root
# setwd("/Users/scott/Handley Lab Dropbox/Scott Handley/restructure/code/R/viromePower")

# Create man directory if it doesn't exist
if (!dir.exists("man")) {
  dir.create("man")
  cat("Created man directory\n")
}

# Generate documentation
cat("Generating documentation...\n")
roxygen2::roxygenize(".")

cat("\nDocumentation generation complete. Man pages should be available in the 'man' directory.\n")