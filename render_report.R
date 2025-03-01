# Install required packages if not already installed
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  install.packages("rmarkdown", repos = "https://cloud.r-project.org")
}

if (!requireNamespace("knitr", quietly = TRUE)) {
  install.packages("knitr", repos = "https://cloud.r-project.org")
}

# Render the report
rmarkdown::render(
  input = "report_template.Rmd",
  output_file = "sample_report.html",
  quiet = TRUE
)

# Check if the report was created
if (file.exists("sample_report.html")) {
  cat("Report file created: sample_report.html\n")
} else {
  cat("Report file was not created\n")
}