# Fix end of line issues in R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in r_files) {
  cat("Fixing file:", file, "\n")
  lines <- readLines(file)
  # Write back with proper newline at end
  writeLines(lines, file)
  cat("  Done.\n")
}