# Check for implementation messages
# Source all R files
r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
for (file in r_files) {
  cat("Checking file:", file, "\n")
  lines <- readLines(file)
  message_lines <- grep("Add actual implementation", lines)
  if (length(message_lines) > 0) {
    cat("  Found 'Add actual implementation' message at line(s):", 
        paste(message_lines, collapse=", "), "\n")
    for (line in message_lines) {
      cat("  ", lines[line], "\n")
    }
  } else {
    cat("  No implementation messages found.\n")
  }
}