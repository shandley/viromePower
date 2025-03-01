# viromePower Development Guidelines

## Build/Test Commands
- Install dependencies: `R -e "install.packages(c('devtools', 'roxygen2', 'testthat', 'knitr'))"` 
- Build package: `R CMD build .`
- Check package: `R CMD check --as-cran viromePower_*.tar.gz`
- Run tests: `R -e "devtools::test()"`
- Run single test: `R -e "devtools::test_file('tests/testthat/test_file_name.R')"`
- Document: `R -e "devtools::document()"`
- Install local package: `R -e "devtools::install()"`

## Code Style Guidelines
- Follow tidyverse style guide for R (https://style.tidyverse.org/)
- Functions: Use snake_case for function names
- Variables: Use snake_case for variables
- Use roxygen2 for documentation with @param, @return, @examples tags
- Include tests for all user-facing functions in tests/testthat/
- Keep functions focused on a single task with clear inputs/outputs
- Handle errors with appropriate error messages using stop() or warning()
- Maintain backward compatibility when updating functions
- Document dependencies in DESCRIPTION file