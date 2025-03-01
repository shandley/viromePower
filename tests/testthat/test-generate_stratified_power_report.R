context("Stratified Power Report Generation")

# Skip tests on CRAN to avoid dependencies and file writing issues
test_that("generate_stratified_power_report creates a report file", {
  skip_on_cran()
  
  # Create a temporary file path for the report
  temp_report <- tempfile(fileext = ".html")
  
  # Simulate minimal stratified power results
  mock_strata_power <- list(
    overall_power = 0.475,
    stratum_specific_power = list(stratum_1 = 0.025, stratum_2 = 0.06),
    effective_sample_size = 19.2,
    design_effect_actual = 2.15,
    fdr = 0.23,
    avg_detected = 1.7,
    sim_summary = list(
      true_positives = c(0, 1, 0, 1, 0),
      false_positives = c(1, 0, 2, 0, 1),
      detected_diff = c(1, 1, 2, 1, 1)
    ),
    parameters = list(
      strata_sizes = c(15, 20),
      effect_sizes = c(3.0, 3.5),
      strata_weights = c(0.4, 0.6),
      n_viruses = 20,
      clustering_factor = 0.05,
      stratification_vars = "geography",
      allocation_method = "proportional",
      alpha = 0.05,
      sparsity = 0.6,
      dispersion = 1.5,
      method = "mixed_effects",
      n_sim = 5
    )
  )
  
  # Skip the actual rendering in tests to avoid rmarkdown dependency
  # Just test the function structure and input validation
  
  # Test invalid input
  expect_error(
    generate_stratified_power_report(
      stratified_power_results = list(a = 1, b = 2),
      output_file = temp_report
    ),
    "stratified_power_results must be the output from calc_stratified_power"
  )
  
  # Test valid structure generation (but mock the rendering)
  with_mock(
    `rmarkdown::render` = function(...) TRUE,
    `writeLines` = function(...) TRUE,
    `normalizePath` = function(...) temp_report,
    {
      result <- generate_stratified_power_report(
        stratified_power_results = mock_strata_power,
        output_file = temp_report,
        title = "Test Report",
        include_code = TRUE
      )
      
      # Check result type
      expect_is(result, "character")
      expect_equal(result, temp_report)
    }
  )
})

# Test the Rmd content generation function
test_that("generate_stratified_report_rmd produces valid Rmd content", {
  skip_on_cran()
  
  # Simulate minimal stratified power results
  mock_strata_power <- list(
    overall_power = 0.475,
    stratum_specific_power = list(stratum_1 = 0.025, stratum_2 = 0.06),
    effective_sample_size = 19.2,
    design_effect_actual = 2.15,
    fdr = 0.23,
    avg_detected = 1.7,
    sim_summary = list(
      true_positives = c(0, 1, 0, 1, 0),
      false_positives = c(1, 0, 2, 0, 1),
      detected_diff = c(1, 1, 2, 1, 1)
    ),
    parameters = list(
      strata_sizes = c(15, 20),
      effect_sizes = c(3.0, 3.5),
      strata_weights = c(0.4, 0.6),
      n_viruses = 20,
      clustering_factor = 0.05,
      stratification_vars = "geography",
      allocation_method = "proportional",
      alpha = 0.05,
      sparsity = 0.6,
      dispersion = 1.5,
      method = "mixed_effects",
      n_sim = 5
    )
  )
  
  # Generate Rmd content
  rmd_content <- generate_stratified_report_rmd(
    stratified_power_results = mock_strata_power,
    title = "Test Report",
    include_code = TRUE,
    custom_css = NULL
  )
  
  # Check that expected content is present
  expect_is(rmd_content, "character")
  expect_true(grepl("title: \"Test Report\"", rmd_content))
  expect_true(grepl("library\\(ggplot2\\)", rmd_content))
  expect_true(grepl("Overall statistical power", rmd_content))
  expect_true(grepl("47.5%", rmd_content))
  expect_true(grepl("Power by Stratum", rmd_content))
  expect_true(grepl("Design Effect", rmd_content))
  expect_true(grepl("False Discovery Rate", rmd_content))
  expect_true(grepl("Recommendations", rmd_content))
})