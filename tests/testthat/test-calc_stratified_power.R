context("Stratified Power Analysis")

# Test basic functionality
test_that("calc_stratified_power returns expected structure", {
  # Skip test on CRAN or CI environments to avoid long-running tests
  skip_on_cran()
  
  # Basic test with minimal parameters
  result <- calc_stratified_power(
    strata_sizes = c(5, 7),
    effect_sizes = 2.0,
    n_viruses = 20,
    n_sim = 10  # Use small number for quick testing
  )
  
  # Check structure
  expect_is(result, "list")
  expect_true("overall_power" %in% names(result))
  expect_true("stratum_specific_power" %in% names(result))
  expect_true("effective_sample_size" %in% names(result))
  expect_true("parameters" %in% names(result))
  
  # Check values
  expect_true(result$overall_power >= 0 && result$overall_power <= 1)
  expect_equal(length(result$stratum_specific_power), 2)
})

# Test parameter validation
test_that("calc_stratified_power validates parameters", {
  # Invalid strata sizes
  expect_error(
    calc_stratified_power(
      strata_sizes = c(1, 2),  # Too small
      effect_sizes = 2.0,
      n_viruses = 20
    ),
    "Each stratum must have at least 3 samples per group"
  )
  
  # Invalid effect sizes
  expect_error(
    calc_stratified_power(
      strata_sizes = c(5, 7),
      effect_sizes = c(1.5, 2.0, 2.5),  # Wrong length
      n_viruses = 20
    ),
    "effect_sizes must be either a single value or match the number of strata"
  )
  
  # Invalid strata weights
  expect_error(
    calc_stratified_power(
      strata_sizes = c(5, 7),
      effect_sizes = 2.0,
      strata_weights = c(0.3, 0.4),  # Doesn't sum to 1
      n_viruses = 20
    ),
    "strata_weights must sum to 1"
  )
})

# Test method selection
test_that("calc_stratified_power handles different methods", {
  # Skip test on CRAN or CI environments to avoid long-running tests
  skip_on_cran()
  
  # Try different methods
  methods <- c("stratified_test", "mixed_effects", "survey_adjusted")
  
  for (m in methods) {
    result <- calc_stratified_power(
      strata_sizes = c(5, 7),
      effect_sizes = 2.0,
      n_viruses = 20,
      method = m,
      n_sim = 5  # Very small for quick testing
    )
    
    # Verify the method was used
    expect_equal(result$parameters$method, m)
    expect_true(!is.null(result$overall_power))
  }
})

# Test with real-world parameters similar to README example
test_that("calc_stratified_power works with README example parameters", {
  # Skip test on CRAN or CI environments to avoid long-running tests
  skip_on_cran()
  
  # Test with parameters from README
  result <- calc_stratified_power(
    strata_sizes = c(10, 15),
    effect_sizes = c(1.8, 2.2),
    strata_weights = c(0.4, 0.6),
    n_viruses = 50,
    clustering_factor = 0.1,
    stratification_vars = "geography",
    method = "mixed_effects",
    n_sim = 5  # Very small for quick testing
  )
  
  # Verify the output structure matches what's used in the README
  expect_true("overall_power" %in% names(result))
  expect_true("stratum_specific_power" %in% names(result))
  expect_true("effective_sample_size" %in% names(result))
})