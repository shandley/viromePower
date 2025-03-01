test_that("calc_viral_diversity_power works with Shannon diversity", {
  result <- calc_viral_diversity_power(
    n_samples = 5,
    effect_size = 1.5,
    n_viruses = 50,
    diversity_measure = "shannon",
    n_sim = 5  # Small number for testing speed
  )
  
  expect_type(result, "list")
  expect_true("power" %in% names(result))
  expect_true("parameters" %in% names(result))
  expect_true("significant_tests" %in% names(result))
  expect_true(result$power >= 0 && result$power <= 1)
})

test_that("calc_viral_diversity_power works with Simpson diversity", {
  result <- calc_viral_diversity_power(
    n_samples = 5,
    effect_size = 1.5,
    n_viruses = 50,
    diversity_measure = "simpson",
    n_sim = 5  # Small number for testing speed
  )
  
  expect_type(result, "list")
  expect_true("power" %in% names(result))
  expect_true(result$power >= 0 && result$power <= 1)
})

test_that("calc_viral_diversity_power works with Inverse Simpson diversity", {
  result <- calc_viral_diversity_power(
    n_samples = 5,
    effect_size = 1.5,
    n_viruses = 50,
    diversity_measure = "inv_simpson",
    n_sim = 5  # Small number for testing speed
  )
  
  expect_type(result, "list")
  expect_true("power" %in% names(result))
  expect_true(result$power >= 0 && result$power <= 1)
})

test_that("calc_viral_diversity_power works with Chao1 richness estimator", {
  result <- calc_viral_diversity_power(
    n_samples = 5,
    effect_size = 1.5,
    n_viruses = 50,
    diversity_measure = "chao1",
    n_sim = 5  # Small number for testing speed
  )
  
  expect_type(result, "list")
  expect_true("power" %in% names(result))
  expect_true(result$power >= 0 && result$power <= 1)
})

test_that("calc_viral_diversity_power works with ACE richness estimator", {
  result <- calc_viral_diversity_power(
    n_samples = 5,
    effect_size = 1.5,
    n_viruses = 50,
    diversity_measure = "ace",
    n_sim = 5  # Small number for testing speed
  )
  
  expect_type(result, "list")
  expect_true("power" %in% names(result))
  expect_true(result$power >= 0 && result$power <= 1)
})


test_that("calc_viral_diversity_power works with Fisher's alpha", {
  result <- calc_viral_diversity_power(
    n_samples = 5,
    effect_size = 1.5,
    n_viruses = 50,
    diversity_measure = "fisher_alpha",
    n_sim = 5  # Small number for testing speed
  )
  
  expect_type(result, "list")
  expect_true("power" %in% names(result))
  expect_true(result$power >= 0 && result$power <= 1)
})

test_that("calc_viral_diversity_power works with Bray-Curtis (beta diversity)", {
  result <- calc_viral_diversity_power(
    n_samples = 5,
    effect_size = 0.2,
    n_viruses = 50,
    diversity_measure = "bray",
    n_sim = 5  # Small number for testing speed
  )
  
  expect_type(result, "list")
  expect_true("power" %in% names(result))
  expect_true(result$power >= 0 && result$power <= 1)
})

test_that("calc_viral_diversity_power validates input parameters", {
  # Invalid n_samples
  expect_error(
    calc_viral_diversity_power(
      n_samples = 2,  # Too small
      effect_size = 1.5,
      n_viruses = 50,
      diversity_measure = "shannon"
    ),
    "n_samples must be at least 3 per group"
  )
  
  # Invalid effect_size
  expect_error(
    calc_viral_diversity_power(
      n_samples = 5,
      effect_size = -0.5,  # Negative
      n_viruses = 50,
      diversity_measure = "shannon"
    ),
    "effect_size must be positive"
  )
  
  # Invalid diversity_measure
  expect_error(
    calc_viral_diversity_power(
      n_samples = 5,
      effect_size = 1.5,
      n_viruses = 50,
      diversity_measure = "invalid"  # Not supported
    ),
    "diversity_measure must be one of"
  )
})

test_that("plot_diversity_power_curve returns a ggplot object", {
  # Test with n_samples parameter type
  curve1 <- plot_diversity_power_curve(
    param_range = c(5, 10),
    param_type = "n_samples",
    diversity_measure = "shannon",
    n_sim = 5  # Small number for testing speed
  )
  
  expect_s3_class(curve1, "ggplot")
  
  # Test with effect_size parameter type
  curve2 <- plot_diversity_power_curve(
    param_range = c(0.5, 1.0),
    param_type = "effect_size",
    diversity_measure = "bray",
    n_sim = 5  # Small number for testing speed
  )
  
  expect_s3_class(curve2, "ggplot")
})

test_that("generate_diversity_power_report creates an HTML file", {
  # Skip on CRAN or CI environments
  skip_on_cran()
  skip_if_not_installed("rmarkdown")
  
  # Create a temporary file for the report
  temp_file <- tempfile(fileext = ".html")
  
  # Generate report with minimal parameters
  result <- generate_diversity_power_report(
    n_samples = 5,
    effect_size = 1.0,
    n_viruses = 30,
    diversity_measure = "shannon",
    output_file = temp_file,
    n_sim = 5  # Small number for testing speed
  )
  
  # Check that the file exists and has content
  expect_true(file.exists(temp_file))
  expect_gt(file.size(temp_file), 1000)  # Should be substantial
  
  # Clean up
  if (file.exists(temp_file)) {
    file.remove(temp_file)
  }
})