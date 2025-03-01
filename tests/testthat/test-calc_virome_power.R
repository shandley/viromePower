test_that("calc_virome_power returns expected structure", {
  n_samples <- 10
  effect_size <- 1.5
  n_viruses <- 100
  
  power_result <- calc_virome_power(
    n_samples = n_samples, 
    effect_size = effect_size, 
    n_viruses = n_viruses
  )
  
  expect_type(power_result, "list")
  expect_true("power" %in% names(power_result))
  expect_true("parameters" %in% names(power_result))
  
  expect_type(power_result$power, "double")
  expect_true(power_result$power >= 0 && power_result$power <= 1)
  
  expect_equal(power_result$parameters$n_samples, n_samples)
  expect_equal(power_result$parameters$effect_size, effect_size)
  expect_equal(power_result$parameters$n_viruses, n_viruses)
})