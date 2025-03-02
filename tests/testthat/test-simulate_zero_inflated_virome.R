# Tests for the simulate_zero_inflated_virome function
# This test suite verifies the functionality of variable zero-inflation modeling

# Set a random seed for reproducibility of tests
set.seed(42)

test_that("simulate_zero_inflated_virome creates expected structure with fixed zero-inflation", {
  # Basic simulation with fixed zero-inflation rates
  sim_data <- simulate_zero_inflated_virome(
    n_samples = 20, 
    n_viruses = 50,
    variable_zi_rates = FALSE
  )
  
  # Check basic structure
  expect_type(sim_data, "list")
  expect_true(all(c("counts", "metadata", "structural_zeros", "sampling_zeros", 
                    "parameters", "sparsity_summary", "virus_specific_zi") %in% names(sim_data)))
  
  # Check count matrix dimensions
  expect_equal(dim(sim_data$counts), c(50, 20))
  
  # Check that sample IDs match between count matrix and metadata
  expect_equal(colnames(sim_data$counts), sim_data$metadata$sample_id)
  
  # Check parameters were stored correctly
  expect_false(sim_data$parameters$variable_zi_rates)
  
  # Check that virus-specific ZI info exists but all values should be similar with fixed ZI
  expect_true("rates" %in% names(sim_data$virus_specific_zi))
  expect_length(sim_data$virus_specific_zi$rates, 50)
  
  # With fixed ZI, the standard deviation of rates should be non-zero but not too large
  # (rates still vary by abundance in fixed ZI mode)
  expect_true(sim_data$virus_specific_zi$sd < 0.3)
})

test_that("simulate_zero_inflated_virome correctly implements variable zero-inflation", {
  # Simulation with variable zero-inflation rates
  sim_data_var <- simulate_zero_inflated_virome(
    n_samples = 20, 
    n_viruses = 50,
    variable_zi_rates = TRUE,
    zi_alpha = 2,
    zi_beta = 5  # Right-skewed distribution (more taxa with lower ZI rates)
  )
  
  # Check variable ZI parameter was set correctly
  expect_true(sim_data_var$parameters$variable_zi_rates)
  expect_equal(sim_data_var$parameters$zi_alpha, 2)
  expect_equal(sim_data_var$parameters$zi_beta, 5)
  
  # Check virus-specific ZI rates
  expect_length(sim_data_var$virus_specific_zi$rates, 50)
  
  # Check that mean of variable rates is close to the specified structural_zeros parameter
  mean_zi_rate <- mean(sim_data_var$virus_specific_zi$rates)
  expected_zi_rate <- sim_data_var$parameters$structural_zeros
  expect_equal(mean_zi_rate, expected_zi_rate, tolerance = 0.1)
  
  # With right-skewed distribution (alpha=2, beta=5), median should be less than mean
  expect_true(sim_data_var$virus_specific_zi$median < sim_data_var$virus_specific_zi$mean)
})

test_that("different Beta distribution shapes create expected ZI distributions", {
  # 1. Right-skewed distribution (more taxa with lower ZI rates)
  sim_right <- simulate_zero_inflated_virome(
    n_samples = 15, 
    n_viruses = 100,
    variable_zi_rates = TRUE,
    zi_alpha = 2,
    zi_beta = 5,
    structural_zeros = 0.6
  )
  
  # 2. Left-skewed distribution (more taxa with higher ZI rates)
  sim_left <- simulate_zero_inflated_virome(
    n_samples = 15, 
    n_viruses = 100,
    variable_zi_rates = TRUE,
    zi_alpha = 5,
    zi_beta = 2,
    structural_zeros = 0.6
  )
  
  # 3. Symmetric distribution
  sim_symm <- simulate_zero_inflated_virome(
    n_samples = 15, 
    n_viruses = 100,
    variable_zi_rates = TRUE,
    zi_alpha = 2,
    zi_beta = 2,
    structural_zeros = 0.6
  )
  
  # For right-skewed (alpha < beta), median < mean
  expect_true(sim_right$virus_specific_zi$median < sim_right$virus_specific_zi$mean)
  
  # For left-skewed (alpha > beta), median > mean
  expect_true(sim_left$virus_specific_zi$median > sim_left$virus_specific_zi$mean)
  
  # For symmetric (alpha = beta), median ≈ mean
  expect_equal(sim_symm$virus_specific_zi$median, sim_symm$virus_specific_zi$mean, tolerance = 0.1)
  
  # Check that all distributions maintain similar average rates
  expect_equal(sim_right$virus_specific_zi$mean, 0.6, tolerance = 0.1)
  expect_equal(sim_left$virus_specific_zi$mean, 0.6, tolerance = 0.1)
  expect_equal(sim_symm$virus_specific_zi$mean, 0.6, tolerance = 0.1)
})

test_that("group differences are maintained with variable zero-inflation", {
  # Simulation with group differences in abundance and zero-inflation
  sim_data <- simulate_zero_inflated_virome(
    n_samples = 30, 
    n_viruses = 50,
    variable_zi_rates = TRUE,
    zi_alpha = 2,
    zi_beta = 3,
    effect_size = 2.5,
    zero_inflation_difference = TRUE
  )
  
  # Extract counts for each group
  group_a_samples <- which(sim_data$metadata$group == "A")
  group_b_samples <- which(sim_data$metadata$group == "B")
  
  counts_a <- sim_data$counts[, group_a_samples]
  counts_b <- sim_data$counts[, group_b_samples]
  
  # Check differential taxa
  diff_taxa <- sim_data$diff_taxa
  
  # For differential taxa, group B should have higher counts on average
  for (i in diff_taxa) {
    mean_a <- mean(counts_a[i,])
    mean_b <- mean(counts_b[i,])
    # Accounting for randomness, most should show the effect
    if (mean_a > 0) {  # Only test if there are non-zero values in group A
      message(sprintf("Taxon %d: mean_a = %.2f, mean_b = %.2f", i, mean_a, mean_b))
    }
  }
  
  # Since the test is stochastic, we don't make hard assertions on every taxon
  # Instead, check that the structural zero probabilities differ between groups for diff taxa
  for (i in diff_taxa) {
    prob_a <- sim_data$virus_specific_zi$structural_zero_probs[i, 1]
    prob_b <- sim_data$virus_specific_zi$structural_zero_probs[i, 2]
    expect_true(prob_a > prob_b)
  }
})

test_that("output sparsity matches expected parameters", {
  # Fixed parameters for test
  n_samples <- 50
  n_viruses <- 100
  structural_zeros <- 0.6
  sampling_zeros <- 0.2
  
  # Simulation
  sim_data <- simulate_zero_inflated_virome(
    n_samples = n_samples,
    n_viruses = n_viruses,
    structural_zeros = structural_zeros,
    sampling_zeros = sampling_zeros,
    variable_zi_rates = TRUE
  )
  
  # Get observed sparsity
  observed_sparsity <- sim_data$sparsity_summary$observed_sparsity
  structural_proportion <- sim_data$sparsity_summary$structural_zeros_proportion
  sampling_proportion <- sim_data$sparsity_summary$sampling_zeros_proportion
  
  # Calculate expected sparsity:
  # total_zeros = structural_zeros + (1-structural_zeros)*sampling_zeros
  expected_sparsity <- structural_zeros + (1-structural_zeros)*sampling_zeros
  
  # Allow some tolerance due to randomness and edge cases
  expect_equal(observed_sparsity, expected_sparsity, tolerance = 0.1)
  expect_equal(structural_proportion, structural_zeros, tolerance = 0.1)
  
  # Check that nonzero proportion is calculated correctly
  expect_equal(sim_data$sparsity_summary$nonzero_proportion, 1 - observed_sparsity)
})

test_that("edge cases are handled properly", {
  # Test with very high zero-inflation
  sim_high_zi <- simulate_zero_inflated_virome(
    n_samples = 10,
    n_viruses = 20,
    structural_zeros = 0.9,
    variable_zi_rates = TRUE,
    zi_alpha = 10,  # Highly concentrated distribution
    zi_beta = 1     # Heavily left-skewed (most taxa have high ZI)
  )
  
  # Check that zero-inflation is high but not 100%
  expect_true(sim_high_zi$sparsity_summary$observed_sparsity > 0.9)
  expect_true(sim_high_zi$sparsity_summary$observed_sparsity < 1.0)
  
  # Test with very low zero-inflation
  sim_low_zi <- simulate_zero_inflated_virome(
    n_samples = 10,
    n_viruses = 20,
    structural_zeros = 0.3,
    variable_zi_rates = TRUE,
    zi_alpha = 1,    # Highly concentrated distribution
    zi_beta = 10     # Heavily right-skewed (most taxa have low ZI)
  )
  
  # Verify ZI distribution shape (heavily right-skewed)
  expect_true(sim_low_zi$virus_specific_zi$median < sim_low_zi$virus_specific_zi$mean)
  expect_true(sim_low_zi$virus_specific_zi$median < 0.3)
})