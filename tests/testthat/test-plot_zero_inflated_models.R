context("Zero-Inflated Model Visualization Functions")

test_that("plot_observed_vs_expected_zeros returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  
  # Simulate zero-inflated data
  n_samples <- 20
  n_viruses <- 50
  
  # Create mock data matrix
  set.seed(123)
  obs_data <- matrix(rnbinom(n_samples * n_viruses, mu = 10, size = 0.5), 
                     nrow = n_samples, ncol = n_viruses)
  # Add some zeros to simulate zero-inflation
  zero_indices <- sample(1:(n_samples * n_viruses), 
                         size = round(0.3 * n_samples * n_viruses))
  obs_data[zero_indices] <- 0
  
  # Define model parameters
  standard_model_params <- list(mu = 10, size = 0.5)
  zinb_model_params <- list(mu = 10, size = 0.5, zi_prob = 0.3)
  
  # Test without group factor
  plot_result <- plot_observed_vs_expected_zeros(
    obs_data = obs_data,
    standard_model_params = standard_model_params,
    zinb_model_params = zinb_model_params
  )
  
  # Check if gridExtra is installed and check result type accordingly
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    expect_s3_class(plot_result, "arrangelist")
  } else {
    expect_s3_class(plot_result, "ggplot")
  }
  
  # Test with group factor
  group_factor <- rep(c("A", "B"), each = n_samples/2)
  plot_result_with_groups <- plot_observed_vs_expected_zeros(
    obs_data = obs_data,
    standard_model_params = standard_model_params,
    zinb_model_params = zinb_model_params,
    group_factor = group_factor
  )
  
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    expect_s3_class(plot_result_with_groups, "arrangelist")
  } else {
    expect_s3_class(plot_result_with_groups, "ggplot")
  }
  
  # Test with custom title
  plot_result_custom_title <- plot_observed_vs_expected_zeros(
    obs_data = obs_data,
    standard_model_params = standard_model_params,
    zinb_model_params = zinb_model_params,
    title = "Custom Title"
  )
  
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    expect_s3_class(plot_result_custom_title, "arrangelist")
  } else {
    expect_s3_class(plot_result_custom_title, "ggplot")
  }
})

test_that("plot_zero_inflation_distribution returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  
  # Create test data
  n_samples <- 30
  n_viruses <- 60
  
  # Create mock ZINB model fit (matrix of zero-inflation probabilities)
  set.seed(456)
  zinb_model_fit <- list(
    pi = matrix(rbeta(n_samples * n_viruses, 2, 5), nrow = n_samples),
    taxa_names = paste0("Taxon_", 1:n_viruses)
  )
  
  # Test basic functionality
  plot_result <- plot_zero_inflation_distribution(
    zinb_model_fit = zinb_model_fit,
    by_taxa = FALSE
  )
  
  expect_s3_class(plot_result, "ggplot")
  
  # Test with by_taxa = TRUE
  plot_result_by_taxa <- plot_zero_inflation_distribution(
    zinb_model_fit = zinb_model_fit,
    by_taxa = TRUE
  )
  
  expect_s3_class(plot_result_by_taxa, "ggplot")
  
  # Test with group factor
  group_factor <- rep(c("A", "B"), each = n_samples/2)
  plot_result_with_groups <- plot_zero_inflation_distribution(
    zinb_model_fit = zinb_model_fit,
    by_taxa = FALSE,
    group_factor = group_factor
  )
  
  expect_s3_class(plot_result_with_groups, "ggplot")
  
  # Test with highlight threshold
  plot_result_with_threshold <- plot_zero_inflation_distribution(
    zinb_model_fit = zinb_model_fit,
    by_taxa = FALSE,
    highlight_threshold = 0.5
  )
  
  expect_s3_class(plot_result_with_threshold, "ggplot")
})

test_that("compare_power_curves returns a ggplot object", {
  skip_if_not_installed("ggplot2")
  
  # Create mock power results
  n_points <- 10
  sample_sizes <- seq(10, 100, length.out = n_points)
  std_power <- 0.3 + 0.5 * (1:n_points) / n_points
  zinb_power <- 0.2 + 0.7 * (1:n_points) / n_points
  
  # Standard model results
  standard_power_results <- list(
    n_samples = sample_sizes,
    power = std_power,
    effect_size = rep(2.0, n_points)
  )
  
  # Zero-inflated model results
  zinb_power_results <- list(
    n_samples = sample_sizes,
    power = zinb_power,
    effect_size = rep(2.0, n_points)
  )
  
  # Test with n_samples as x-variable
  plot_result <- compare_power_curves(
    standard_power_results = standard_power_results,
    zinb_power_results = zinb_power_results,
    x_variable = "n_samples"
  )
  
  expect_s3_class(plot_result, "ggplot")
  
  # Test with highlight threshold
  plot_result_with_thresh <- compare_power_curves(
    standard_power_results = standard_power_results,
    zinb_power_results = zinb_power_results,
    x_variable = "n_samples",
    highlight_thresh = 0.8
  )
  
  expect_s3_class(plot_result_with_thresh, "ggplot")
  
  # Data frame input format
  std_df <- data.frame(
    n_samples = sample_sizes,
    power = std_power,
    effect_size = rep(2.0, n_points)
  )
  
  zinb_df <- data.frame(
    n_samples = sample_sizes,
    power = zinb_power,
    effect_size = rep(2.0, n_points)
  )
  
  # Test with data frame inputs
  plot_result_df <- compare_power_curves(
    standard_power_results = std_df,
    zinb_power_results = zinb_df,
    x_variable = "n_samples"
  )
  
  expect_s3_class(plot_result_df, "ggplot")
  
  # Test with effect_size as x-variable
  effect_sizes <- seq(1.2, 3.0, length.out = n_points)
  std_effect_results <- list(
    effect_size = effect_sizes,
    power = std_power,
    n_samples = rep(50, n_points)
  )
  
  zinb_effect_results <- list(
    effect_size = effect_sizes,
    power = zinb_power,
    n_samples = rep(50, n_points)
  )
  
  plot_result_effect <- compare_power_curves(
    standard_power_results = std_effect_results,
    zinb_power_results = zinb_effect_results,
    x_variable = "effect_size"
  )
  
  expect_s3_class(plot_result_effect, "ggplot")
})

test_that("plot_zinb_diagnostics returns a list of ggplot objects", {
  skip_if_not_installed("ggplot2")
  
  # Create test data
  n_samples <- 25
  n_viruses <- 40
  
  # Create mock data and model fit
  set.seed(789)
  obs_data <- matrix(rnbinom(n_samples * n_viruses, mu = 10, size = 0.5), 
                     nrow = n_samples, ncol = n_viruses)
  # Add some zeros
  zero_indices <- sample(1:(n_samples * n_viruses), 
                        size = round(0.3 * n_samples * n_viruses))
  obs_data[zero_indices] <- 0
  
  # Create mock ZINB model fit
  zinb_model_fit <- list(
    mu = matrix(rgamma(n_samples * n_viruses, 5, 0.5), nrow = n_samples),
    size = 0.5,
    zi_prob = matrix(rbeta(n_samples * n_viruses, 2, 5), nrow = n_samples)
  )
  
  # Run diagnostics with default parameters
  diagnostics <- plot_zinb_diagnostics(
    obs_data = obs_data,
    zinb_model_fit = zinb_model_fit
  )
  
  # Check returned object type and structure
  expect_type(diagnostics, "list")
  expect_named(diagnostics, c("zero_proportion_plot", "qq_plot", "taxa_fit_plots"))
  
  # Check that each element is a ggplot
  expect_s3_class(diagnostics$zero_proportion_plot, "ggplot")
  expect_s3_class(diagnostics$qq_plot, "ggplot")
  expect_s3_class(diagnostics$taxa_fit_plots, "ggplot")
  
  # Test with custom parameters
  custom_diagnostics <- plot_zinb_diagnostics(
    obs_data = obs_data,
    zinb_model_fit = zinb_model_fit,
    n_taxa_to_plot = 6,
    seed = 123
  )
  
  expect_s3_class(custom_diagnostics$zero_proportion_plot, "ggplot")
  expect_s3_class(custom_diagnostics$qq_plot, "ggplot")
  expect_s3_class(custom_diagnostics$taxa_fit_plots, "ggplot")
})

test_that("zero-inflated visualization functions handle edge cases", {
  skip_if_not_installed("ggplot2")
  
  # Small dataset
  n_samples <- 5
  n_viruses <- 10
  
  # Generate minimal test data
  set.seed(101)
  obs_data <- matrix(rnbinom(n_samples * n_viruses, mu = 10, size = 0.5), 
                    nrow = n_samples, ncol = n_viruses)
  obs_data[sample(1:(n_samples * n_viruses), size = 15)] <- 0
  
  # Test plot_observed_vs_expected_zeros with minimal data
  standard_model_params <- list(mu = 10, size = 0.5)
  zinb_model_params <- list(mu = 10, size = 0.5, zi_prob = 0.3)
  
  small_plot <- plot_observed_vs_expected_zeros(
    obs_data = obs_data,
    standard_model_params = standard_model_params,
    zinb_model_params = zinb_model_params
  )
  
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    expect_s3_class(small_plot, "arrangelist")
  } else {
    expect_s3_class(small_plot, "ggplot")
  }
  
  # Test plot_zero_inflation_distribution with minimal data
  small_zinb_fit <- list(
    pi = matrix(rbeta(n_samples * n_viruses, 2, 5), nrow = n_samples),
    taxa_names = paste0("Taxon_", 1:n_viruses)
  )
  
  small_zi_plot <- plot_zero_inflation_distribution(
    zinb_model_fit = small_zinb_fit,
    by_taxa = FALSE
  )
  
  expect_s3_class(small_zi_plot, "ggplot")
  
  # Test compare_power_curves with minimal data
  n_points <- 3
  sample_sizes <- c(10, 20, 30)
  std_power <- c(0.3, 0.5, 0.7)
  zinb_power <- c(0.4, 0.6, 0.8)
  
  small_std_results <- list(
    n_samples = sample_sizes,
    power = std_power,
    effect_size = rep(2.0, n_points)
  )
  
  small_zinb_results <- list(
    n_samples = sample_sizes,
    power = zinb_power,
    effect_size = rep(2.0, n_points)
  )
  
  small_compare_plot <- compare_power_curves(
    standard_power_results = small_std_results,
    zinb_power_results = small_zinb_results,
    x_variable = "n_samples"
  )
  
  expect_s3_class(small_compare_plot, "ggplot")
  
  # Test plot_zinb_diagnostics with minimal data
  small_model_fit <- list(
    mu = matrix(rgamma(n_samples * n_viruses, 5, 0.5), nrow = n_samples),
    size = 0.5,
    zi_prob = matrix(rbeta(n_samples * n_viruses, 2, 5), nrow = n_samples)
  )
  
  small_diagnostics <- plot_zinb_diagnostics(
    obs_data = obs_data,
    zinb_model_fit = small_model_fit,
    n_taxa_to_plot = 2
  )
  
  expect_s3_class(small_diagnostics$zero_proportion_plot, "ggplot")
  expect_s3_class(small_diagnostics$qq_plot, "ggplot")
  expect_s3_class(small_diagnostics$taxa_fit_plots, "ggplot")
})