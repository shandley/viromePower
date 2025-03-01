test_that("simulate_virome_data creates data with expected dimensions", {
  n_samples <- 20
  n_viruses <- 100
  
  sim_data <- simulate_virome_data(n_samples = n_samples, n_viruses = n_viruses)
  
  expect_type(sim_data, "list")
  expect_named(sim_data, c("counts", "metadata"))
  
  expect_equal(dim(sim_data$counts), c(n_viruses, n_samples))
  expect_equal(nrow(sim_data$metadata), n_samples)
  expect_true("sample_id" %in% colnames(sim_data$metadata))
  expect_true("group" %in% colnames(sim_data$metadata))
})