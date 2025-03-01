# Test script for viral diversity power analysis functions

# Just source the R files directly without using devtools
source("R/simulate_virome_data.R")
source("R/calc_viral_diversity_power.R")
source("R/plot_diversity_power_curve.R")
source("R/generate_diversity_power_report.R")

# Test the alpha diversity power calculation
# Try with different diversity measures
cat("Testing Shannon diversity power calculation...\n")
power_shannon <- calc_viral_diversity_power(
  n_samples = 20, 
  effect_size = 1.5, 
  n_viruses = 150,
  diversity_measure = "shannon",
  sparsity = 0.75,
  dispersion = 1.8,
  n_sim = 30  # Using fewer simulations for quick testing
)
print(power_shannon$power)

cat("\nTesting Good's Coverage power calculation...\n")
power_coverage <- calc_viral_diversity_power(
  n_samples = 20, 
  effect_size = 2.5, 
  n_viruses = 150,
  diversity_measure = "goods_coverage",
  sparsity = 0.7,  # Less sparse = more difference in coverage
  n_sim = 30
)
print(power_coverage$power)

cat("\nTesting Berger-Parker dominance index power calculation...\n")
power_bp <- calc_viral_diversity_power(
  n_samples = 20, 
  effect_size = 2.0, 
  n_viruses = 100,
  diversity_measure = "berger_parker",
  sparsity = 0.7,
  n_sim = 30
)
print(power_bp$power)

cat("\nTesting Species Richness power calculation...\n")
power_richness <- calc_viral_diversity_power(
  n_samples = 10, 
  effect_size = 1.5, 
  n_viruses = 100,
  diversity_measure = "richness",
  n_sim = 30
)
print(power_richness$power)

# Test beta diversity power calculation
cat("\nTesting Bray-Curtis dissimilarity power calculation...\n")
power_bray <- calc_viral_diversity_power(
  n_samples = 18, 
  effect_size = 0.25, 
  n_viruses = 200,
  diversity_measure = "bray",
  sparsity = 0.75,
  dispersion = 2.0,
  n_sim = 30
)
print(power_bray$power)

cat("\nTesting Morisita-Horn similarity power calculation...\n")
power_morisita <- calc_viral_diversity_power(
  n_samples = 15, 
  effect_size = 0.3, 
  n_viruses = 150,
  diversity_measure = "morisita_horn",
  sparsity = 0.7,
  dispersion = 1.8,
  n_sim = 30
)
print(power_morisita$power)

cat("\nTesting Ružička distance power calculation...\n")
power_ruzicka <- calc_viral_diversity_power(
  n_samples = 20, 
  effect_size = 0.2, 
  n_viruses = 180,
  diversity_measure = "ruzicka",
  sparsity = 0.75,
  n_sim = 30
)
print(power_ruzicka$power)

# Test power curve plotting
cat("\nGenerating power curve for sample size...\n")
sample_curve <- plot_diversity_power_curve(
  param_range = seq(5, 20, by = 5),
  diversity_measure = "shannon",
  n_sim = 20
)
print("Power curve generated")

# Skip report generation as it requires pandoc
cat("\nSkipping report generation (requires pandoc)...\n")

cat("\nAll tests completed successfully!\n")