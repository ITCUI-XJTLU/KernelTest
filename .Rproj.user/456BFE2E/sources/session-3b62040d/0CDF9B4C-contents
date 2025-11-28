# This demo shows how to use the KernelTest package with the updated
# TS_twosample function that includes visualization capabilities.

# Load development version
# Note: After updating the code, make sure to run:
#   1. devtools::document()  - to regenerate documentation
#   2. devtools::load_all()  - to load the updated package
devtools::load_all()

# Load example data
data(data1)  # Condition 1 (case) - 5 genes × 280 bins
data(data4)  # Condition 2 (control) - 5 genes × 280 bins

# Step 1: Estimate bias term (tao) from stable genes
# max1, max4: genes with maximum values below these thresholds are considered "stable"
tao <- est.c(data1, data4, max1 = 5, max4 = 5)
cat("Estimated bias term (tao):", tao, "\n\n")

# Step 2: Run two-sample kernel test with basic visualization
results_basic <- TS_twosample(
  data1 = data1,
  data4 = data4,
  tao = tao,
  band = 180,              # Kernel bandwidth
  quant = c(0.01, 0.01, 0.01),  # Quantile shrinkage for variance estimates
  alpha = 0.05,            # Significance level
  plot = TRUE              # Generate summary plots
)

# Step 3: Access and examine results
cat("Significant genes (indices):", results_basic$significant_genes, "\n")
cat("P-values:", round(results_basic$p_values, 4), "\n")
cat("TS_kn statistics:", round(results_basic$TS_kn, 4), "\n")
cat("Critical value (±):", round(results_basic$z_critical, 4), "\n\n")

# Step 4: Plot specific genes with custom settings
# Example: Plot specific genes with custom names
results_custom <- TS_twosample(
  data1 = data1,
  data4 = data4,
  tao = tao,
  band = 180,
  quant = c(0.01, 0.01, 0.01),
  alpha = 0.05,
  plot = TRUE,
  plot_genes = c(1, 3),    # Plot genes 1 and 3 in detail
  gene_names = c("GENE_A", "GENE_B", "GENE_C", "GENE_D", "GENE_E"),
  condition_names = c("Treatment", "Control")
)

# Step 5: Examine variance and test statistics
cat("Variance (Condition 1):", round(results_basic$sigma1, 4), "\n")
cat("Variance (Condition 2):", round(results_basic$sigma4, 4), "\n")
cat("Equal variance test (Deql):", round(results_basic$Deql, 4), "\n")
cat("Unequal variance test (Dnun):", round(results_basic$Dnun, 4), "\n\n")


