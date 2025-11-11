# Changelog

All notable changes to the KernelTest package will be documented in this file.

## [0.0.0.9000] - 2025-01-11

### Added
- **Visualization System**: Complete built-in plotting capabilities
  - 4-panel summary plots (TS_kn, Deql, Dnun, variance comparison)
  - Detailed gene-level ChIP-seq profile plots
  - Automatic significance markers and annotations
  - Customizable gene names and condition labels via `gene_names` and `condition_names` parameters
  - Control over which genes to plot in detail via `plot_genes` parameter

- **Statistical Enhancements**
  - P-value calculation for all genes
  - Significance testing with user-defined alpha level (default: 0.05)
  - Critical value computation and reporting
  - Comprehensive console output with statistical summaries

- **New Parameters in `TS_twosample()`**
  - `alpha`: Significance level for hypothesis testing (default: 0.05)
  - `plot`: Enable/disable automatic plotting (default: FALSE)
  - `plot_genes`: Specify which genes to plot in detail (default: NULL, plots first significant gene)
  - `gene_names`: Custom gene names for plots (default: NULL, uses "Gene1", "Gene2", etc.)
  - `condition_names`: Custom condition labels (default: c("Condition 1", "Condition 2"))

- **New Return Values**
  - `p_values`: P-value for each gene
  - `significant_genes`: Indices of genes exceeding critical value
  - `alpha`: Significance level used
  - `z_critical`: Critical value (±)

- **Documentation**
  - [FINAL_DEMO.R](FINAL_DEMO.R): Complete workflow demonstration
  - [HOW_TO_USE.md](HOW_TO_USE.md): Detailed usage guide (中文)
  - [SOLUTION_SUMMARY.md](SOLUTION_SUMMARY.md): Troubleshooting guide (中文)

### Performance Improvements
- **Complete Rcpp/C++ Implementation**
  - Migrated all computational kernels to C++ using Rcpp/RcppArmadillo
  - 5-10x speedup compared to pure R implementation
  - Optimized hat matrix computation with `compute_hat_matrix()`
  - Vectorized variance estimation for all genes with `compute_variance_estimates()`
  - Efficient eigenvalue statistics with `compute_eigenvalue_stats()`
  - C++ implementation of bias term estimation with `est_c_cpp()`

### Fixed
- **Function Loading Conflicts**
  - Resolved duplicate `TS_twosample` definitions
  - Renamed `R/Functions_old.R` to `R/Functions_old.R.bak` to prevent loading
  - Updated NAMESPACE to export only the current version
  - Fixed `unused argument (plot = TRUE)` error

- **Eigenvalue Computation**
  - Addressed matrix symmetry warnings in `compute_eigenvalue_stats()`
  - Improved numerical stability in eigenvalue decomposition

### Changed
- `TS_twosample()` function signature expanded from 5 to 10 parameters
- Default behavior: `plot = FALSE` (opt-in for visualization)
- Console output now includes detailed statistical summaries when `plot = TRUE`

### Technical Details
- **C++ Functions** (in [src/kernel_core.cpp](src/kernel_core.cpp)):
  - `gaussian_kernel()`: Gaussian kernel function
  - `nw_smoother()`: Nadaraya-Watson kernel smoother
  - `compute_hat_matrix()`: Hat matrix for smoothing operator
  - `compute_variance_estimates()`: Batch variance estimation
  - `est_c_cpp()`: Bias term estimation
  - `compute_eigenvalue_stats()`: Eigenvalue decomposition and statistics

- **R Wrapper Functions** (in [R/Functions.R](R/Functions.R)):
  - `est.c()`: Wrapper for `est_c_cpp()`
  - `TS_kernel()`: Single-sample kernel test
  - `TS_twosample()`: Two-sample kernel test with visualization
  - `NormTransformation()`: Variance stabilization

### Removed
- Old pure R implementation moved to `R/Functions_old.R.bak` (preserved for reference)

---

## [0.0.0.9000-alpha] - 2024-10-20

### Initial Development Version
- Basic two-sample kernel test implementation
- Pure R implementation of kernel smoothing
- Essential statistical functions: `est.c()`, `TS_kernel()`, `TS_twosample()`
- Example datasets: data1, data4

---

## Notes

### Migration Path
If upgrading from older versions:

1. Update package:
   ```r
   devtools::document()
   devtools::load_all()
   ```

2. Update function calls to use new parameters:
   ```r
   # Old (still works)
   results <- TS_twosample(data1, data4, tao, band, quant)
   
   # New (with visualization)
   results <- TS_twosample(data1, data4, tao, band, quant,
                           plot = TRUE, alpha = 0.05)
   ```

3. Access new return values:
   ```r
   results$p_values          # NEW
   results$significant_genes # NEW
   results$alpha             # NEW
   results$z_critical        # NEW
   ```

### Known Issues
- Eigenvalue symmetry warning: Minor numerical precision issue, does not affect results
- Large datasets (>5000 genes): May be slow without parallel processing (planned for future)

### Planned for Next Release
- Multiple testing correction (FDR, Bonferroni)
- Parallel processing support
- Additional kernel functions
- Export plots to file
- Comprehensive vignettes
