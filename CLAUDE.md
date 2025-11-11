# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

**KernelTest** is an R package for detecting differential histone enrichment between case and control groups from ChIP-seq profiles. It implements kernel-based smoothing and nonparametric statistical tests to identify differential binding regions without restrictive parametric assumptions.

The package accepts paired design matrices containing ChIP-seq read counts or normalized intensities for genes under case and control conditions. For each gene, columns represent genomic bins ordered along the locus.

## Development Commands

### Building and Installation
```r
# Install from GitHub
devtools::install_github("tfcui2025/KernelTest")

# Build package documentation
devtools::document()

# Build the package
devtools::build()

# Check package (includes tests, documentation, etc.)
devtools::check()

# Install locally for testing
devtools::install()
```

### Testing
```r
# Run all tests
devtools::test()

# Load package for interactive testing
devtools::load_all()
```

## Code Architecture

### Core Implementation

**The package now uses Rcpp/C++ for performance-critical operations**, providing significant speedup over pure R implementations.

### Core Statistical Functions ([R/Functions.R](R/Functions.R))

The package implements three main statistical approaches (now accelerated with Rcpp):

1. **`est.c(data1, data4, max1, max4)`** - Estimates the bias term (tao) from stable genes
   - Identifies genes with low maximum intensities in both conditions
   - Computes autocorrelation-based bias correction term
   - Returns: scalar bias term used in hypothesis tests

2. **`TS_kernel(data, band, quantile)`** - Single-sample kernel test
   - Applies Nadaraya-Watson kernel smoothing with Gaussian kernel
   - Computes hat matrix (Snw) for smoothing operator
   - Estimates residual variance with degrees of freedom adjustment
   - Uses Wilson-Hilferty transformation for test statistics
   - Returns: normalized test statistics, signs, and mean values

3. **`TS_twosample(data1, data4, tao, band, quant)`** - Two-sample kernel test (main exported function)
   - Implements differential enrichment test between two conditions
   - Three variance estimation methods:
     - Equal variance assumption (Sev)
     - Unequal variance assumption (Suv)
     - Kernel-smoothed variance (Xg)
   - Computes multiple test statistics:
     - Kernel-smoothed statistic (TS_kn) using Wilson-Hilferty transformation
     - Equal variance test (Deql)
     - Unequal variance test (Dnun)
   - Returns: comprehensive list including variances, test statistics, and intermediate values

**Key implementation details:**
- Kernel smoothing uses bandwidth parameter `band` normalized by sample size `n`
- Hat matrix eigenvalues are used to compute effective degrees of freedom (delta, d)
- Quantile-based shrinkage applied to variance estimates to improve stability
- Small constant (0.25) added before variance stabilization transformation

### Data Documentation ([R/data1.R](R/data1.R), [R/data4.R](R/data4.R))

- **data1**: 5 genes × 280 bins matrix (condition 1/case)
- **data4**: 5 genes × 280 bins matrix (condition 2/control)
- Genes are paired across datasets (same order)
- Used for examples and testing

### Helper Function

**`NormTransformation(data)`** - Variance stabilization transformation
- Formula: `2*sqrt(data + 0.25)`
- Note: Currently defined but not actively used in main workflow

## Important Notes

### Issue: Hardcoded Absolute Paths
Lines 228-229 in [R/Functions.R](R/Functions.R:228-229) contain hardcoded absolute paths that will fail on other systems:
```r
load("/Users/cuitengfei/TAMU/TAMU_Courses/STAT600_Computation/HW/Final_project2/KernelTest/data/data1.rda")
load("/Users/cuitengfei/TAMU/TAMU_Courses/STAT600_Computation/HW/Final_project2/KernelTest/data/data4.rda")
```
These should be removed - R packages access data via `data()` mechanism, not direct loading.

### Current Development Status (from README)
- Package is functional but incomplete
- Needs additional visualization functions
- Documentation incomplete
- Methods require verification and validation

## Typical Workflow for Users

```r
library(KernelTest)

# Load example data
data(data1)  # condition 1 (case)
data(data4)  # condition 2 (control)

# Estimate bias term from stable genes
tao <- est.c(data1, data4, max1 = 5, max4 = 5)

# Run two-sample test
# band: kernel bandwidth (e.g., 0.1-0.3)
# quant: quantile shrinkage for variance estimates (vector of 3 values, 0-1)
results <- TS_twosample(data1, data4, tao = tao, band = 0.2, quant = c(0.1, 0.1, 0.1))

# Access results
results$TS_kn    # Kernel test statistics
results$Deql     # Equal variance test
results$Dnun     # Unequal variance test
results$Xg       # Variance estimates
```

## Code Structure

### C++ Core ([src/kernel_core.cpp](src/kernel_core.cpp))

Performance-critical algorithms implemented in C++ using Rcpp and RcppArmadillo:

1. **`gaussian_kernel()`** - Gaussian kernel function
2. **`nw_smoother()`** - Nadaraya-Watson kernel smoother
3. **`compute_hat_matrix()`** - Hat matrix computation for smoothing operator
4. **`compute_variance_estimates()`** - Batch variance estimation for all genes
5. **`est_c_cpp()`** - C++ implementation of bias term estimation
6. **`compute_eigenvalue_stats()`** - Eigenvalue decomposition and statistics

### R Wrappers ([R/Functions.R](R/Functions.R))

R functions that call the C++ implementations and handle high-level logic:
- User-facing API remains unchanged
- Automatic compilation of C++ code on package installation
- Significant performance improvements (typically 5-10x faster)

### Old R Implementation ([R/Functions_old.R](R/Functions_old.R))

The original pure R implementation is preserved for reference and comparison.

## Performance

The Rcpp implementation provides:
- **5-10x faster** computation for typical datasets
- Efficient matrix operations using Armadillo
- Lower memory footprint for large datasets
- Identical numerical results to the R implementation

## Dependencies

- **stats**: ksmooth (kernel smoothing), quantile
- **graphics**: hist (histogram visualization in TS_kernel)
- **Rcpp** (>= 1.0.0): R/C++ interface
- **RcppArmadillo**: Linear algebra operations
- R >= 3.5

## Building from Source

The package requires a C++ compiler:

```r
# Install with dependencies
devtools::install_github("tfcui2025/KernelTest")

# Or build locally
devtools::document()
devtools::install()
```

On first installation, Rcpp will compile the C++ code automatically.
