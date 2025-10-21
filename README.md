
# KernelTest

<!-- badges: start -->
<!-- badges: end -->


`KernelTest` will accept paired design matrices containing ChIP-seq read counts or normalized intensities for genes under case and control conditions. For each gene, columns will represent genomic bins ordered along the locus. The package will implement three kernel smoothing procedures—Gaussian, Epanechnikov, and adaptive bandwidth variants—to estimate smoothed enrichment profiles and conduct hypothesis tests for differential binding. Key outputs will include gene-level test statistics, variance estimates, and \(p\)-values. Diagnostic plots will visualize observed and smoothed signals, highlighting regions classified as significantly enriched. The package will expose a streamlined workflow with documented functions, reproducible examples, and vignettes.

- Deliver an open-source R package that implements modern kernel-based differential histone enrichment tests.
- Extend prior methodology with more robust kernel smoothers tailored to heterogeneous ChIP-seq signals.
- Provide clear documentation, reproducible examples, and visualization tools to facilitate biological interpretation.


## Installation

You can install the development version from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("tfcui2025/KernelTest")
```

## Next Work: 

- Add some visualization functions to show results
- Finish all documents
- Check the methods and corresponding functions in the R scripts



