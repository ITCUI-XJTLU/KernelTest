#' Estimate bias term using Rcpp (wrapper for est_c_cpp)
#'
#' @param data1 First condition data matrix
#' @param data4 Second condition data matrix
#' @param max1 Maximum threshold for condition 1
#' @param max4 Maximum threshold for condition 2
#' @return Bias term (tao)
#' @export
est.c <- function(data1, data4, max1 = 5, max4 = 5) {
  # Use Rcpp version
  return(est_c_cpp(data1, data4, max1, max4))
}

#' Normalize transformation
#' @param data Input data
#' @return Normalized data
#' @export
NormTransformation <- function(data) {
  data.normal <- 2 * sqrt(data + 0.25)
  return(data.normal)
}

#' Single-sample kernel test using Rcpp acceleration
#'
#' @param data Data matrix (genes x bins)
#' @param band Kernel bandwidth
#' @param quantile Quantile for variance shrinkage
#' @return List containing test statistics, signs, and mean values
#' @export
TS_kernel <- function(data, band, quantile) {
  n <- ncol(data)
  hwidth <- band / n

  # Compute hat matrix using Rcpp
  Snw <- compute_hat_matrix(n, hwidth)
  ad_df <- sum(diag(Snw))

  # Initialize result vectors
  n_genes <- nrow(data)
  Xg <- numeric(n_genes)
  x <- seq(1, n) / n

  # Compute variance for each gene
  for (k in 1:n_genes) {
    diff41 <- data[k, ]
    d41_y <- nw_smoother(x, diff41, x, hwidth)
    sigma <- sqrt(sum((d41_y - diff41)^2) / (length(diff41) - ad_df))
    Xg[k] <- sigma
  }

  # Handle zero variance estimates
  pos_Xg <- Xg[which(Xg > 0)]
  if (length(pos_Xg) > 0) {
    Xg[which(Xg == 0)] <- min(pos_Xg)
  } else {
    Xg[which(Xg == 0)] <- 1e-6
  }

  # Handle NaN values
  Xg[is.nan(Xg)] <- 1e-6

  # Histogram
  hist(Xg^2, xlab = expression(paste(hat(sigma)^2)), main = "")
  CHQBC_1_adjB <- Xg

  if (quantile > 0) {
    CHQBC_1_adjB <- Xg + quantile(Xg, quantile, na.rm = TRUE)
  }

  Ts_yvec <- numeric(n_genes)
  Sign <- numeric(n_genes)

  # Start Kernel Smoothing
  for (k in 1:n_genes) {
    # Normalized variance
    diff <- data[k, ] / CHQBC_1_adjB[k]
    d41_y <- nw_smoother(x, diff, x, hwidth)
    diffy4y1 <- d41_y

    # Build the test statistics
    ytest.vec <- diffy4y1
    Ts_yvec[k] <- mean(ytest.vec^2)
    Sign[k] <- sign(sum(ytest.vec))
  }

  # Compute eigenvalue statistics using Rcpp
  eig_stats <- compute_eigenvalue_stats(Snw, n)
  d <- eig_stats$d
  delta <- eig_stats$delta

  Test.adj <- Ts_yvec
  TS_kn <- (((Test.adj / (delta * d))^(1/3) - (1 - 2 / (9 * d))) / sqrt(2 / (9 * d)))

  return(list("TS" = TS_kn, "TS_sign" = Sign, "Tmean" = Ts_yvec))
}

#' Two-sample Kernel Test using Rcpp acceleration
#'
#' Compute a kernel-smoothed statistic for comparing two expression data from ChIPseq.
#' This function applies three kernel smoothing methods to each sample and compares
#' their distributions using a nonparametric statistic.
#'
#' @param data1 Sequence 1 (genes x bins matrix)
#' @param data4 Sequence 2 (genes x bins matrix, same dimensions as data1)
#' @param tao Bias term (from est.c)
#' @param band Kernel bandwidth
#' @param quant Shrinkage variance (vector of 3 values)
#' @param alpha Significance level for hypothesis testing (default: 0.05)
#' @param plot Logical, whether to generate summary plots (default: FALSE)
#' @param plot_genes Integer vector of gene indices to plot in detail (default: NULL, plots first significant gene)
#' @param gene_names Optional character vector of gene names for plotting (default: NULL, uses indices)
#' @param condition_names Optional character vector of length 2 for condition names (default: c("Condition 1", "Condition 2"))
#' @return List containing statistics, p-values, significant genes, and optionally plots
#' @export
#' @examples
#' \dontrun{
#' data(data1)
#' data(data4)
#' tao <- est.c(data1, data4)
#'
#' # Basic usage
#' results <- TS_twosample(data1, data4, tao, band=180, quant=c(0.01,0.01,0.01))
#'
#' # With visualization
#' results <- TS_twosample(data1, data4, tao, band=180, quant=c(0.01,0.01,0.01),
#'                         plot=TRUE, alpha=0.05)
#'
#' # Plot specific genes with custom names
#' results <- TS_twosample(data1, data4, tao, band=180, quant=c(0.01,0.01,0.01),
#'                         plot=TRUE, plot_genes=c(1,3),
#'                         gene_names=c("GeneA","GeneB","GeneC","GeneD","GeneE"),
#'                         condition_names=c("Control", "Treatment"))
#' }
TS_twosample <- function(data1, data4, tao, band, quant,
                         alpha = 0.05,
                         plot = TRUE,
                         plot_genes = NULL,
                         gene_names = NULL,
                         condition_names = c("Condition 1", "Condition 2")) {
  n <- ncol(data1)
  hwidth <- band / n

  # browser()

  # Compute hat matrix using Rcpp
  Snw <- compute_hat_matrix(n, hwidth)

  # Compute all variance estimates using Rcpp
  var_estimates <- compute_variance_estimates(data1, data4, band, Snw)

  # Extract results
  Xg <- var_estimates$Xg
  Ts_yvec <- var_estimates$Ts_yvec
  M <- var_estimates$M
  sigma1 <- var_estimates$sigma1
  sigma4 <- var_estimates$sigma4
  Sev <- var_estimates$Sev
  Suv <- var_estimates$Suv

  # Handle non-positive variance estimates
  pos_Suv <- Suv[which(Suv > 0)]
  pos_Sev <- Sev[which(Sev > 0)]
  pos_Xg <- Xg[which(Xg > 0 & !is.nan(Xg))]

  if (length(pos_Suv) > 0) {
    Suv[which(Suv <= 0)] <- min(pos_Suv)
  } else {
    Suv[which(Suv <= 0)] <- 1e-6
  }

  if (length(pos_Sev) > 0) {
    Sev[which(Sev <= 0)] <- min(pos_Sev)
  } else {
    Sev[which(Sev <= 0)] <- 1e-6
  }

  if (length(pos_Xg) > 0) {
    Xg[which(Xg <= 0 | is.nan(Xg))] <- min(pos_Xg)
  } else {
    Xg[which(Xg <= 0 | is.nan(Xg))] <- 1e-6
  }

  Sev_a0 <- sqrt(Sev)
  Suv_a0 <- sqrt(Suv)
  CHQBC_1_adjB <- Xg

  if (quant[1] > 0) {
    Sev_a0 <- sqrt(Sev) + quantile(sqrt(Sev), quant[1], na.rm = TRUE)
  }
  if (quant[2] > 0) {
    Suv_a0 <- sqrt(Suv) + quantile(sqrt(Suv), quant[2], na.rm = TRUE)
  }
  if (quant[3] > 0) {
    CHQBC_1_adjB <- Xg + quantile(Xg, quant[3], na.rm = TRUE)
  }

  Dsum <- M
  Deql <- (Dsum - tao) / (Sev_a0 / sqrt(n - 1))
  Dnun <- (Dsum - tao) / (Suv_a0 / sqrt(n - 1))

  # Compute eigenvalue statistics using Rcpp
  eig_stats <- compute_eigenvalue_stats(Snw, n)
  d <- eig_stats$d
  delta <- eig_stats$delta

  Test.adj <- Ts_yvec / (CHQBC_1_adjB^2)
  TS_kn <- (((Test.adj / (delta * d))^(1/3) - (1 - 2 / (9 * d))) / sqrt(2 / (9 * d)))

  # Identify significant genes
  z_critical <- qnorm(1 - alpha/2)
  significant_genes <- which(abs(TS_kn) > z_critical)
  p_values <- 2 * (1 - pnorm(abs(TS_kn)))

  # Set up gene names if not provided
  n_genes <- nrow(data1)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene", 1:n_genes)
  }

  # Create results list
  results <- list(
    "sigma1" = sigma1,
    "sigma4" = sigma4,
    "TS_kn" = TS_kn,
    "Ts_yvec" = Ts_yvec,
    "Dsum" = Dsum,
    "Deql" = Deql,
    "Dnun" = Dnun,
    "Sev" = Sev,
    "Suv" = Suv,
    "Xg" = Xg,
    "p_values" = p_values,
    "significant_genes" = significant_genes,
    "alpha" = alpha,
    "z_critical" = z_critical
  )

  # Generate plots if requested
  if (plot) {
    cat("\n=============================================================================\n")
    cat("Generating Visualizations\n")
    cat("=============================================================================\n\n")

    # Summary plot (4 panels)
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))

    par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))

    # Panel 1: TS_kn statistics
    plot(TS_kn, type = "h", lwd = 2, col = "steelblue",
         main = "Kernel Test Statistics (TS_kn)",
         xlab = "Gene Index", ylab = "Test Statistic",
         ylim = range(c(TS_kn, -z_critical, z_critical)),
         xaxt = "n")
    axis(1, at = 1:n_genes, labels = 1:n_genes)
    abline(h = c(-z_critical, z_critical), col = "red", lty = 2, lwd = 1.5)
    abline(h = 0, col = "gray50")
    if (length(significant_genes) > 0) {
      points(significant_genes, TS_kn[significant_genes], pch = 19, col = "red", cex = 1.5)
    }
    legend("topright", legend = c("Test statistic", "Critical value", "Significant"),
           col = c("steelblue", "red", "red"), lty = c(1, 2, NA),
           pch = c(NA, NA, 19), lwd = c(2, 1.5, NA), cex = 0.8)

    # Panel 2: Equal variance test
    plot(Deql, type = "h", lwd = 2, col = "forestgreen",
         main = "Equal Variance Test (Deql)",
         xlab = "Gene Index", ylab = "Test Statistic",
         xaxt = "n")
    axis(1, at = 1:n_genes, labels = 1:n_genes)
    abline(h = c(-z_critical, z_critical), col = "red", lty = 2)
    abline(h = 0, col = "gray50")

    # Panel 3: Unequal variance test
    plot(Dnun, type = "h", lwd = 2, col = "darkorange",
         main = "Unequal Variance Test (Dnun)",
         xlab = "Gene Index", ylab = "Test Statistic",
         xaxt = "n")
    axis(1, at = 1:n_genes, labels = 1:n_genes)
    abline(h = c(-z_critical, z_critical), col = "red", lty = 2)
    abline(h = 0, col = "gray50")

    # Panel 4: Variance comparison
    plot(sigma1, sigma4, pch = 19, col = "purple",
         main = "Variance Comparison",
         xlab = paste("Variance (", condition_names[1], ")", sep=""),
         ylab = paste("Variance (", condition_names[2], ")", sep=""))
    abline(0, 1, col = "gray50", lty = 2)
    text(sigma1, sigma4, labels = 1:length(sigma1), pos = 3, cex = 0.8)

    # Print significance summary
    cat("\nSignificance Summary (alpha =", alpha, "):\n")
    cat("  Critical value:", round(z_critical, 4), "\n")
    cat("  Number of significant genes:", length(significant_genes), "out of", n_genes, "\n")
    if (length(significant_genes) > 0) {
      cat("  Significant gene indices:", significant_genes, "\n")
      cat("  Gene names:", gene_names[significant_genes], "\n")
    }

    # Detailed gene plots
    if (!is.null(plot_genes)) {
      genes_to_plot <- plot_genes
    } else if (length(significant_genes) > 0) {
      genes_to_plot <- significant_genes[1]  # Plot first significant gene
    } else {
      genes_to_plot <- 1  # Default to first gene
    }

    # Plot each selected gene in detail
    for (gene_idx in genes_to_plot) {
      if (gene_idx > n_genes || gene_idx < 1) {
        warning(paste("Gene index", gene_idx, "is out of range. Skipping."))
        next
      }

      par(mfrow = c(1, 1), mar = c(4, 4, 3, 1))

      cat("\n-----------------------------------------------------------------------------\n")
      cat("Detailed Analysis:", gene_names[gene_idx], "(Index:", gene_idx, ")\n")
      cat("-----------------------------------------------------------------------------\n")
      cat("  TS_kn statistic:", round(TS_kn[gene_idx], 4), "\n")
      cat("  P-value:", format.pval(p_values[gene_idx], digits = 4), "\n")
      cat("  Significant:", ifelse(gene_idx %in% significant_genes, "YES", "NO"), "\n")
      cat("  Deql statistic:", round(Deql[gene_idx], 4), "\n")
      cat("  Dnun statistic:", round(Dnun[gene_idx], 4), "\n")
      cat("  Variance (", condition_names[1], "):", round(sigma1[gene_idx], 4), "\n")
      cat("  Variance (", condition_names[2], "):", round(sigma4[gene_idx], 4), "\n")
      cat("  Mean squared difference:", round(Ts_yvec[gene_idx], 4), "\n")

      # Plot ChIP-seq profile
      plot(1:ncol(data1), data1[gene_idx, ], type = "l", col = "blue", lwd = 2,
           main = paste(gene_names[gene_idx], "ChIP-seq Profile"),
           xlab = "Genomic Position (bin)", ylab = "Signal Intensity",
           ylim = range(c(data1[gene_idx, ], data4[gene_idx, ])))
      lines(1:ncol(data4), data4[gene_idx, ], col = "red", lwd = 2)
      legend("topright",
             legend = c(paste(condition_names[1], "(Case)"),
                        paste(condition_names[2], "(Control)")),
             col = c("blue", "red"), lwd = 2, cex = 0.9)

      # Add significance annotation
      if (gene_idx %in% significant_genes) {
        mtext(paste("*** Significant (p =", format.pval(p_values[gene_idx], digits = 3), ")"),
              side = 3, line = 0.5, col = "red", cex = 0.8)
      }
    }

    cat("\n=============================================================================\n")
    cat("Visualization Complete!\n")
    cat("=============================================================================\n\n")
  }

  return(results)
}
