#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' Gaussian kernel function
 //' @param x Distance
 //' @param bandwidth Bandwidth parameter
 // [[Rcpp::export]]
 double gaussian_kernel(double x, double bandwidth) {
   return exp(-0.5 * pow(x / bandwidth, 2)) / (bandwidth * sqrt(2.0 * M_PI));
 }

 //' Nadaraya-Watson kernel smoother (single point)
 //' @param x_points All x coordinates
 //' @param y_points All y coordinates
 //' @param x_target Target x coordinate to evaluate
 //' @param bandwidth Bandwidth parameter
 // [[Rcpp::export]]
 double nw_smoother_point(const arma::vec& x_points,
                          const arma::vec& y_points,
                          double x_target,
                          double bandwidth) {
   int n = x_points.n_elem;
   double sum_weights = 0.0;
   double sum_weighted_y = 0.0;

   for (int i = 0; i < n; i++) {
     double w = gaussian_kernel(x_target - x_points(i), bandwidth);
     sum_weights += w;
     sum_weighted_y += w * y_points(i);
   }

   if (sum_weights > 1e-10) {
     return sum_weighted_y / sum_weights;
   } else {
     return 0.0;
   }
 }

 //' Nadaraya-Watson kernel smoother (full vector)
 //' @param x_points All x coordinates
 //' @param y_points All y coordinates
 //' @param x_eval Evaluation points
 //' @param bandwidth Bandwidth parameter
 // [[Rcpp::export]]
 arma::vec nw_smoother(const arma::vec& x_points,
                       const arma::vec& y_points,
                       const arma::vec& x_eval,
                       double bandwidth) {
   int n_eval = x_eval.n_elem;
   arma::vec result(n_eval);

   for (int i = 0; i < n_eval; i++) {
     result(i) = nw_smoother_point(x_points, y_points, x_eval(i), bandwidth);
   }

   return result;
 }

 //' Compute the hat matrix for Nadaraya-Watson estimator
 //' @param n Sample size
 //' @param bandwidth Bandwidth parameter
 // [[Rcpp::export]]
 arma::mat compute_hat_matrix(int n, double bandwidth) {
   arma::mat Snw(n, n, arma::fill::zeros);
   arma::vec x = arma::linspace(1.0/n, 1.0, n);

   // For each column of identity matrix
   for (int j = 0; j < n; j++) {
     arma::vec y = arma::zeros(n);
     y(j) = 1.0;

     // Compute smoothed version
     for (int i = 0; i < n; i++) {
       Snw(i, j) = nw_smoother_point(x, y, x(i), bandwidth);
     }
   }

   return Snw;
 }

 //' Compute variance estimates for all genes
 //' @param data1 First condition data matrix (genes x bins)
 //' @param data4 Second condition data matrix (genes x bins)
 //' @param bandwidth Bandwidth parameter
 //' @param Snw Hat matrix
 // [[Rcpp::export]]
 List compute_variance_estimates(const arma::mat& data1,
                                 const arma::mat& data4,
                                 double bandwidth,
                                 const arma::mat& Snw) {
   int n_genes = data1.n_rows;
   int n = data1.n_cols;
   double hwidth = bandwidth / n;

   arma::vec x = arma::linspace(1.0/n, 1.0, n);
   double ad_df = arma::trace(Snw);

   // Initialize result vectors
   arma::vec Xg(n_genes);
   arma::vec Ts_yvec(n_genes);
   arma::vec M(n_genes);
   arma::vec sigma1(n_genes);
   arma::vec sigma4(n_genes);
   arma::vec sigma41(n_genes);
   arma::vec sigma_1_sq(n_genes);
   arma::vec sigma_4_sq(n_genes);
   arma::vec sigma1_fg(n_genes);
   arma::vec sigma4_fg(n_genes);
   arma::vec sig4sig1(n_genes);
   arma::vec sigma_unequal(n_genes);

   for (int k = 0; k < n_genes; k++) {
     arma::rowvec L1 = data1.row(k);
     arma::rowvec L4 = data4.row(k);
     arma::rowvec Diff = L4 - L1;

     // Convert to vec for smoothing
     arma::vec Diff_vec = arma::conv_to<arma::vec>::from(Diff);
     arma::vec L1_vec = arma::conv_to<arma::vec>::from(L1);
     arma::vec L4_vec = arma::conv_to<arma::vec>::from(L4);

     // Kernel smoothing
     arma::vec d41_y = nw_smoother(x, Diff_vec, x, hwidth);

     // Kernel smoothed variance estimate
     double residual_ss = arma::sum(arma::square(d41_y - Diff_vec));
     Xg(k) = sqrt(residual_ss / (n - ad_df));
     Ts_yvec(k) = arma::mean(arma::square(d41_y));

     // Autocorrelation term M
     double M_sum = 0.0;
     for (int i = 0; i < n - 1; i++) {
       M_sum += Diff(i+1) * Diff(i);
     }
     M(k) = M_sum / (n - 1);

     // Equal variance estimates
     double s1 = 0.0, s4 = 0.0;
     for (int i = 0; i < n - 1; i++) {
       s1 += pow(L1(i+1) - L1(i), 2);
       s4 += pow(L4(i+1) - L4(i), 2);
     }
     sigma1(k) = s1 / (2.0 * (n - 1));
     sigma4(k) = s4 / (2.0 * (n - 1));
     sigma41(k) = pow(sigma1(k) + sigma4(k), 2) + 4.0 * M(k) * (sigma1(k) + sigma4(k));

     // Unequal variance estimates
     if (n >= 4) {
       double s1_sq = 0.0, s4_sq = 0.0, s1_fg = 0.0, s4_fg = 0.0, s4s1 = 0.0;

       for (int i = 0; i < n - 3; i++) {
         s1_sq += pow(L1(i+1) - L1(i), 2) * pow(L1(i+3) - L1(i+2), 2);
         s4_sq += pow(L4(i+1) - L4(i), 2) * pow(L4(i+3) - L4(i+2), 2);
       }
       sigma_1_sq(k) = s1_sq / (4.0 * (n - 3));
       sigma_4_sq(k) = s4_sq / (4.0 * (n - 3));

       for (int i = 0; i < n - 2; i++) {
         double prev_L1 = (i == 0) ? L1(0) : L1(i-1);
         double prev_L4 = (i == 0) ? L4(0) : L4(i-1);
         s1_fg += (L4(i) - L1(i)) * (prev_L4 - prev_L1) * pow(L1(i+2) - L1(i+1), 2);
         s4_fg += (L4(i) - L1(i)) * (prev_L4 - prev_L1) * pow(L4(i+2) - L4(i+1), 2);
       }
       sigma1_fg(k) = s1_fg / (2.0 * (n - 3));
       sigma4_fg(k) = s4_fg / (2.0 * (n - 3));

       for (int i = 0; i < n - 1; i++) {
         s4s1 += pow(L1(i+1) - L1(i), 2) * pow(L4(i+1) - L4(i), 2);
       }
       sig4sig1(k) = s4s1 / (4.0 * (n - 3));

       sigma_unequal(k) = (sigma_1_sq(k) + 4.0 * sigma1_fg(k)) +
         (sigma_4_sq(k) + 4.0 * sigma4_fg(k)) +
         2.0 * sig4sig1(k);
     } else {
       sigma_1_sq(k) = 0.0;
       sigma_4_sq(k) = 0.0;
       sigma1_fg(k) = 0.0;
       sigma4_fg(k) = 0.0;
       sig4sig1(k) = 0.0;
       sigma_unequal(k) = 0.0;
     }
   }

   return List::create(
     Named("Xg") = Xg,
     Named("Ts_yvec") = Ts_yvec,
     Named("M") = M,
     Named("sigma1") = sigma1,
     Named("sigma4") = sigma4,
     Named("Sev") = sigma41,
     Named("Suv") = sigma_unequal
   );
 }

 //' Estimate bias term (tao) from stable genes
 //' @param data1 First condition data matrix
 //' @param data4 Second condition data matrix
 //' @param max1 Maximum threshold for condition 1
 //' @param max4 Maximum threshold for condition 2
 // [[Rcpp::export]]
 double est_c_cpp(const arma::mat& data1,
                  const arma::mat& data4,
                  double max1,
                  double max4) {
   int n_genes = data1.n_rows;
   int n = data1.n_cols;
   arma::vec M(n_genes);

   // Find maximum values for each gene
   arma::vec dat1_max = arma::max(data1, 1);
   arma::vec dat4_max = arma::max(data4, 1);

   // Calculate M for all genes
   for (int k = 0; k < n_genes; k++) {
     arma::rowvec L1 = data1.row(k);
     arma::rowvec L4 = data4.row(k);
     arma::rowvec Diff = L4 - L1;

     double M_sum = 0.0;
     for (int i = 0; i < n - 1; i++) {
       M_sum += Diff(i+1) * Diff(i);
     }
     M(k) = M_sum / n;
   }

   // Find stable genes
   double tao_sum = 0.0;
   int count = 0;
   for (int k = 0; k < n_genes; k++) {
     if (dat1_max(k) <= max1 && dat4_max(k) <= max4) {
       tao_sum += M(k);
       count++;
     }
   }

   if (count > 0) {
     return tao_sum / count;
   } else {
     return 0.0;
   }
 }

 //' Compute eigenvalues and related statistics
 //' @param Snw Hat matrix
 //' @param n Sample size
 // [[Rcpp::export]]
 List compute_eigenvalue_stats(const arma::mat& Snw, int n) {
   arma::mat Amax = Snw * Snw;
   arma::vec eigenvalues = arma::eig_sym(Amax);

   double sum_eig = arma::sum(eigenvalues);
   double sum_eig_sq = arma::sum(arma::square(eigenvalues));

   double d = pow(sum_eig, 2) / sum_eig_sq;
   double delta = sum_eig / (n * d);

   return List::create(
     Named("eigenvalues") = eigenvalues,
     Named("d") = d,
     Named("delta") = delta
   );
 }
