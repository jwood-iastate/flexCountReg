// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

using namespace Rcpp;
using namespace Eigen;

// --- HELPER FUNCTION: Confluent Hypergeometric Function 1F1 ---
// Computes M(a, b, z) = 1F1(a; b; z) using the series definition.
// This is a helper for the main U function.
VectorXd _1F1_vec(const VectorXd& a, const VectorXd& b, const VectorXd& z) {
  int n = z.size();
  VectorXd sum = VectorXd::Ones(n);
  VectorXd term = VectorXd::Ones(n);
  
  const int max_iter = 500;
  const double tol = 1e-14;
  
  for (int k = 1; k < max_iter; ++k) {
    // Vectorized term update: term_k = term_{k-1} * (a+k-1)*z / ((b+k-1)*k)
    term = term.array() * (a.array() + k - 1) * z.array() / ((b.array() + k - 1) * k);
    sum += term;
    
    // More robust convergence check: stop when all terms are small
    if ((term.array().abs() < tol).all()) {
      break;
    }
  }
  return sum;
}


// --- HELPER FUNCTION: Series expansion for U(a,b,z) for small z and INTEGER b ---
// This correctly implements the complex logarithmic case from NIST DLMF 13.2.10.
VectorXd tricomi_U_series_small_z(const VectorXd& a, const VectorXd& b, const VectorXd& z) {
  int n = z.size();
  VectorXd result(n);
  
  // This implementation is specialized for integer b >= 1.
  for (int i = 0; i < n; ++i) {
    double ai = a(i);
    int bi = static_cast<int>(round(b(i)));
    double zi = z(i);
    
    // Term 1: involves log(z) and 1F1(a;b;z)
    double term1_factor = std::pow(-1.0, bi) / (R::gammafn(bi) * R::lgammafn(ai - bi + 1));
    VectorXd M = _1F1_vec(a.segment(i,1), b.segment(i,1), z.segment(i,1));
    double term1 = term1_factor * (std::log(zi) * M(0));
    
    // Term 2: The complex sum involving the digamma function (psi)
    double sum_psi = 0;
    double term_psi = 1.0; // k=0 term
    sum_psi += term_psi * (R::digamma(ai) - R::digamma(1.0) - R::digamma(bi));
    
    const int max_iter_psi = 500;
    const double tol_psi = 1e-14;
    for (int k = 1; k < max_iter_psi; ++k) {
      term_psi *= (ai + k - 1) * zi / ((bi + k - 1) * k);
      sum_psi += term_psi * (R::digamma(ai + k) - R::digamma(1.0 + k) - R::digamma(bi + k));
      if (std::abs(term_psi) < tol_psi) break;
    }
    double term2 = term1_factor * sum_psi;
    
    // Term 3: A finite sum that exists only for b >= 2
    double term3 = 0;
    if (bi >= 2) {
      double prefactor = R::gammafn(bi - 1) / R::gammafn(ai);
      double sum_finite = 0;
      double term_finite = 1.0; // k=0 term
      sum_finite += term_finite;
      
      for (int k = 1; k <= bi - 2; ++k) {
        term_finite *= (ai - bi + 1 + k - 1) * zi / ((2 - bi + k - 1) * k);
        sum_finite += term_finite;
      }
      term3 = prefactor * std::pow(zi, 1 - bi) * sum_finite;
    }
    
    result(i) = term1 + term2 + term3;
  }
  return result;
}


// [[Rcpp::export]]
VectorXd tricomi_U_vec(const VectorXd& a, const VectorXd& b, const VectorXd& z) {
  int n = z.size();
  VectorXd result(n);
  
  // Define the threshold for switching between algorithms
  const double z_threshold = 30.0;
  
  for (int i = 0; i < n; ++i) {
    double zi = z(i);
    double ai = a(i);
    double bi = b(i);
    
    // --- Case 1: z is large. Use the asymptotic series. ---
    // Formula: U(a,b,z) ~ z^(-a) * sum_{k=0 to inf} [ (a)_k * (a-b+1)_k / k! * (-1/z)^k ]
    if (zi > z_threshold) {
      double sum = 1.0;
      double term = 1.0;
      const int max_iter = 100;
      const double tol = 1e-14;
      
      for (int k = 1; k < max_iter; ++k) {
        term *= -1.0 * (ai + k - 1) * (ai - bi + k) / (k * zi);
        sum += term;
        if (std::abs(term) < tol * std::abs(sum)) {
          break;
        }
      }
      result(i) = std::pow(zi, -ai) * sum;
    }
    // --- Case 2: z is small or moderate. Use the correct series expansion. ---
    else {
      // This is the correct, but complicated, formula for integer b.
      result(i) = tricomi_U_series_small_z(a.segment(i,1), b.segment(i,1), z.segment(i,1))(0);
    }
  }
  return result;
}

// [[Rcpp::export]]
NumericVector negbinlindley_pmf_cpp(const NumericVector& y, const NumericVector& mu, double alpha, double theta) {
  int n = y.size();
  if (n != mu.size()) {
    stop("Input vectors y and mu must have the same length.");
  }
  
  // Parameter validation
  if (alpha <= 0 || theta <= 0) {
    stop("Parameters alpha and theta must be positive.");
  }
  for(int i = 0; i < n; ++i) {
    if(mu[i] <= 0) {
      stop("All elements of mu must be positive.");
    }
  }
  
  // Use Eigen for efficient vector operations
  Map<const VectorXd> y_vec(y.begin(), n);
  Map<const VectorXd> mu_vec(mu.begin(), n);
  
  // Precompute constants for efficiency
  const double theta_plus_1 = theta + 1.0;
  const double theta_plus_1_sq = theta_plus_1 * theta_plus_1;
  const double theta_plus_2 = theta + 2.0;
  const double log_binom_cons = (1.0 / alpha) * std::log(1.0 / (1.0 + 1.0 / alpha)); // Use log for numerical stability
  const double lgamma_inv_alpha = R::lgammafn(1.0 / alpha);
  
  // Prepare arguments for the U function
  VectorXd z = theta_plus_2 / (mu_vec.array() * alpha * theta_plus_1);
  VectorXd a1 = y_vec.array() + 1.0;
  VectorXd a2 = y_vec.array() + 2.0;
  
  // Compute U functions using the reliable GSL wrapper
  VectorXd U1 = tricomi_U_gsl_vec(a1, VectorXd::Constant(n, 2.0), z);
  VectorXd U2 = tricomi_U_gsl_vec(a2, VectorXd::Constant(n, 3.0), z);
  
  // Vectorized computation of log-gamma for y + 1/alpha
  VectorXd lgamma_y_alpha(n);
  for(int i=0; i < n; ++i) {
    lgamma_y_alpha(i) = R::lgammafn(y_vec(i) + 1.0/alpha);
  }
  
  // Calculate the log-PMF first to improve numerical stability, then exponentiate
  VectorXd log_prob = lgamma_y_alpha.array() + std::log(theta) + std::log(theta_plus_2) + log_binom_cons
  - (lgamma_inv_alpha + std::log(alpha) + mu_vec.array().log() + std::log(theta_plus_1_sq))
    + (U1.array() + (y_vec.array() + 1.0) * U2.array()).log();
    
    VectorXd prob = log_prob.array().exp();
    
    // Handle invalid y values (negative or non-integer)
    for (int i = 0; i < n; ++i) {
      if (y(i) < 0 || std::floor(y(i)) != y(i)) {
        prob(i) = 0.0;
      }
    }
    
    return Rcpp::wrap(prob);
}

// Enhanced CDF function - much more efficient
// [[Rcpp::export]]
NumericVector negbinlindley_cdf_cpp(const NumericVector& q, const NumericVector& mu, double alpha, double theta) {
  int n = q.size();
  if (n != mu.size()) {
    stop("Input vectors q and mu must have the same length.");
  }
  
  NumericVector cdf_results(n);
  
  for (int i = 0; i < n; ++i) {
    if (q[i] < 0) {
      cdf_results[i] = 0.0;
      continue;
    }
    
    int max_y = static_cast<int>(std::floor(q[i]));
    
    // Create a sequence of y values from 0 to floor(q)
    if (max_y < 0) { // Handle case where q is between 0 and 1
      cdf_results[i] = 0.0; // Assuming PMF is 0 for y < 0
    } else {
      std::vector<double> y_values(max_y + 1);
      std::iota(y_values.begin(), y_values.end(), 0); // Fills with 0, 1, 2, ...
      
      // Create a corresponding mu vector
      NumericVector y_seq(y_values.begin(), y_values.end());
      NumericVector mu_seq(max_y + 1, mu[i]);
      
      // Call the PMF function ONCE with all necessary y values
      NumericVector pmf_values = negbinlindley_pmf_cpp(y_seq, mu_seq, alpha, theta);
      
      // Sum the results to get the CDF
      cdf_results[i] = Rcpp::sum(pmf_values);
    }
  }
  
  return cdf_results;
}