#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector dplindgamma_cpp(const NumericVector& x, 
                              const NumericVector& mean, 
                              const NumericVector& theta, 
                              const NumericVector& alpha, 
                              const NumericVector& h) {
  
  // Get sizes
  const int nx = x.size();
  const int nmean = mean.size();
  const int nalpha = alpha.size();
  const int ntheta = theta.size();
  const int nh = h.size();
  
  // Validate input lengths
  if (nmean != 1 && nmean != nx) {
    stop("mean must be of length 1 or the same length as x");
  }
  if (nalpha != 1 && nalpha != nx) {
    stop("alpha must be of length 1 or the same length as x");
  }
  if (ntheta != 1 && ntheta != nx) {
    stop("theta must be of length 1 or the same length as x");
  }
  if (nh == 0) {
    stop("h must have at least one element");
  }
  
  // Pre-allocate result vector
  NumericVector result(nx);
  
  // Precompute 1/nh for efficiency
  const double inv_nh = 1.0 / nh;
  
  for (int i = 0; i < nx; ++i) {
    // Get current parameter values
    const double current_mean = (nmean == 1) ? mean[0] : mean[i];
    const double current_alpha = (nalpha == 1) ? alpha[0] : alpha[i];
    const double current_theta = (ntheta == 1) ? theta[0] : theta[i];
    const double current_x = x[i];
    
    // Check for invalid parameters or negative x
    if (current_x < 0 || current_mean <= 0 || current_alpha <= 0 || current_theta <= 0) {
      result[i] = 0.0;
      continue;
    }
    
    // Check for non-finite values
    if (!R_FINITE(current_x) || !R_FINITE(current_mean) || 
        !R_FINITE(current_alpha) || !R_FINITE(current_theta)) {
        result[i] = NA_REAL;
      continue;
    }
    
    double p = 0.0;
    
    // Precompute values that don't depend on j
    const double theta_squared = current_theta * current_theta;
    const double theta_plus_one = current_theta + 1.0;
    const double scale_factor = 1.0 / current_alpha;
    
    for (int j = 0; j < nh; ++j) {
      const double current_h = h[j];
      
      // Skip invalid h values
      if (current_h <= 0.0 || current_h >= 1.0 || !R_FINITE(current_h)) {
        continue;
      }
      
      // Compute quantile
      const double quantile = R::qgamma(current_h, current_alpha, scale_factor, 1, 0);
      const double mu_j = current_mean * quantile;
      
      // Check for numerical issues
      if (!R_FINITE(mu_j) || mu_j <= 0) {
        continue;
      }
      
      const double lambda = mu_j * current_theta * theta_plus_one / (current_theta + 2.0);
      
      // Check lambda validity
      if (!R_FINITE(lambda) || lambda <= 0) {
        continue;
      }
      
      // Compute probability contribution with numerical stability
      const double theta_plus_lambda = current_theta + lambda;
      const double log_numerator = 2.0 * std::log(current_theta) + 
        current_x * std::log(lambda) + 
        std::log(theta_plus_lambda + current_x + 1.0);
      const double log_denominator = std::log(theta_plus_one) + 
        (current_x + 2.0) * std::log(theta_plus_lambda);
      
      const double log_prob = log_numerator - log_denominator;
      
      // Check for numerical overflow/underflow
      if (R_FINITE(log_prob)) {
        p += std::exp(log_prob);
      }
    }
    
    result[i] = p * inv_nh;
  }
  
  return result;
}