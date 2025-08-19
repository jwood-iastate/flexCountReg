#include <Rcpp.h>
#include <cmath> // Include the cmath header

using namespace Rcpp;
using namespace std; // Optional
using Rcpp::NumericVector;


// [[Rcpp::export]]
double dplind_cpp(int x, double mean, double theta) {
  double lambda = mean * theta * (theta + 1) / (theta + 2);
  double theta_sq = theta * theta;
  double theta_plus_lambda = theta + lambda;
  
  // Use exp and log for numerical stability with large x
  double log_result = std::log(theta_sq) + x * std::log(lambda) + 
    std::log(theta + lambda + x + 1) - 
    std::log(theta + 1) - (x + 2) * std::log(theta_plus_lambda);
  
  return std::exp(log_result);
}

// [[Rcpp::export]]
NumericVector dplindlogn_cpp(NumericVector x, NumericVector mean, NumericVector theta, 
                             NumericVector sigma, NumericVector h) {
  int nx = x.size();
  int nmean = mean.size();
  int nsigma = sigma.size();
  int ntheta = theta.size();
  int nh = h.size();
  
  // Input validation
  if (nmean != 1 && nmean != nx) {
    Rcpp::stop("mean must be of length 1 or the same length as x");
  }
  if (nsigma != 1 && nsigma != nx) {
    Rcpp::stop("sigma must be of length 1 or the same length as x");
  }
  if (ntheta != 1 && ntheta != nx) {
    Rcpp::stop("theta must be of length 1 or the same length as x");
  }
  
  NumericVector result(nx);
  
  // Pre-compute exp(h * sigma) for all unique sigma values
  NumericVector unique_sigmas;
  std::vector<NumericVector> lnormdist_cache;
  
  if (nsigma == 1) {
    // Single sigma - compute once
    unique_sigmas.push_back(sigma[0]);
    lnormdist_cache.push_back(exp(h * sigma[0]));
  } else {
    // Multiple sigmas - cache unique values
    for (int i = 0; i < nx; ++i) {
      bool found = false;
      for (int j = 0; j < unique_sigmas.size(); ++j) {
        if (std::abs(sigma[i] - unique_sigmas[j]) < 1e-10) {
          found = true;
          break;
        }
      }
      if (!found) {
        unique_sigmas.push_back(sigma[i]);
        lnormdist_cache.push_back(exp(h * sigma[i]));
      }
    }
  }
  
  // Main computation loop
#pragma omp parallel for if(nx > 1000)
  for (int i = 0; i < nx; ++i) {
    double current_mean = (nmean == 1) ? mean[0] : mean[i];
    double current_sigma = (nsigma == 1) ? sigma[0] : sigma[i];
    double current_theta = (ntheta == 1) ? theta[0] : theta[i];
    
    // Early exit for invalid inputs
    if (x[i] < 0 || current_mean <= 0 || current_sigma <= 0 || current_theta <= 0) {
      result[i] = 0.0;
      continue;
    }
    
    // Find cached lnormdist
    NumericVector* lnormdist = nullptr;
    if (nsigma == 1) {
      lnormdist = &lnormdist_cache[0];
    } else {
      for (int j = 0; j < unique_sigmas.size(); ++j) {
        if (std::abs(current_sigma - unique_sigmas[j]) < 1e-10) {
          lnormdist = &lnormdist_cache[j];
          break;
        }
      }
    }
    
    double p = 0.0;
    
    // Vectorized summation
    for (int j = 0; j < nh; ++j) {
      double mu_j = current_mean * (*lnormdist)[j];
      p += dplind_cpp(x[i], mu_j, current_theta);
    }
    
    result[i] = p / nh;
  }
  
  return result;
}

double poisson_pmf(const double k, const double lambda) {
  return exp(k * log(lambda) - lgamma(k + 1.0) - lambda);
}

// Function to calculate the probability of a Poisson-Lognormal for multiple values
// [[Rcpp::export]]
NumericVector dpLnorm_cpp(NumericVector x, NumericVector mean, NumericVector sigma, NumericVector h) {
  int nx = x.size();
  int nmean = mean.size();
  int nsigma = sigma.size();
  int nh = h.size();
  
  // Input validation
  if (nx == 0) {
    Rcpp::stop("x cannot be empty");
  }
  if (nh == 0) {
    Rcpp::stop("h cannot be empty");
  }
  
  // Ensure mean and sigma are either of length 1 or the same length as x
  if (nmean != 1 && nmean != nx) {
    Rcpp::stop("mean must be of length 1 or the same length as x");
  }
  if (nsigma != 1 && nsigma != nx) {
    Rcpp::stop("sigma must be of length 1 or the same length as x");
  }
  
  NumericVector result(nx);
  
  // Pre-compute exp(h * sigma) values to avoid repeated computation
  NumericVector exp_h_sigma(nh);
  for (int k = 0; k < nh; ++k) {
    exp_h_sigma[k] = exp(h[k]);
  }
  
  for (int i = 0; i < nx; ++i) {
    double current_mean = (nmean == 1) ? mean[0] : mean[i];
    double current_sigma = (nsigma == 1) ? sigma[0] : sigma[i];
    
    // Input validation for current iteration
    if (x[i] < 0 || current_mean <= 0 || current_sigma <= 0) {
      result[i] = 0.0;
      continue;
    }
    
    // Check for non-integer x values (since this is for Poisson)
    if (x[i] != floor(x[i])) {
      result[i] = 0.0;
      continue;
    }
    
    double p = 0.0;
    
    // Vectorized computation of lognormal random effects
    for (int j = 0; j < nh; ++j) {
      double lognormal_effect = pow(exp_h_sigma[j], current_sigma);
      double mu_j = current_mean * lognormal_effect;
      
      // Add overflow protection
      if (mu_j > 1e10) {
        Rcpp::warning("Very large mu detected, results may be unreliable");
        continue;
      }
      
      p += poisson_pmf(x[i], mu_j);
    }
    
    result[i] = p / nh;
  }
  
  return result;
}

// Function to calculate the probability of a Poisson-Weibull mixture for multiple values
// [[Rcpp::export]]
NumericVector dpWeib_cpp(NumericVector x, NumericVector mean, NumericVector alpha, 
                         NumericVector sigma, NumericVector h) {
  int nx = x.size();
  int nmean = mean.size();
  int nsigma = sigma.size();
  int nalpha = alpha.size();
  int nh = h.size();
  
  // Input validation with consistent parameter names
  if (nmean != 1 && nmean != nx) {
    Rcpp::stop("mean must be of length 1 or the same length as x");
  }
  if (nsigma != 1 && nsigma != nx) {
    Rcpp::stop("sigma must be of length 1 or the same length as x");
  }
  if (nalpha != 1 && nalpha != nx) {
    Rcpp::stop("alpha must be of length 1 or the same length as x");
  }
  if (nh == 0) {
    Rcpp::stop("h must have at least one element");
  }
  
  NumericVector result(nx);
  
  // Pre-compute gamma values if alpha is constant
  double gamma_factor = 0.0;
  bool alpha_constant = (nalpha == 1);
  if (alpha_constant) {
    if (alpha[0] <= 0) {
      Rcpp::stop("alpha must be positive");
    }
    gamma_factor = std::tgamma(1.0 + 1.0 / alpha[0]);
  }
  
  for (int i = 0; i < nx; ++i) {
    double current_x = x[i];
    double current_mean = (nmean == 1) ? mean[0] : mean[i];
    double current_alpha = (nalpha == 1) ? alpha[0] : alpha[i];
    double current_sigma = (nsigma == 1) ? sigma[0] : sigma[i];
    
    // Enhanced input validation
    if (current_x < 0) {
      result[i] = 0.0;
      continue;
    }
    if (current_mean <= 0 || current_alpha <= 0 || current_sigma <= 0) {
      result[i] = R_NaN;
      continue;
    }
    
    // Use pre-computed gamma factor when possible
    double current_gamma = alpha_constant ? gamma_factor : std::tgamma(1.0 + 1.0 / current_alpha);
    double lambda = current_mean / (current_sigma * current_gamma);
    
    // Calculate mixture probability with improved numerical stability
    double log_prob_sum = R_NegInf;
    
    for (int j = 0; j < nh; ++j) {
      double current_h = h[j];
      
      // Input validation for h values
      if (current_h <= 0 || current_h >= 1) {
        continue; // Skip invalid h values
      }
      
      // Calculate Weibull quantile
      double weibull_quantile = R::qweibull(current_h, current_alpha, current_sigma, true, false);
      double mu_j = lambda * weibull_quantile;
      
      // Check for numerical issues
      if (mu_j <= 0 || !R_finite(mu_j)) {
        continue;
      }
      
      // Use log-sum-exp trick for better numerical stability
      double log_prob = R::dpois(current_x, mu_j, true); // log=TRUE
      
      if (R_finite(log_prob)) {
        if (log_prob_sum == R_NegInf) {
          log_prob_sum = log_prob;
        } else {
          // Log-sum-exp: log(exp(a) + exp(b)) = a + log(1 + exp(b-a)) when a >= b
          double max_log = std::max(log_prob_sum, log_prob);
          double min_log = std::min(log_prob_sum, log_prob);
          log_prob_sum = max_log + std::log(1.0 + std::exp(min_log - max_log));
        }
      }
    }
    
    // Convert back from log space and normalize
    if (log_prob_sum == R_NegInf) {
      result[i] = 0.0;
    } else {
      result[i] = std::exp(log_prob_sum - std::log(static_cast<double>(nh)));
    }
  }
  
  return result;
}


// Helper function for numerically stable log-sum-exp
double log_sum_exp(double log_a, double log_b, bool subtract = false) {
  if (log_a == -std::numeric_limits<double>::infinity()) {
    return subtract ? -log_b : log_b;
  }
  if (log_b == -std::numeric_limits<double>::infinity()) {
    return log_a;
  }
  
  double max_log = std::max(log_a, log_b);
  double min_log = std::min(log_a, log_b);
  
  if (subtract) {
    return max_log + std::log(std::exp(log_a - max_log) - std::exp(log_b - max_log));
  } else {
    return max_log + std::log(std::exp(log_a - max_log) + std::exp(log_b - max_log));
  }
}

// Function to calculate the probability of a Negative Binomial - Generalized Exponential Model for multiple values
// [[Rcpp::export]]
NumericVector dnbGE_cpp(IntegerVector x, NumericVector mean, NumericVector alpha, NumericVector beta) {
  int nx = x.size();
  int nmean = mean.size();
  int nalpha = alpha.size();
  int nbeta = beta.size();
  
  // Vectorization parameter validation
  if (nmean != 1 && nmean != nx) {
    Rcpp::stop("mean must be of length 1 or the same length as x");
  }
  if (nalpha != 1 && nalpha != nx) {
    Rcpp::stop("alpha must be of length 1 or the same length as x");
  }
  if (nbeta != 1 && nbeta != nx) {
    Rcpp::stop("beta must be of length 1 or the same length as x");
  }
  
  NumericVector result(nx);
  
  for (int i = 0; i < nx; ++i) {
    double mu = (nmean == 1) ? mean[0] : mean[i];
    double alpha_val = (nalpha == 1) ? alpha[0] : alpha[i];
    double beta_val = (nbeta == 1) ? beta[0] : beta[i];
    int x_val = x[i];
    
    // Input validation
    if (x_val < 0 || mu <= 0 || alpha_val <= 0 || beta_val <= 0) {
      result[i] = 0.0;
      continue;
    }
    
    // Check for numerical stability conditions
    double inv_beta = 1.0 / beta_val;
    if (alpha_val <= inv_beta || inv_beta >= 1.0) {
      result[i] = 0.0;
      continue;
    }
    
    // Calculate r parameter using log-space arithmetic for stability
    double log_r_num = std::log(mu) + std::lgamma(alpha_val - inv_beta + 1.0);
    double log_r_den = std::log(std::tgamma(alpha_val + 1.0) * std::tgamma(1.0 - inv_beta) - 
                                std::tgamma(alpha_val - inv_beta + 1.0));
    double r = std::exp(log_r_num - log_r_den);
    
    // Ensure r is positive and finite
    if (r <= 0.0 || !std::isfinite(r)) {
      result[i] = 0.0;
      continue;
    }
    
    // Calculate negative binomial coefficient in log space
    double log_nb_coeff = std::lgamma(r + x_val) - std::lgamma(x_val + 1) - std::lgamma(r);
    
    // Calculate the alternating sum
    double log_sum = -std::numeric_limits<double>::infinity();
    
    for (int j = 0; j <= x_val; ++j) {
      // Binomial coefficient in log space
      double log_binom_coeff = std::lgamma(x_val + 1) - std::lgamma(j + 1) - std::lgamma(x_val - j + 1);
      
      // Gamma ratio in log space
      double log_gamma_ratio = std::lgamma(alpha_val + 1.0) + std::lgamma(1.0 + (r + j) / beta_val) - 
        std::lgamma(alpha_val + (r + j) / beta_val + 1.0);
      
      double log_term = log_binom_coeff + log_gamma_ratio;
      
      // Handle alternating signs using log-sum-exp trick
      if (j % 2 == 0) {
        // Positive term
        log_sum = log_sum_exp(log_sum, log_term);
      } else {
        // Negative term
        log_sum = log_sum_exp(log_sum, log_term, true);
      }
    }
    
    // Combine results
    result[i] = std::exp(log_nb_coeff + log_sum);
    
    // Final check for numerical validity
    if (!std::isfinite(result[i])) {
      result[i] = 0.0;
    }
  }
  
  return result;
}



// Helper function to validate parameter vector sizes
void validateParameterSizes(int nx, int nmean, int nr, int ntheta, int ngamma) {
  if (nmean != 1 && nmean != nx) {
    Rcpp::stop("mean must be of length 1 or the same length as x");
  }
  if (nr != 1 && nr != nx) {
    Rcpp::stop("r must be of length 1 or the same length as x");
  }
  if (ntheta != 1 && ntheta != nx) {
    Rcpp::stop("theta must be of length 1 or the same length as x");
  }
  if (ngamma != 1 && ngamma != nx) {
    Rcpp::stop("gamma must be of length 1 or the same length as x");
  }
}

// Helper function to get the appropriate parameter value
template<typename T>
double getParameter(const T& param, int param_size, int index) {
  return (param_size == 1) ? param[0] : param[index];
}

// Helper function to validate input parameters
bool isValidInput(int x, double r, double theta, double gamma) {
  return (x >= 0 && 
          r > 0 && 
          theta > -0.5 && 
          theta < 0.5 &&
          !std::isnan(gamma) && 
          std::isfinite(gamma));
}

// Helper function to compute discriminant and check validity
struct DiscriminantResult {
  double value;
  bool valid;
  
  DiscriminantResult(double val) : value(val), valid(val > 0.0) {}
};

// Helper function to compute lambda parameter
struct LambdaResult {
  double value;
  bool valid;
  
  LambdaResult(double mean, double r, double theta, double gamma) {
    const double discriminant = 1.0 - 2.0 * theta;
    
    if (discriminant <= 0.0) {
      valid = false;
      return;
    }
    
    const double sqrt_term = std::sqrt(discriminant);
    const double denom = 1.0 - sqrt_term;
    
    if (std::abs(denom) < 1e-12) {  // More robust zero check
      valid = false;
      return;
    }
    
    const double gamma_term = 1.0 - gamma * denom;
    const double numerator = (mean + r) * sqrt_term;
    const double denominator = r * gamma_term;
    
    if (denominator <= 0.0 || numerator <= 0.0) {
      valid = false;
      return;
    }
    
    value = std::log(numerator / denominator) / denom;
    valid = std::isfinite(value);
  }
};

// Helper function to compute binomial coefficient using log-gamma
double logBinomialCoeff(int n, int k) {
  if (k > n || k < 0) return -INFINITY;
  if (k == 0 || k == n) return 0.0;
  
  return std::lgamma(n + 1) - std::lgamma(k + 1) - std::lgamma(n - k + 1);
}

// Helper function to compute negative binomial coefficient using log-gamma
double logNegBinomialCoeff(double r, int x) {
  return std::lgamma(r + x) - std::lgamma(x + 1) - std::lgamma(r);
}

// Numerically stable log-sum-exp function
double logSumExp(double log_a, double log_b) {
  if (log_a == -INFINITY) return log_b;
  if (log_b == -INFINITY) return log_a;
  
  const double max_val = std::max(log_a, log_b);
  const double min_val = std::min(log_a, log_b);
  
  if (max_val - min_val > 30.0) {  // Prevent underflow
    return max_val;
  }
  
  return max_val + std::log1p(std::exp(min_val - max_val));
}



// Numerically stable log-diff-exp function
double logDiffExp(double log_a, double log_b) {
  if (log_a == -INFINITY) return log_b;
  if (log_b == -INFINITY) return log_a;
  
  if (log_a < log_b) {
    return log_b + std::log1p(-std::exp(log_a - log_b));
  } else {
    return log_a + std::log1p(-std::exp(log_b - log_a));
  }
}
// Main density computation function
double computeDensity(int x, double mean, double r, double theta, double gamma) {
  // Validate input parameters
  if (!isValidInput(x, r, theta, gamma)) {
    return NA_REAL;
  }
  
  // Compute lambda parameter
  LambdaResult lambda_result(mean, r, theta, gamma);
  if (!lambda_result.valid) {
    return NA_REAL;
  }
  
  const double lambda = lambda_result.value;
  
  // Compute the negative binomial coefficient in log space
  const double log_nb_coeff = logNegBinomialCoeff(r, x);
  if (!std::isfinite(log_nb_coeff)) {
    return NA_REAL;
  }
  
  // Compute the sum using more stable numerical methods
  double log_sum = -INFINITY;  // Start with log(0)
  
  for (int j = 0; j <= x; ++j) {
    // Compute binomial coefficient in log space
    const double log_binom_coeff = logBinomialCoeff(x, j);
    if (!std::isfinite(log_binom_coeff)) {
      continue;
    }
    
    const double temp = r + j;
    const double discriminant_j = 1.0 + 2.0 * theta * temp;
    
    if (discriminant_j <= 0.0) {
      return NA_REAL;
    }
    
    const double sqrt_term_j = std::sqrt(discriminant_j);
    const double exp_term = lambda * (1.0 - sqrt_term_j);
    const double gamma_term_j = 1.0 - gamma * (1.0 - sqrt_term_j);
    
    if (sqrt_term_j == 0.0 || gamma_term_j <= 0.0) {
      continue;
    }
    
    // Compute term in log space
    const double log_term = log_binom_coeff + exp_term + 
      std::log(std::abs(gamma_term_j)) - std::log(sqrt_term_j);
    
    if (!std::isfinite(log_term)) {
      continue;
    }
    
    // Handle alternating signs and accumulate using log-sum-exp
    const double sign = (j % 2 == 0) ? 1.0 : -1.0;
    if (sign > 0) {
      log_sum = logSumExp(log_sum, log_term);
    } else {
      log_sum = logDiffExp(log_sum, log_term);
    }
  }
  
  if (!std::isfinite(log_sum)) {
    return NA_REAL;
  }
  
  const double log_result = log_nb_coeff + log_sum;
  return std::isfinite(log_result) ? std::exp(log_result) : NA_REAL;
}




// Negative Binomial Crack Distribution - Improved Version
// [[Rcpp::export]]
NumericVector dnbCrack_cpp(const IntegerVector& x, const NumericVector& mean, 
                           const NumericVector& r, const NumericVector& theta, 
                           const NumericVector& gamma) {
  const int nx = x.size();
  const int nmean = mean.size();
  const int nr = r.size();
  const int ntheta = theta.size();
  const int ngamma = gamma.size();
  
  // Validate parameter vector sizes
  validateParameterSizes(nx, nmean, nr, ntheta, ngamma);
  
  NumericVector result(nx);
  
  for (int i = 0; i < nx; ++i) {
    result[i] = computeDensity(x[i], 
                               getParameter(mean, nmean, i),
                               getParameter(r, nr, i),
                               getParameter(theta, ntheta, i),
                               getParameter(gamma, ngamma, i));
  }
  
  return result;
}



