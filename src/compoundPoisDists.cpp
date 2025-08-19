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