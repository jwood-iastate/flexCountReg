// Functions for the Conway-Maxwell-Poisson distribution
#include <Rcpp.h>
using namespace Rcpp;

// Calculate the normalizing constant for the Conway-Maxwell-Poisson distribution
// [[Rcpp::export]]
double cmp_normalizer_cpp(double lambda, double nu, int maxval = 1000) {
  double sum = 0.0;
  double log_lambda = std::log(lambda);
  double log_fact = 0.0; // log(0!) = 0
  double max_log_term = -INFINITY;
  std::vector<double> log_terms(maxval + 1);
  
  // Compute log terms and find the maximum for numerical stability
  for (int k = 0; k <= maxval; ++k) {
    if (k > 0) {
      log_fact += std::log(k); // Update log(k!) incrementally
    }
    double log_term = k * log_lambda - nu * log_fact;
    log_terms[k] = log_term;
    if (log_term > max_log_term) {
      max_log_term = log_term;
    }
  }
  
  // Compute the sum using the log-sum-exp trick for numerical stability
  for (int k = 0; k <= maxval; ++k) {
    sum += std::exp(log_terms[k] - max_log_term);
  }
  
  sum *= std::exp(max_log_term);
  return sum;
}

// Calculate the expected value of the Conway-Maxwell-Poisson distribution
// [[Rcpp::export]]
double com_expect_cpp(double lambda, double nu, int maxval = 1000) {
  double sum = 0.0;
  double sum_expect = 0.0;
  double log_lambda = std::log(lambda);
  double log_fact = 0.0; // log(0!) = 0
  double max_log_term = -INFINITY;
  std::vector<double> log_terms(maxval + 1);
  
  // Compute log terms and find the maximum for numerical stability
  for (int k = 0; k <= maxval; ++k) {
    if (k > 0) {
      log_fact += std::log(k);
    }
    double log_term = k * log_lambda - nu * log_fact;
    log_terms[k] = log_term;
    if (log_term > max_log_term) {
      max_log_term = log_term;
    }
  }
  
  // Compute the sums using the log-sum-exp trick
  for (int k = 0; k <= maxval; ++k) {
    double term = std::exp(log_terms[k] - max_log_term);
    sum += term;
    sum_expect += k * term;
  }
  
  sum *= std::exp(max_log_term);
  sum_expect *= std::exp(max_log_term);
  
  return sum_expect / sum;
}

// Calculate the variance of the Conway-Maxwell-Poisson distribution
// [[Rcpp::export]]
double com_var_cpp(double lambda, double nu, int maxval = 1000) {
  double sum = 0.0;
  double sum_expect = 0.0;
  double sum_expect2 = 0.0;
  double log_lambda = std::log(lambda);
  double log_fact = 0.0;
  double max_log_term = -INFINITY;
  std::vector<double> log_terms(maxval + 1);
  
  // Compute log terms and find the maximum for numerical stability
  for (int k = 0; k <= maxval; ++k) {
    if (k > 0) {
      log_fact += std::log(k);
    }
    double log_term = k * log_lambda - nu * log_fact;
    log_terms[k] = log_term;
    if (log_term > max_log_term) {
      max_log_term = log_term;
    }
  }
  
  // Compute the sums using the log-sum-exp trick
  for (int k = 0; k <= maxval; ++k) {
    double term = std::exp(log_terms[k] - max_log_term);
    sum += term;
    sum_expect += k * term;
    sum_expect2 += k * k * term;
  }
  
  sum *= std::exp(max_log_term);
  sum_expect *= std::exp(max_log_term);
  sum_expect2 *= std::exp(max_log_term);
  
  double mean = sum_expect / sum;
  double mean2 = sum_expect2 / sum;
  return mean2 - mean * mean; // Variance = E[X^2] - (E[X])^2
}
