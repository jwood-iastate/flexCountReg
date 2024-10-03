#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;


// [[Rcpp::export]]
double dplind_cpp(int x, double mean, double theta) {
  double lambda = mean * theta * (theta + 1) / (theta + 2);
  return (std::pow(theta, 2) * std::pow(lambda, x) * (theta + lambda + x + 1) / ((theta + 1) * std::pow(theta + lambda, x + 2)));
}

// [[Rcpp::export]]
NumericVector dplindlogn_cpp(NumericVector x, NumericVector mean, NumericVector theta, NumericVector sigma, NumericVector h) {
  int nx = x.size();
  int nmean = mean.size();
  int nsigma = sigma.size();
  int ntheta = theta.size();
  int nh = h.size();
  
  // Ensure mean and sigma are either of length 1 or the same length as x
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
  
  for (int i = 0; i < nx; ++i) {
    double current_mean = (nmean == 1) ? mean[0] : mean[i];
    double current_sigma = (nsigma == 1) ? sigma[0] : sigma[i];
    double current_theta = (ntheta == 1) ? theta[0] : theta[i];
    
    if (x[i] < 0 || current_mean <= 0 || current_sigma <= 0 || current_theta <= 0) {
      result[i] = 0.0;
      continue;
    }
    
    NumericVector lnormdist = exp(h * current_sigma);
    double p = 0.0;
    
    for (NumericVector::iterator j = lnormdist.begin(); j != lnormdist.end(); ++j) {
      double mu_j = current_mean * (*j);
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
  
  // Ensure mean and sigma are either of length 1 or the same length as x
  if (nmean != 1 && nmean != nx) {
    Rcpp::stop("mean must be of length 1 or the same length as x");
  }
  if (nsigma != 1 && nsigma != nx) {
    Rcpp::stop("sigma must be of length 1 or the same length as x");
  }
  
  NumericVector result(nx);
  
  for (int i = 0; i < nx; ++i) {
    
    double current_mean = (nmean == 1) ? mean[0] : mean[i];
    double current_sigma = (nsigma == 1) ? sigma[0] : sigma[i];
    
    if (x[i] < 0 || current_mean <= 0 || current_sigma <= 0) {
      result[i] = 0.0;
      continue;
    }
    
    NumericVector lnormdist = exp(h * current_sigma);
    double p = 0.0;
    
    for (NumericVector::iterator j = lnormdist.begin(); j != lnormdist.end(); ++j) {
      double mu_j = current_mean * (*j);
      p += poisson_pmf(x[i], mu_j);
    }
    
    result[i] = p / nh;
  }
  
  return result;
}

// Function to calculate the probability of a Poisson-Weibull for multiple values
// [[Rcpp::export]]
NumericVector dpWeib_cpp(NumericVector x, NumericVector mean, NumericVector alpha,  NumericVector sigma, NumericVector h) {
  int nx = x.size();
  int nmean = mean.size();
  int nsigma = sigma.size();
  int nalpha = alpha.size();
  int nh = h.size();
  
  // Ensure mean and sigma are either of length 1 or the same length as x
  if (nmean != 1 && nmean != nx) {
    Rcpp::stop("lambda must be of length 1 or the same length as x");
  }
  if (nsigma != 1 && nsigma != nx) {
    Rcpp::stop("sigma must be of length 1 or the same length as x");
  }
  if (nalpha != 1 && nalpha != nx) {
    Rcpp::stop("alpha must be of length 1 or the same length as x");
  }
  
  NumericVector result(nx);
  
  for (int i = 0; i < nx; ++i) {
    double current_mean = (nmean == 1) ? mean[0] : mean[i];
    double current_alpha = (nalpha == 1) ? alpha[0] : alpha[i];
    double current_sigma = (nsigma == 1) ? sigma[0] : sigma[i];
    double p = 0.0;
    
    if (x[i] < 0 || current_mean <= 0 || current_alpha <= 0 || current_sigma <= 0) {
      result[i] = 0.0;
      continue;
    }
    
    double lambda = current_mean / (current_sigma*std::tgamma(1+1/current_alpha));
    
    for (int j = 0; j < nh; ++j) {
      double current_h = h[j];
      double mu_j = lambda * R::qweibull(current_h, current_alpha, current_sigma, true, false);
      p += R::dpois(x[i], mu_j, false);
    }
    
    result[i] = p / nh;
  }
  
  return result;
}

// Function to calculate the probability of a Negative Binomial - Generalized Exponential Model for multiple values
// [[Rcpp::export]]
NumericVector dnbGE_cpp(IntegerVector x, NumericVector mean, NumericVector alpha, NumericVector beta) {
  int nx = x.size();
  int nmean = mean.size();
  int nbeta = beta.size();
  int nalpha = alpha.size();
  
  // Ensure mean and parameters are either of length 1 or the same length as x
  if (nmean != 1 && nmean != nx) {
    Rcpp::stop("mean must be of length 1 or the same length as x");
  }
  if (nbeta != 1 && nbeta != nx) {
    Rcpp::stop("beta must be of length 1 or the same length as x");
  }
  if (nalpha != 1 && nalpha != nx) {
    Rcpp::stop("alpha must be of length 1 or the same length as x");
  }
  
  NumericVector result(nx);
  
  for (int i = 0; i < nx; ++i) {
    double current_mean = (nmean == 1) ? mean[0] : mean[i];
    double current_alpha = (nalpha == 1) ? alpha[0] : alpha[i];
    double current_beta = (nbeta == 1) ? beta[0] : beta[i];
    
    if (x[i] < 0 || current_mean <= 0 || current_alpha <= 0 || current_beta <= 0) {
      result[i] = 0.0;
      continue;
    }
    
    int current_x = x[i];
    
    // Calculate r
    double r_numerator = current_mean * std::tgamma(current_alpha - (1.0 / current_beta) + 1.0);
    double r_denominator = std::tgamma(current_alpha + 1.0) * std::tgamma(1.0 - (1.0 / current_beta)) - std::tgamma(current_alpha - (1.0 / current_beta) + 1.0);
    double r = r_numerator / r_denominator;
    
    // Ensure r is positive
    if (r <= 0.0) {
      result[i] = 0.0;
      continue;
    }
    
    double nb_coeff = std::exp(std::lgamma(r + current_x) - std::lgamma(current_x + 1) - std::lgamma(r));
    
    double remaining = 0.0;
    
    for (int j = 0; j <= current_x; ++j) {
      double binom_coeff = std::exp(std::lgamma(current_x) - std::lgamma(j + 1) - std::lgamma(current_x - j));
      double sign = (j % 2 == 0) ? 1.0 : -1.0;
      double numerator_remaining = std::tgamma(current_alpha + 1.0) * std::tgamma(1.0 + (r + j) / current_beta);
      double denominator_remaining = std::tgamma(current_alpha + (r + j) / current_beta + 1.0);
      remaining += binom_coeff * sign * numerator_remaining / denominator_remaining;
    }
    
    result[i] = nb_coeff * remaining;
  }
  
  return result;
}

// Negative Binomial Crack Distribution
// [[Rcpp::export]]
NumericVector dnbCrack_cpp(IntegerVector x, NumericVector mean, NumericVector r, NumericVector theta, NumericVector gamma) {
  int nx = x.size();
  int nmean = mean.size();
  int nr = r.size();
  int ntheta = theta.size();
  int ngamma = gamma.size();
  
  // Ensure parameters are either of length 1 or the same length as x
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
  
  NumericVector result(nx);
  
  for (int i = 0; i < nx; ++i) {
    int current_x = x[i];
    double current_mean = (nmean == 1) ? mean[0] : mean[i];
    double current_r = (nr == 1) ? r[0] : r[i];
    double current_theta = (ntheta == 1) ? theta[0] : theta[i];
    double current_gamma = (ngamma == 1) ? gamma[0] : gamma[i];
    
    // Validate input parameters
    if (current_x < 0 || current_r <= 0 || current_theta <= -0.5 || current_theta >= 0.5) {
      result[i] = NA_REAL;
      continue;
    }
    
    double discriminant = 1.0 - 2.0 * current_theta;
    
    // Check for negative discriminant
    if (discriminant <= 0.0) {
      result[i] = NA_REAL;
      continue;
    }
    
    double sqrt_term = std::sqrt(discriminant);
    double denom = 1.0 - sqrt_term;
    
    // Check for division by zero in denom
    if (denom == 0.0) {
      result[i] = NA_REAL;
      continue;
    }
    
    double gamma_term = 1.0 - current_gamma * denom;
    
    // Compute lambda using the closed-form solution
    double numerator = (current_mean + current_r) * sqrt_term;
    double denominator = current_r * gamma_term;
    if (denominator <= 0.0 || numerator <= 0.0) {
      result[i] = NA_REAL;
      continue;
    }
    
    double lambda = std::log(numerator / denominator) / denom;
    
    // Compute the negative binomial coefficient
    double nb_coeff = std::exp(std::lgamma(current_r + current_x) - std::lgamma(current_x + 1) - std::lgamma(current_r));
    
    double sum = 0.0;
    for (int j = 0; j <= current_x; ++j) {
      // Binomial coefficient
      double binom_coeff = std::exp(std::lgamma(current_x + 1) - std::lgamma(j + 1) - std::lgamma(current_x - j + 1));
      double sign = (j % 2 == 0) ? 1.0 : -1.0;
      double temp = current_r + j;
      double discriminant_j = 1.0 + 2.0 * current_theta * temp;
      
      // Check for negative discriminant
      if (discriminant_j <= 0.0) {
        result[i] = NA_REAL;
        break;
      }
      
      double sqrt_term_j = std::sqrt(discriminant_j);
      double exp_term = std::exp(lambda * (1.0 - sqrt_term_j));
      double denom_j = sqrt_term_j;
      double gamma_term_j = 1.0 - current_gamma * (1.0 - sqrt_term_j);
      
      // Check for division by zero
      if (denom_j == 0.0) {
        result[i] = NA_REAL;
        break;
      }
      
      double term = binom_coeff * sign * exp_term * gamma_term_j / denom_j;
      sum += term;
    }
    
    result[i] = nb_coeff * sum;
  }
  
  return result;
}