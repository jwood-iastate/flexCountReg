#include <Rcpp.h>
#include <cmath> // Include the cmath header

using namespace Rcpp;
using namespace std; // Optional
using Rcpp::NumericVector;

// [[Rcpp::export]]
NumericVector dplindgamma_cpp(NumericVector x, NumericVector mean, NumericVector theta, NumericVector alpha, NumericVector h) {
  int nx = x.size();
  int nmean = mean.size();
  int nalpha = alpha.size();
  int ntheta = theta.size();
  int nh = h.size();
  
  // Ensure mean and sigma are either of length 1 or the same length as x
  if (nmean != 1 && nmean != nx) {
    Rcpp::stop("mean must be of length 1 or the same length as x");
  }
  if (nalpha != 1 && nalpha != nx) {
    Rcpp::stop("alpha must be of length 1 or the same length as x");
  }
  if (ntheta != 1 && ntheta != nx) {
    Rcpp::stop("theta must be of length 1 or the same length as x");
  }
  
  NumericVector result(nx);
  
  for (int i = 0; i < nx; ++i) {
    double current_mean = (nmean == 1) ? mean[0] : mean[i];
    double current_alpha = (nalpha == 1) ? alpha[0] : alpha[i];
    double current_theta = (ntheta == 1) ? theta[0] : theta[i];
    double p = 0.0;
    
    if (x[i] < 0 || current_mean <= 0 || current_alpha <= 0 || current_theta <= 0) {
      result[i] = 0.0;
      continue;
    }
    
    for (int j = 0; j < nh; ++j) {
      double current_h = h[j];
      double mu_j = current_mean * R::qgamma(current_h, current_alpha, 1.0 / current_alpha, 1, 0);
      double lambda = mu_j * current_theta * (current_theta + 1) / (current_theta + 2);
      
      p += (std::pow(current_theta, 2) * std::pow(lambda, x(i)) * (current_theta + lambda + x(i) + 1) / ((current_theta + 1) * std::pow(current_theta + lambda, x(i)+ 2)));
    }
    
    result[i] = p / nh;
  }
  
  return result;
}

