#include <Rcpp.h>
#include <cmath>
using Rcpp::NumericVector;

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
