#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

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
    stop("mean must be of length 1 or the same length as x");
  }
  if (nsigma != 1 && nsigma != nx) {
    stop("sigma must be of length 1 or the same length as x");
  }
  
  NumericVector result(nx);
  
  for (int i = 0; i < nx; ++i) {
    double current_mean = (nmean == 1) ? mean[0] : mean[i];
    double current_sigma = (nsigma == 1) ? sigma[0] : sigma[i];
    
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