#include <Rcpp.h>
#include <cmath>

using Rcpp::NumericVector;

// PMF of the generalized Waring distribution
// [[Rcpp::export]]
NumericVector genWaring_cpp(NumericVector x, NumericVector mean, NumericVector k, NumericVector p, bool log_prob = false) {
  int nx = x.size();
  int nmean = mean.size();
  int nk = k.size();
  int np = p.size();
  
  if (nmean != 1 && nmean != nx) {
    Rcpp::stop("mean must be of length 1 or the same length as x");
  }
  if (nk != 1 && nk != nx) {
    Rcpp::stop("k must be of length 1 or the same length as x");
  }
  if (np != 1 && np != nx) {
    Rcpp::stop("p must be of length 1 or the same length as x");
  }
  
  NumericVector result(nx);
  
  for (int i = 0; i < nx; ++i) {
    double current_mean = (nmean == 1) ? mean[0] : mean[i];
    double current_k = (nk == 1) ? k[0] : k[i];
    double current_p = (np == 1) ? p[0] : p[i];
    
    // Update conditions to account for invalid parameters
    // or likelihood
    if (!R_finite(current_mean) || current_mean <= 0 ||
        !R_finite(current_k) || current_k <= 0 ||
        !R_finite(current_p) || current_p <= 1) {
        result[i] = R_NaN;
      continue;
    }
    if (Rcpp::NumericVector::is_na(x[i])) {
      result[i] = NA_REAL;
      continue;
    }
    if (!R_finite(x[i]) || x[i] < 0 || x[i] != std::floor(x[i])) {
      result[i] = log_prob ? R_NegInf : 0.0;
      continue;
    }
    
    double a = current_mean * (current_p - 1.0) / current_k;
    
    if (a <= 0) {
      result[i] = 0.0;
      continue;
    }
    
    double ay = a + x[i];
    double ky = current_k + x[i];
    double pk = current_k + current_p;
    double ap = a + current_p;
    double akpy = ay + pk;
    
    double log_numerator = std::lgamma(ay) + std::lgamma(ky) + 
      std::lgamma(pk) + std::lgamma(ap);
    double log_denominator = std::lgamma(a) +
      std::lgamma(current_k) + std::lgamma(current_p) +
      std::lgamma(akpy) + std::lgamma(x[i] + 1.0);
    
    double lp = log_numerator - log_denominator;
    result[i] = log_prob ? lp : std::exp(lp);
  }
  
  return result;
}
