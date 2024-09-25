#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector p_poisweibull_cpp(NumericVector beta,
                                NumericVector x,
                                NumericMatrix X_Fixed,
                                NumericMatrix X_alpha,  
                                NumericMatrix X_sigma, 
                                NumericVector h) {
  int nx = x.size();
  int N_fixed = X_Fixed.ncol(); // Number of fixed effect parameters
  int Nrows = X_Fixed.nrow(); // Number of rows in X_Fixed
  int N_alpha = X_alpha.ncol(); // Number of alpha parameters
  int N_sigma = X_sigma.ncol(); // Number of sigma parameters
  int nh = h.size();
  
  int pars = beta.size(); // Total number of coefficients
  NumericVector coefs = beta[Range(0, N_fixed - 1)];
  NumericVector alphas = beta[Range(N_fixed, N_fixed + N_alpha - 1)];
  NumericVector sigmas = beta[Range(N_fixed + N_alpha, pars - 1)];
  
  // Compute mu, alpha, and sigma
  NumericVector mu(Nrows);
  NumericVector alpha(Nrows);
  NumericVector sigma(Nrows);
  for (int i = 0; i < Nrows; ++i) {
    double dot_product = 0.0;
    double dot_product_alpha = 0.0;
    double dot_product_sigma = 0.0;
    for (int j = 0; j < N_fixed; ++j) {
      dot_product += X_Fixed(i, j) * coefs[j];
    }
    for (int j = 0; j < N_alpha; ++j) {
      dot_product_alpha += X_alpha(i, j) * alphas[j];
    }
    for (int j = 0; j < N_sigma; ++j) {
      dot_product_sigma += X_sigma(i, j) * sigmas[j];
    }
    mu[i] = std::exp(dot_product);
    alpha[i] = std::exp(dot_product_alpha);
    sigma[i] = std::exp(dot_product_sigma);
  }
  
  NumericVector result(nx);
  
  for (int i = 0; i < nx; ++i) {
    double current_alpha = alpha[i];
    double current_sigma = sigma[i];
    double p = 0.0;
    
    // Calculate lambda
    double lambda = mu[i] / (current_sigma * std::tgamma(1 + 1 / current_alpha));
    
    for (int j = 0; j < nh; ++j) {
      double current_h = h[j];
      double mu_j = lambda * R::qweibull(current_h, current_alpha, current_sigma, true, false);
      p += R::dpois(x[i], mu_j, false);
    }
    
    result[i] = p / nh;
  }
  
  return log(result);
}
