#include <Rcpp.h>
#include <cmath>

// [[Rcpp::export]]
double dtri_cpp(double x, double mode = 0, double sigma = 1, double upper = NA_REAL, double lower = NA_REAL, bool log = false) {
  double a, b, c = mode;
  
  if (!Rcpp::NumericVector::is_na(upper) && !Rcpp::NumericVector::is_na(lower)) {
    if (lower >= upper) {
      Rcpp::stop("The value of `lower` must be smaller than the value of `upper`");
    } else if (lower > c) {
      Rcpp::stop("The value of `mode` must be greater than or equal to the value of `lower`");
    } else if (c > upper) {
      Rcpp::stop("The value of `mode` must be smaller than or equal to the value of `upper`");
    }
    a = lower;
    b = upper;
  } else {
    a = c - sigma;
    b = c + sigma;
  }
  
  double p;
  if (x <= a || x >= b) {
    p = 0;
  } else if (x < c) {
    p = 2 * (x - a) / ((b - a) * (c - a));
  } else if (x == c) {
    p = 2 / (b - a);
  } else {
    p = 2 * (b - x) / ((b - a) * (b - c));
  }
  
  if (log) return std::log(p);
  else return p;
}

// [[Rcpp::export]]
double ptri_cpp(double q, double mode = 0, double sigma = 1, double upper = NA_REAL, double lower = NA_REAL, bool lower_tail = true, bool log_p = false) {
  double a, b, c = mode;
  
  if (!Rcpp::NumericVector::is_na(upper) && !Rcpp::NumericVector::is_na(lower)) {
    if (lower >= upper) {
      Rcpp::stop("The value of `lower` must be smaller than the value of `upper`");
    } else if (lower > c) {
      Rcpp::stop("The value of `mode` must be greater than or equal to the value of `lower`");
    } else if (c > upper) {
      Rcpp::stop("The value of `mode` must be smaller than or equal to the value of `upper`");
    }
    a = lower;
    b = upper;
  } else {
    a = c - sigma;
    b = c + sigma;
  }
  
  double p;
  if (q < a) {
    p = 0;
  } else if (q > b) {
    p = 1;
  } else if (q <= c) {
    p = std::pow(q - a, 2) / ((b - a) * (c - a));
  } else {
    p = 1 - std::pow(b - q, 2) / ((b - a) * (b - c));
  }
  
  if (!lower_tail) p = 1 - p;
  if (log_p) return std::log(p);
  else return p;
}

// [[Rcpp::export]]
double qtri_cpp(double p, double mode = 0, double sigma = 1, double upper = NA_REAL, double lower = NA_REAL) {
  double a, b, c = mode;
  double p_mode;
  
  if (!Rcpp::NumericVector::is_na(upper) && !Rcpp::NumericVector::is_na(lower)) {
    if (lower >= upper) {
      Rcpp::stop("The value of `lower` must be smaller than the value of `upper`");
    } else if (lower > c) {
      Rcpp::stop("The value of `mode` must be greater than or equal to the value of `lower`");
    } else if (c > upper) {
      Rcpp::stop("The value of `mode` must be smaller than or equal to the value of `upper`");
    }
    a = lower;
    b = upper;
    p_mode = ptri_cpp(mode, mode, sigma, upper, lower);
  } else {
    a = c - sigma;
    b = c + sigma;
    p_mode = ptri_cpp(mode, mode, sigma);
  }
  
  double q;
  if (p <= p_mode) {
    q = a + std::sqrt(p * (b - a) * (c - a));
  } else {
    q = b - std::sqrt((1 - p) * (b - a) * (b - c));
  }
  return q;
}

// [[Rcpp::export]]
Rcpp::NumericVector rtri_cpp(int n, double mode = 0, double sigma = 1, double upper = NA_REAL, double lower = NA_REAL) {
  Rcpp::NumericVector result(n);
  Rcpp::NumericVector u = Rcpp::runif(n);
  
  for (int i = 0; i < n; ++i) {
    result[i] = qtri_cpp(u[i], mode, sigma, upper, lower);
  }
  
  return result;
}
