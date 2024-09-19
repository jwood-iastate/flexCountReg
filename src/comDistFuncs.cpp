// Functions for the Conway-Maxwell-Poisson distribution
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
# include <math.h>
using namespace Rcpp;
using Rcpp::NumericVector;


double com_adjustFactor_cpp(double lambda, double nu, double rel_tol = 1e-10, int max_iter = 1000){
  int x = 0;
  double sumval = 0;
  double iterval = 0;
  bool notconverged = true;
  
  while (notconverged) {
    iterval = std::pow(lambda,x) / std::pow(std::tgamma(x+1),nu);
    sumval += iterval;
    if (iterval<rel_tol){
      notconverged = false;
    }
    x += 1;
    if (x>=max_iter){
      Rcpp::warning("Maximum number of iterations reached in com_adjustFactor_cpp.");
      notconverged = false;
    }
  }
  
  return(sumval);
}

// derivative of the log-transformed adjustment factor
double log_adjust_deriv(double lambda, double nu, double rel_tol = 1e-10, int max_iter = 1000){
  double h = 1e-6 * lambda; // Small increment proportional to lambda
  
  // Ensure h is positive and not too small
  if (h <= 0) {
    h = 1e-8;
  }
  
  double log_adjust0 = std::log(com_adjustFactor_cpp(lambda - h, nu, rel_tol, max_iter));
  double log_adjust1 = std::log(com_adjustFactor_cpp(lambda + h, nu, rel_tol, max_iter));
  double numderiv = (log_adjust1 - log_adjust0) / (2 * h);
  
  // Check if numderiv is valid
  if (numderiv < 1e-12) {
    Rcpp::warning("Numerical derivative is too small or negative.");
    numderiv = 1e-12; // Assign a minimal positive value to avoid division by zero later
  }
  
  return(numderiv);
}

// com_expect_cpp with adaptive summation
// [[Rcpp::export]]
double com_expect_cpp(double lambda, double nu, double rel_tol = 1e-10, int max_iter = 1000) {
  
  return(lambda*log_adjust_deriv(lambda, nu, rel_tol, max_iter));
}


// Secant method implemented in Rcpp to find lambda values
// [[Rcpp::export]]
double find_lambda_cpp(double mu, double nu, double tol = 1e-8, int max_iter = 1000) {
  double lambda0 = mu * 0.5;
  double lambda1 = mu * 1.5;
  double f0 = com_expect_cpp(lambda0, nu) - mu;
  double f1 = com_expect_cpp(lambda1, nu) - mu;
  double lambda_new;
  
  for (int i = 0; i < max_iter; ++i) {
    // Check for division by zero
    if (std::abs(f1 - f0) < 1e-12) {
      Rcpp::warning("Division by zero encountered in find_lambda_cpp.");
      return lambda1;
    }
    // Update lambda using secant formula
    lambda_new = lambda1 - f1 * (lambda1 - lambda0) / (f1 - f0);
    
    // Ensure lambda_new is positive
    if (lambda_new <= 0) {
      lambda_new = (lambda1 + lambda0) / 2.0;
    }
    
    // Check for convergence
    if (std::abs(lambda_new - lambda1) < tol) {
      return lambda_new;
    }
    
    // Prepare for next iteration
    lambda0 = lambda1;
    f0 = f1;
    lambda1 = lambda_new;
    f1 = com_expect_cpp(lambda1, nu) - mu;
  }
  
  return lambda1;
}


// Find lambda for a vector
// [[Rcpp::export]]
NumericVector find_lambda_vec_cpp(NumericVector mu, NumericVector nu, double tol = 1e-8, int max_iter = 1000) {
  int n = mu.size();
  NumericVector results(n);
  
  for (int i = 0; i < n; ++i) {
    results[i] = find_lambda_cpp(mu[i], nu[i], tol, max_iter);
  }
  
  return results;
}

// Calculate the normalizing constant for the Conway-Maxwell-Poisson distribution
// [[Rcpp::export]]
NumericVector cmp_normalizer_vec_cpp(NumericVector lambda_vec, NumericVector nu_vec, double rel_tol=1e-8, int max_iter = 1000) {
  int n = lambda_vec.size();
  NumericVector results(n);
  
  for (int i = 0; i < n; ++i) {
    double lambda = lambda_vec[i];
    double nu = nu_vec[i];
    
    results[i] = com_adjustFactor_cpp(lambda, nu, rel_tol, max_iter);
  }
  
  return results;
}


// [[Rcpp::export]]
NumericVector dcom_vec_cpp(IntegerVector x_vec, NumericVector lambda_values, NumericVector nu_vec, bool log_prob = false) {
  int n = x_vec.size();
  NumericVector probabilities(n);
  
  // Compute the normalizers
  NumericVector normalizers = cmp_normalizer_vec_cpp(lambda_values, nu_vec);
  
  // Compute probabilities
  for (int i = 0; i < n; ++i) {
    double lambda = lambda_values[i];
    double nu = nu_vec[i];
    int x = x_vec[i];
    double normfactor = normalizers[i];
    
    // Compute log probability using lgamma for numerical stability
    double log_p = x * std::log(lambda) - nu * std::lgamma(x + 1) - std::log(normfactor);
    
    if (log_prob) {
      probabilities[i] = log_p;
    } else {
      probabilities[i] = std::exp(log_p);
    }
  }
  
  return probabilities;
}

double pcom_cpp(int q, double lambda, double nu) {
  double adjuster = com_adjustFactor_cpp(lambda, nu);
  double cum_prob = 0.0;
  
  for (int x = 0; x <= q; ++x) {
    double log_p = x * std::log(lambda) - nu * std::lgamma(x + 1) - std::log(adjuster);
    cum_prob += std::exp(log_p);
  }
  
  return cum_prob;
}

// double pcom_cpp(int q, double lambda, double nu) {
//   double adjuster = com_adjustFactor_cpp(lambda, nu);
//   double prob_x = 0.0;
//   
//   for (int x = 0; x <= q; ++x) {
//     double term = std::pow(lambda, x) / (std::pow(std::tgamma(x + 1), nu) * adjuster);
//     prob_x += term;
//   }
//   
//   return prob_x;
// }


// [[Rcpp::export]]
NumericVector pcom_vec_cpp(IntegerVector q_vec, NumericVector lambda_values, NumericVector nu_vec, bool lower_tail = true, bool log_p = false) {
  int n = q_vec.size();
  NumericVector cumulative_probs(n);

  
  for (int i = 0; i < n; ++i) {
    int q = q_vec[i];
    double lambda = lambda_values[i];
    double nu = nu_vec[i];
    
    cumulative_probs[i] = pcom_cpp(q, lambda, nu);
  }
  
  return cumulative_probs;
}

double qcom_cpp(double p, double lambda, double nu) {
  double adjuster = com_adjustFactor_cpp(lambda, nu);
  int lower = 0;
  int upper = 1000; // Adjust upper limit as needed
  int x = -1;
  
  // Precompute cumulative probabilities up to upper
  std::vector<double> cum_probs(upper + 1);
  double cum_prob = 0.0;
  for (int k = 0; k <= upper; ++k) {
    double log_p = k * std::log(lambda) - nu * std::lgamma(k + 1) - std::log(adjuster);
    cum_prob += std::exp(log_p);
    cum_probs[k] = cum_prob;
  }
  
  // Binary search
  while (lower <= upper) {
    int mid = (lower + upper) / 2;
    if (cum_probs[mid] < p) {
      lower = mid + 1;
    } else {
      upper = mid - 1;
      x = mid;
    }
  }
  
  if (x == -1) {
    Rcpp::warning("Quantile not found within the upper limit.");
    x = upper;
  }
  
  return x;
}

// double qcom_cpp(double p, double lambda, double nu) {
//   double adjuster = com_adjustFactor_cpp(lambda, nu);
//   double cum_prob = 0.0;
//   int x = -1;
//   
//   while (cum_prob < p) {
//     x += 1;
//     double term = std::pow(lambda, x) / (std::pow(std::tgamma(x + 1), nu) * adjuster);
//     cum_prob += term;
//     
//     if (x >= 1e6) { // Prevent infinite loop
//       Rcpp::warning("Reached maximum x in qcom_cpp.");
//       break;
//     }
//   }
//   
//   return x;
// }


// [[Rcpp::export]]
NumericVector qcom_vec_cpp(NumericVector p_vec, NumericVector lambda_values, NumericVector nu_vec) {
  int n = p_vec.size();
  NumericVector quantiles(n);
  
  for (int i = 0; i < n; ++i) {
    double p = p_vec[i];
    double lambda = lambda_values[i];
    double nu = nu_vec[i];
    
    quantiles[i] = qcom_cpp(p, lambda, nu);
  }
  
  return quantiles;
}


// [[Rcpp::export]]
NumericVector rcom_cpp(int n, double lambda, double nu) {
  NumericVector quantiles = Rcpp::runif(n);
  NumericVector lambda_vect(n, lambda);
  NumericVector nu_vect(n, nu);
  
  NumericVector random_numbers = qcom_vec_cpp(quantiles, lambda_vect, nu_vect);

  return random_numbers;
}


// 
// // Calculate the variance of the Conway-Maxwell-Poisson distribution
// // [[Rcpp::export]]
// double com_var_cpp(double lambda, double nu, int maxval = 1000) {
//   double sum = 0.0;
//   double sum_expect = 0.0;
//   double sum_expect2 = 0.0;
//   double log_lambda = std::log(lambda);
//   double log_fact = 0.0;
//   double max_log_term = -INFINITY;
//   std::vector<double> log_terms(maxval + 1);
//   
//   // Compute log terms and find the maximum for numerical stability
//   for (int k = 0; k <= maxval; ++k) {
//     if (k > 0) {
//       log_fact += std::log(k);
//     }
//     double log_term = k * log_lambda - nu * log_fact;
//     log_terms[k] = log_term;
//     if (log_term > max_log_term) {
//       max_log_term = log_term;
//     }
//   }
//   
//   // Compute the sums using the log-sum-exp trick
//   for (int k = 0; k <= maxval; ++k) {
//     double term = std::exp(log_terms[k] - max_log_term);
//     sum += term;
//     sum_expect += k * term;
//     sum_expect2 += k * k * term;
//   }
//   
//   sum *= std::exp(max_log_term);
//   sum_expect *= std::exp(max_log_term);
//   sum_expect2 *= std::exp(max_log_term);
//   
//   double mean = sum_expect / sum;
//   double mean2 = sum_expect2 / sum;
//   return mean2 - mean * mean; // Variance = E[X^2] - (E[X])^2
// }




