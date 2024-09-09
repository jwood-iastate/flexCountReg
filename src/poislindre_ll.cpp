#include <Rcpp.h>
using namespace Rcpp;

// Function to convert a SEXP object to an IntegerVector of group indices
IntegerVector convert_to_integer_group(SEXP group) {
  if (TYPEOF(group) == INTSXP) {
    return as<IntegerVector>(group);
  } else if (TYPEOF(group) == STRSXP) {
    CharacterVector group_char = as<CharacterVector>(group);
    return match(group_char, sort_unique(group_char)) - 1;
  } else if (Rf_isFactor(group)) {
    IntegerVector group_factor = as<IntegerVector>(group);
    return group_factor - 1; // Convert factor levels to 0-based indices
  } else {
    stop("Group variable must be of type integer, character, or factor.");
  }
}

// [[Rcpp::export]]
double pollind_i_group(NumericVector mu, IntegerVector y, double theta) {
  
  // Validate inputs: check if vectors 'mu' and 'y' have the same length
  if (mu.size() != y.size()) {
    stop("Vectors 'mu' and 'y' must have the same length.");
  }
  
  // Number of observations in the group
  int n_i = y.size();
  
  // Compute constants for the probability formula
  double theta_factor = (theta * (theta + 1)) / (theta + 2);
  
  // Compute the first part of the probability
  double prob_part1 = (theta * theta) / (theta + 1);
  
  // Compute the second part: the product term
  double prob_product = 1.0;
  for (int i = 0; i < n_i; ++i) {
    prob_product *= pow(mu[i] * theta_factor, y[i]) / R::gammafn(y[i] + 1); // Using gamma function for factorial
  }
  
  // Sum terms for the third part
  int sum_y = std::accumulate(y.begin(), y.end(), 0);
  double sum_mu = std::inner_product(mu.begin(), mu.end(), NumericVector(n_i, theta_factor).begin(), 0.0);
  
  // Compute the third part of the probability
  double prob_part3 = R::gammafn(sum_y + 1) * (sum_mu + theta + sum_y + 1) / pow(sum_mu + theta, sum_y + 2);
  
  // Combine all parts to get the total probability
  double probability = prob_part1 * prob_product * prob_part3;
  
  return probability;
}

// [[Rcpp::export]]
double reg_run_RE(NumericVector beta, IntegerVector y, NumericMatrix X, SEXP group, NumericVector weights) {
  
  int pars = beta.size() - 1; // Number of coefficients (excluding theta)
  NumericVector coefs = beta[Range(0, pars - 1)];
  double theta = exp(beta[pars]); // Exponentiate the last element of beta to get theta
  
  // Compute mu as exp(X %*% coefs)
  NumericVector mu(X.nrow());
  for (int i = 0; i < X.nrow(); ++i) {
    double dot_product = 0.0;
    for (int j = 0; j < X.ncol(); ++j) {
      dot_product += X(i, j) * coefs[j];
    }
    mu[i] = exp(dot_product);
  }
  
  // Convert the group variable to an integer vector of indices
  IntegerVector group_indices = convert_to_integer_group(group);
  
  // Compute log-likelihood contributions per group
  IntegerVector unique_groups = sort_unique(group_indices);
  int n_groups = unique_groups.size();
  
  NumericVector log_probs(n_groups);
  NumericVector group_weights(n_groups);
  NumericVector group_n(n_groups);
  
  for (int g = 0; g < n_groups; ++g) {
    int current_group = unique_groups[g];
    group_n[g] = std::count(group_indices.begin(), group_indices.end(), current_group);
    
    // Extract y and mu for the current group
    std::vector<int> y_group;
    std::vector<double> mu_group;
    for (int i = 0; i < group_indices.size(); ++i) {
      if (group_indices[i] == current_group) {
        y_group.push_back(y[i]);
        mu_group.push_back(mu[i]);
      }
    }
    
    // Convert std::vector to Rcpp vectors
    IntegerVector y_group_rcpp = wrap(y_group); // Ensure correct type for y (IntegerVector)
    NumericVector mu_group_rcpp = wrap(mu_group); // Ensure correct type for mu (NumericVector)
    
    // Calculate the probability for the group
    double prob = pollind_i_group(mu_group_rcpp, y_group_rcpp, theta);
    log_probs[g] = log(prob);
    
    // Calculate the sum of weights for the group
    double sum_weights = 0.0;
    for (int i = 0; i < group_indices.size(); ++i) {
      if (group_indices[i] == current_group) {
        sum_weights += weights[i];
      }
    }
    group_weights[g] = sum_weights;
  }
  
  // Compute the weighted sum of log-likelihood contributions
  double ll = 0.0;
  for (int g = 0; g < n_groups; ++g) {
    ll += log_probs[g] * group_weights[g] / group_n[g];
  }
  
  return ll;
}
