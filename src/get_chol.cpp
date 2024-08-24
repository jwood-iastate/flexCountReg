#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix get_chol_cpp(NumericVector pars, int Nvars) {
  NumericMatrix Ch(Nvars, Nvars);
  int counter = 0;
  
  for (int i = 0; i < Nvars; ++i) {
    for (int j = 0; j <= i; ++j) {
      Ch(j, i) = pars[counter];
      counter++;
    }
  }
  
  return Ch;
}