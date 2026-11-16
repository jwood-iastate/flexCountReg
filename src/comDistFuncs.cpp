// Conway-Maxwell-Poisson (COM) distribution:
//   P(Y = y) = lambda^y / (y!)^nu / Z(lambda, nu),  y = 0, 1, 2, ...
// Everything is computed on the log scale. All functions use the same
// support rule: terms are summed from y = 0 until y is past the mode
// lambda^(1/nu) and the next term is negligible (relative 1e-17) compared
// with the running sum, so the normalizing constant, moments, CDF and
// quantiles are consistent with each other.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

static const long COM_MAX_TERMS = 10000000L;
static const double COM_NEGLIGIBLE = -39.14; // log(1e-17)

static inline double com_logterm(double y, double loglam, double nu) {
  return y * loglam - nu * std::lgamma(y + 1.0);
}

// log(exp(a) + exp(b))
static inline double com_logadd(double a, double b) {
  if (a == R_NegInf) return b;
  if (b == R_NegInf) return a;
  double m = std::max(a, b);
  return m + std::log1p(std::exp(std::min(a, b) - m));
}

static inline bool com_valid(double lambda, double nu) {
  return !ISNAN(lambda) && !ISNAN(nu) && lambda > 0 && nu > 0 &&
    R_FINITE(lambda) && R_FINITE(nu);
}

// Log of the sum of the unnormalized terms for y in [from, to] (to < 0: to
// the end of the support).
static double com_log_sum(long from, long to, double lambda, double nu) {
  double loglam = std::log(lambda);
  double mode = std::pow(lambda, 1.0 / nu);
  double total = R_NegInf;
  for (long y = from; to < 0 || y <= to; ++y) {
    double lt = com_logterm((double) y, loglam, nu);
    total = com_logadd(total, lt);
    if (to < 0 && y > mode && lt - total < COM_NEGLIGIBLE) return total;
    if (y - from > COM_MAX_TERMS) return R_NaN;  // reported by callers
    if ((y & 0xFFFF) == 0) Rcpp::checkUserInterrupt();
  }
  return total;
}

// log Z(lambda, nu)
static double com_log_z(double lambda, double nu) {
  if (!com_valid(lambda, nu)) return R_NaN;
  return com_log_sum(0, -1, lambda, nu);
}

// Mean and variance from the normalized log weights
static void com_moments(double lambda, double nu, double &mean,
                        double &var) {
  mean = R_NaN;
  var = R_NaN;
  double log_z = com_log_z(lambda, nu);
  if (ISNAN(log_z)) return;
  double loglam = std::log(lambda);
  double mode = std::pow(lambda, 1.0 / nu);
  double m1 = 0.0, m2 = 0.0;
  for (long y = 0; y <= COM_MAX_TERMS; ++y) {
    double lt = com_logterm((double) y, loglam, nu) - log_z;
    double p = std::exp(lt);
    m1 += y * p;
    m2 += (double) y * y * p;
    if (y > mode && lt < COM_NEGLIGIBLE && y * p < 1e-17 * (m1 + 1e-300))
      break;
  }
  mean = m1;
  var = m2 - m1 * m1;
}

// [[Rcpp::export]]
double com_expect_cpp(double lambda, double nu, double rel_tol = 1e-10,
                      int max_iter = 1000) {
  double m, v;
  com_moments(lambda, nu, m, v);
  return m;
}

// Mean and variance for vectors of parameters: a matrix with columns mean
// and variance.
// [[Rcpp::export]]
NumericMatrix com_moments_vec_cpp(NumericVector lambda, NumericVector nu) {
  int n = lambda.size();
  NumericMatrix out(n, 2);
  for (int i = 0; i < n; ++i) {
    double m, v;
    com_moments(lambda[i], nu[i], m, v);
    out(i, 0) = m;
    out(i, 1) = v;
  }
  return out;
}

// Rate lambda with E[Y] = mu, solved in log(lambda) (the mean increases
// with lambda) by a bracketed regula falsi (Illinois) method. NaN, with a
// warning, when the solution is not found.
// [[Rcpp::export]]
double find_lambda_cpp(double mu, double nu, double tol = 1e-8,
                       int max_iter = 1000) {
  if (ISNAN(mu) || ISNAN(nu) || !(mu > 0) || !(nu > 0) || !R_FINITE(mu) ||
      !R_FINITE(nu)) {
    return R_NaN;
  }
  // A mean that cannot be computed (support too long) is treated as too
  // large, so the bracket moves to smaller rates
  auto f = [&](double t) {
    double m = com_expect_cpp(std::exp(t), nu) - mu;
    return ISNAN(m) ? R_PosInf : m;
  };
  // Starting guess from the approximation E[Y] ~ lambda^(1/nu) - (nu-1)/(2nu)
  double g = mu + (nu - 1.0) / (2.0 * nu);
  double t0 = nu * std::log(g > 0 ? g : mu);
  // Narrow first steps: d E[Y] / d log(lambda) = Var(Y) can be very large
  double step = 0.1 * std::min(1.0, nu);
  double lo = t0 - step, hi = t0 + step;
  double flo = f(lo), fhi = f(hi);
  for (int k = 0; k < 200 && !(flo < 0); ++k) {
    hi = lo; fhi = flo;
    lo -= step * std::pow(2.0, k + 1);
    flo = f(lo);
  }
  for (int k = 0; k < 200 && !(fhi > 0); ++k) {
    lo = hi; flo = fhi;
    hi += step * std::pow(2.0, k + 1);
    fhi = f(hi);
  }
  if (fhi == R_PosInf) {
    // Shrink towards lo until the mean can be computed
    for (int k = 0; k < 200 && fhi == R_PosInf; ++k) {
      hi = (lo + hi) / 2;
      fhi = f(hi);
    }
    if (!(fhi > 0)) { lo = hi; flo = fhi; hi = lo; }
  }
  if (!(flo < 0 && fhi > 0)) {
    Rcpp::warning("find_lambda_cpp: could not bracket the COM rate; NaN "
                  "returned.");
    return R_NaN;
  }
  int side = 0;
  double tol_f = 1e-10 * std::max(1.0, mu);
  for (int it = 0; it < 500; ++it) {
    double t = (lo * fhi - hi * flo) / (fhi - flo);
    double ft = f(t);
    if (ISNAN(ft)) break;
    if (std::fabs(ft) < tol_f || (hi - lo) < 1e-13 * std::max(1.0,
                                                               std::fabs(t)))
      return std::exp(t);
    if (ft < 0) {
      lo = t; flo = ft;
      if (side == -1) fhi /= 2;
      side = -1;
    } else {
      hi = t; fhi = ft;
      if (side == 1) flo /= 2;
      side = 1;
    }
  }
  Rcpp::warning("find_lambda_cpp: the COM rate did not converge; NaN "
                "returned.");
  return R_NaN;
}

// [[Rcpp::export]]
NumericVector find_lambda_vec_cpp(NumericVector mu, NumericVector nu,
                                  double tol = 1e-8, int max_iter = 1000) {
  int n = mu.size();
  NumericVector results(n);
  for (int i = 0; i < n; ++i) {
    results[i] = find_lambda_cpp(mu[i], nu[i], tol, max_iter);
  }
  return results;
}

// Normalizing constants Z (kept for compatibility; may overflow to Inf, so
// the other functions use log Z).
// [[Rcpp::export]]
NumericVector cmp_normalizer_vec_cpp(NumericVector lambda_vec,
                                     NumericVector nu_vec,
                                     double rel_tol = 1e-8,
                                     int max_iter = 1000) {
  int n = lambda_vec.size();
  NumericVector results(n);
  for (int i = 0; i < n; ++i) {
    results[i] = std::exp(com_log_z(lambda_vec[i], nu_vec[i]));
  }
  return results;
}

// [[Rcpp::export]]
NumericVector dcom_vec_cpp(IntegerVector x_vec, NumericVector lambda_values,
                           NumericVector nu_vec, bool log_prob = false) {
  int n = x_vec.size();
  NumericVector out(n);
  double prev_lambda = R_NaN, prev_nu = R_NaN, log_z = R_NaN;
  for (int i = 0; i < n; ++i) {
    double lambda = lambda_values[i];
    double nu = nu_vec[i];
    if (!com_valid(lambda, nu)) {
      out[i] = R_NaN;
      continue;
    }
    if (lambda != prev_lambda || nu != prev_nu) {
      log_z = com_log_z(lambda, nu);
      prev_lambda = lambda;
      prev_nu = nu;
    }
    double lp;
    if (x_vec[i] == NA_INTEGER) {
      out[i] = NA_REAL;
      continue;
    } else if (x_vec[i] < 0) {
      lp = R_NegInf;
    } else {
      lp = com_logterm((double) x_vec[i], std::log(lambda), nu) - log_z;
    }
    out[i] = log_prob ? lp : std::exp(lp);
  }
  return out;
}

// log P(Y <= q) (lower) or log P(Y > q) (upper), each tail summed directly
// so that small tails keep their relative precision.
double com_log_tail_cpp(int q, double lambda, double nu, bool lower_tail) {
  if (!com_valid(lambda, nu)) return R_NaN;
  if (q < 0) return lower_tail ? R_NegInf : 0.0;
  double log_z = com_log_z(lambda, nu);
  double log_t = lower_tail ? com_log_sum(0, q, lambda, nu)
                            : com_log_sum((long) q + 1, -1, lambda, nu);
  return std::min(0.0, log_t - log_z);
}

// [[Rcpp::export]]
NumericVector pcom_vec_cpp(IntegerVector q_vec, NumericVector lambda_values,
                           NumericVector nu_vec, bool lower_tail = true,
                           bool log_p = false) {
  int n = q_vec.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    if (q_vec[i] == NA_INTEGER) {
      out[i] = NA_REAL;
      continue;
    }
    double lp = com_log_tail_cpp(q_vec[i], lambda_values[i], nu_vec[i],
                                 lower_tail);
    out[i] = log_p ? lp : std::exp(lp);
  }
  return out;
}

// Smallest y with P(Y <= y) >= p, using the same support rule as the
// normalizing constant.
static double qcom_cpp(double p, double lambda, double nu) {
  if (ISNAN(p)) return NA_REAL;
  if (!com_valid(lambda, nu) || p < 0 || p > 1) return R_NaN;
  if (p == 1.0) return R_PosInf;
  double log_z = com_log_z(lambda, nu);
  if (ISNAN(log_z)) return R_NaN;
  double loglam = std::log(lambda);
  double mode = std::pow(lambda, 1.0 / nu);
  double cum = 0.0;
  for (long y = 0; y <= COM_MAX_TERMS; ++y) {
    double lt = com_logterm((double) y, loglam, nu) - log_z;
    cum += std::exp(lt);
    if (cum >= p) return (double) y;
    // Past the end of the support (p within rounding of 1)
    if (y > mode && lt < COM_NEGLIGIBLE) return (double) y;
    if ((y & 0xFFFF) == 0) Rcpp::checkUserInterrupt();
  }
  return R_NaN;
}

// [[Rcpp::export]]
NumericVector qcom_vec_cpp(NumericVector p_vec, NumericVector lambda_values,
                           NumericVector nu_vec) {
  int n = p_vec.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    out[i] = qcom_cpp(p_vec[i], lambda_values[i], nu_vec[i]);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector rcom_cpp(int n, double lambda, double nu) {
  NumericVector u = Rcpp::runif(n);
  NumericVector out(n);
  for (int i = 0; i < n; ++i) out[i] = qcom_cpp(u[i], lambda, nu);
  return out;
}
