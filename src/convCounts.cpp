// Log probabilities of count distributions that are sums of independent
// count variables, computed by direct summation of positive terms on the log
// scale (no subtraction, so tiny tail probabilities keep their relative
// precision).
//
//   Delaporte        Y = A + B,       A ~ Poisson(lam), B ~ NB(size, mean nbmu)
//   Hermite (order m) Y = N1 + m * Nm, N1 ~ Poisson(lam1), Nm ~ Poisson(lamm)
//   Polya-Aeppli     Y = C_1 + ... + C_M, M ~ Poisson(lam), the C_i i.i.d.
//                    geometric on 1, 2, ... with P(C = c) = (1 - q) q^(c - 1);
//                    given M = j >= 1, Y - j ~ NB(size = j, prob = 1 - q).
//
// For each count the sum runs over the value j of the Poisson component. The
// terms are products of a Poisson probability and a probability (or tail
// probability) that is at most 1, so once j is past the Poisson mode and the
// Poisson log probability is 60 below the running total, the remaining terms
// (bounded by the Poisson upper tail) cannot change the result in double
// precision and the sum stops.
//
// The tails use the same decomposition:
//   P(Y <= y) = sum_j P(J = j) P(rest <= y - c j)
//   P(Y >  y) = sum_j P(J = j) P(rest >  y - c j) + P(J > floor(y / c))
// with c the multiplier of the Poisson component (1, m, or 1 for the
// Polya-Aeppli, where every jump is at least 1).
//
// Inputs are validated (and recycled) in R: y holds non-negative whole
// numbers and the parameters are inside their domains.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <limits>
using namespace Rcpp;

namespace {

const double NEG_INF = -std::numeric_limits<double>::infinity();

// Running log-sum-exp
struct LogSum {
  double m = NEG_INF;
  double s = 0.0;
  void add(double t) {
    if (t == NEG_INF || std::isnan(t)) return;
    if (t > m) {
      s = s * std::exp(m - t) + 1.0;
      m = t;
    } else {
      s += std::exp(t - m);
    }
  }
  double value() const { return m == NEG_INF ? NEG_INF : m + std::log(s); }
};

// Can the sum over j stop? (past the Poisson mode and negligible)
inline bool negligible(double lp_pois, double j, double lam,
                       const LogSum& acc) {
  return j > lam && acc.m != NEG_INF && lp_pois < acc.value() - 60.0;
}

// 0 = log pmf, 1 = log lower tail P(Y <= y), 2 = log upper tail P(Y > y)
double delap_one(double y, double lam, double size, double nbmu, int what) {
  LogSum acc;
  for (double j = 0; j <= y; j += 1.0) {
    double lp = R::dpois(j, lam, 1);
    if (lp == NEG_INF) {
      if (j > lam) break;
      continue;
    }
    double r = y - j, lb;
    if (what == 0) lb = R::dnbinom_mu(r, size, nbmu, 1);
    else lb = R::pnbinom_mu(r, size, nbmu, what == 1, 1);
    acc.add(lp + lb);
    if (negligible(lp, j, lam, acc)) break;
  }
  if (what == 2) acc.add(R::ppois(y, lam, 0, 1));
  return acc.value();
}

double herm_one(double y, double lam1, double lamm, double m, int what) {
  LogSum acc;
  double jmax = std::floor(y / m + 1e-9);
  for (double j = 0; j <= jmax; j += 1.0) {
    double lp = R::dpois(j, lamm, 1);
    if (lp == NEG_INF) {
      if (j > lamm) break;
      continue;
    }
    double r = y - m * j, lb;
    if (what == 0) lb = R::dpois(r, lam1, 1);
    else lb = R::ppois(r, lam1, what == 1, 1);
    acc.add(lp + lb);
    if (negligible(lp, j, lamm, acc)) break;
  }
  if (what == 2) acc.add(R::ppois(jmax, lamm, 0, 1));
  return acc.value();
}

double pa_one(double y, double lam, double q, int what) {
  LogSum acc;
  // j = 0 (no clusters): Y = 0
  double lp0 = -lam;
  if (what == 0) {
    if (y == 0) acc.add(lp0);
  } else if (what == 1) {
    acc.add(lp0);
  }
  for (double j = 1; j <= y; j += 1.0) {
    double lp = R::dpois(j, lam, 1);
    if (lp == NEG_INF) {
      if (j > lam) break;
      continue;
    }
    double r = y - j, lb;
    if (what == 0) lb = R::dnbinom(r, j, 1.0 - q, 1);
    else lb = R::pnbinom(r, j, 1.0 - q, what == 1, 1);
    acc.add(lp + lb);
    if (negligible(lp, j, lam, acc)) break;
  }
  // More than y clusters means Y > y (every cluster has at least one event)
  if (what == 2) acc.add(R::ppois(y, lam, 0, 1));
  return acc.value();
}

}  // namespace

// [[Rcpp::export]]
NumericVector delap_log_cpp(NumericVector y, NumericVector lam,
                            NumericVector size, NumericVector nbmu,
                            int what) {
  R_xlen_t n = y.size();
  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i % 1024 == 0) Rcpp::checkUserInterrupt();
    out[i] = delap_one(y[i], lam[i], size[i], nbmu[i], what);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector herm_log_cpp(NumericVector y, NumericVector lam1,
                           NumericVector lamm, NumericVector m, int what) {
  R_xlen_t n = y.size();
  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i % 1024 == 0) Rcpp::checkUserInterrupt();
    out[i] = herm_one(y[i], lam1[i], lamm[i], m[i], what);
  }
  return out;
}

// [[Rcpp::export]]
NumericVector pa_log_cpp(NumericVector y, NumericVector lam, NumericVector q,
                         int what) {
  R_xlen_t n = y.size();
  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i % 1024 == 0) Rcpp::checkUserInterrupt();
    out[i] = pa_one(y[i], lam[i], q[i], what);
  }
  return out;
}

// Cumulative log-sum-exp: out[k] = log(sum(exp(x[0..k]))).
// [[Rcpp::export]]
NumericVector cumlogsumexp_cpp(NumericVector x) {
  R_xlen_t n = x.size();
  NumericVector out(n);
  LogSum acc;
  for (R_xlen_t i = 0; i < n; ++i) {
    acc.add(x[i]);
    out[i] = acc.value();
  }
  return out;
}
