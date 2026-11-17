// Mean parameterizations of count distributions whose mean has no closed
// form, used by the regression families (countreg(), countreg.rp()): given
// the mean mu and the shape/dispersion parameter, the natural parameter is
// found by a safeguarded Newton iteration (with a bisection bracket) on its
// logarithm, using the analytic derivative of the mean.
//
//   DPO        double Poisson, location lambda:  d E / d log(lambda) = phi Var
//   HYPERPOIS  hyper-Poisson, rate lambda:       d E / d log(lambda) = Var
//   DW         discrete Weibull, c = -log(theta): E = sum_k exp(-c k^b)
//   GAMMACOUNT gamma renewal count, unit window, rate s:
//              E = sum_k P(Gamma(k a, 1) <= s),  d E / d s = sum_k density
//
// DPO and hyper-Poisson probabilities are normalized by summing the weights
// over a window around the mode, grown until the weights at its ends are 50
// log units below the largest one and decrease outward with non-increasing
// ratios (the omitted tails are then below geometric series, which are
// added); the work is proportional to the standard deviation, not the mean. The
// discrete Weibull sums are taken directly up to a cut-off and completed by
// the Euler-Maclaurin formula with closed-form incomplete gamma integrals.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

namespace {

const double NEG_INF = -std::numeric_limits<double>::infinity();

double logaddexp(double a, double b) {
  if (a == NEG_INF) return b;
  if (b == NEG_INF) return a;
  double m = std::max(a, b);
  return m + std::log1p(std::exp(-std::fabs(a - b)));
}

// ---- weight tables (DPO = 0, hyper-Poisson = 1) ----------------------------

inline double tab_logw(int fam, double y, double loglam, double par,
                       double lgpar) {
  if (fam == 0) {
    double ylogy = (y > 0) ? y * std::log(y) : 0.0;
    return par * y * loglam + (1.0 - par) * ylogy + (par - 1.0) * y -
      std::lgamma(y + 1.0);
  }
  return y * loglam - std::lgamma(par + y) + lgpar;
}

// log(1 - exp(x)) for x <= 0
double log1mexp(double x) {
  return (x > -M_LN2) ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x));
}

// Normalized weights on a window [L0, L0 + size - 1] around the mode. The
// window grows until the weights at both ends are 50 log units below the
// largest one and decrease outward with non-increasing ratios; the omitted
// tails on each side are then below geometric series, which are added to
// the normalizing constant. Work is O(sd), not O(mean).
// log w(y + 1) - log w(y), used to fill tables by recurrence (one log per
// term instead of lgamma); y log(1 + 1 / y) is 0 at y = 0
inline double tab_dlogw(int fam, double y, double loglam, double par) {
  if (fam == 0) {
    double ylog = (y > 0) ? y * std::log1p(1.0 / y) : 0.0;
    return par * (loglam - std::log(y + 1.0)) + (1.0 - par) * ylog +
      par - 1.0;
  }
  return loglam - std::log(par + y);
}

struct Table {
  std::vector<double> lw;
  double L0 = 0, logZ = NA_REAL, mean = NA_REAL, var = NA_REAL;
  double ltail_l = NEG_INF, ltail_r = NEG_INF;
  int fam = 0;
  double loglam = 0, par = 1, lgpar = 0;
  bool ok = false;
  double logw(double y) const { return tab_logw(fam, y, loglam, par, lgpar); }
};

void tab_fill(Table &t, int fam, double loglam, double par) {
  t.ok = false;
  t.fam = fam; t.loglam = loglam; t.par = par;
  t.lgpar = (fam == 1) ? std::lgamma(par) : 0.0;
  double lam = std::exp(loglam);
  double m = (fam == 0) ? lam : std::max(lam - par + 1.0, 0.0);
  double sd = (fam == 0) ? std::sqrt(lam / par + 1.0) : std::sqrt(lam + 1.0);
  if (!R_FINITE(m) || !R_FINITE(sd)) return;
  double half = 10.0 * sd + 30.0;
  double L = std::max(0.0, std::floor(m - half)), R = std::ceil(m + half);
  // Work budget: the window never needs more than a few dozen SDs
  double maxW = std::min(200.0 * sd + 2000.0, 2e6);
  double top = NEG_INF;
  while (true) {
    if (R - L > maxW) return;
    t.L0 = L;
    t.lw.resize((std::size_t) (R - L + 1));
    top = NEG_INF;
    t.lw[0] = t.logw(L);
    for (std::size_t i = 0; i < t.lw.size(); ++i) {
      if (i > 0) {
        t.lw[i] = t.lw[i - 1] + tab_dlogw(fam, L + (double) i - 1.0,
                                          loglam, par);
      }
      if (t.lw[i] > top) top = t.lw[i];
    }
    // Log weights this large cannot be normalized to useful precision
    if (!R_FINITE(top) || std::fabs(top) > 1e12) return;
    std::size_t n = t.lw.size();
    double dr1 = t.lw[n - 1] - t.lw[n - 2], dr2 = t.lw[n - 2] - t.lw[n - 3];
    bool okR = t.lw[n - 1] < top - 50.0 && dr1 < 0 && dr1 <= dr2;
    double dl1 = t.lw[0] - t.lw[1], dl2 = t.lw[1] - t.lw[2];
    bool okL = (L == 0) || (t.lw[0] < top - 50.0 && dl1 < 0 && dl1 <= dl2);
    if (okR && okL) {
      double r = std::exp(dr1);
      t.ltail_r = t.lw[n - 1] + std::log(r) - std::log1p(-r);
      if (L > 0) {
        double rl = std::exp(dl1);
        t.ltail_l = t.lw[0] + std::log(rl) - std::log1p(-rl);
      } else {
        t.ltail_l = NEG_INF;
      }
      break;
    }
    double width = R - L;
    if (!okR) R += width;
    if (!okL) L = std::max(0.0, L - width);
  }
  double s0 = 0.0, s1 = 0.0, s2 = 0.0;
  for (std::size_t i = 0; i < t.lw.size(); ++i) {
    double w = std::exp(t.lw[i] - top), k = t.L0 + (double) i;
    s0 += w;
    s1 += w * k;
    s2 += w * k * k;
  }
  t.logZ = logaddexp(logaddexp(top + std::log(s0), t.ltail_r), t.ltail_l);
  t.mean = s1 / s0;
  t.var = std::max(s2 / s0 - t.mean * t.mean, 0.0);
  t.ok = true;
}

// log P(Y <= y) and log P(Y > y) from a filled table
double tab_lower_direct(const Table &t, double y) {
  double acc = NEG_INF;
  for (double k = 0; k <= y; k += 1.0) acc = logaddexp(acc, t.logw(k));
  return std::min(acc - t.logZ, 0.0);
}

double tab_upper(const Table &t, double y);

double tab_lower(const Table &t, double y) {
  double R = t.L0 + (double) t.lw.size() - 1.0;
  if (y < t.L0) return tab_lower_direct(t, y);
  if (y > R) return log1mexp(tab_upper(t, y));
  double acc = t.ltail_l;
  for (std::size_t i = 0; i <= (std::size_t) (y - t.L0); ++i) {
    acc = logaddexp(acc, t.lw[i]);
  }
  return std::min(acc - t.logZ, 0.0);
}

double tab_upper(const Table &t, double y) {
  double R = t.L0 + (double) t.lw.size() - 1.0;
  if (y < t.L0) return log1mexp(tab_lower_direct(t, y));
  if (y < R) {
    double acc = t.ltail_r;
    for (std::size_t i = (std::size_t) (y + 1 - t.L0); i < t.lw.size(); ++i) {
      acc = logaddexp(acc, t.lw[i]);
    }
    return std::min(acc - t.logZ, 0.0);
  }
  // Past the window the weights decrease ever faster: sum until negligible
  double acc = NEG_INF, prev = R_PosInf;
  for (double k = y + 1.0; k <= y + 1e6; k += 1.0) {
    double term = t.logw(k);
    acc = logaddexp(acc, term);
    if (term < acc - 40.0 && term < prev) return std::min(acc - t.logZ, 0.0);
    prev = term;
  }
  return NA_REAL;
}

// log(lambda) with mean mu; the table is left filled at the solution
double tab_solve(Table &t, int fam, double mu, double par) {
  double v = std::log(mu);
  if (fam == 1) v = std::log(std::max(mu * par, mu + par - 1.0));
  double lo = NEG_INF, hi = R_PosInf;
  for (int it = 0; it < 100; ++it) {
    tab_fill(t, fam, v, par);
    if (!t.ok) return NA_REAL;
    double f = t.mean - mu;
    if (std::fabs(f) <= 1e-13 * mu) return v;
    if (f < 0) lo = v; else hi = v;
    if (R_FINITE(lo) && R_FINITE(hi) &&
        hi - lo <= 1e-14 * std::max(1.0, std::fabs(v))) return v;
    double slope = ((fam == 0) ? par : 1.0) * t.var;
    double vn = (slope > 0) ? v - f / slope : NA_REAL;
    double step_max = 5.0;
    if (!R_FINITE(vn) || vn <= lo || vn >= hi ||
        std::fabs(vn - v) > step_max) {
      if (R_FINITE(lo) && R_FINITE(hi)) vn = 0.5 * (lo + hi);
      else vn = v + ((f < 0) ? step_max : -step_max);
    }
    if (std::fabs(vn - v) <= 1e-15 * std::max(1.0, std::fabs(v))) return vn;
    v = vn;
  }
  return NA_REAL;
}

// ---- discrete Weibull --------------------------------------------------------

struct DWSums { double E, dE, E2; };

// Sums over k >= 1 of f(k), k^b f(k) and (2k - 1) f(k), f(k) = exp(-c k^b).
// Terms are formed on the log scale (c k^b = exp(log c + b log k)), so
// extreme shapes give 0 or 1 rather than 0 * Inf.
DWSums dw_sums(double c, double b, bool second) {
  double lc = std::log(c);
  DWSums s = {0.0, 0.0, 0.0};
  if (!R_FINITE(lc) || !R_FINITE(b) || b <= 0) {
    s.E = s.dE = s.E2 = NA_REAL;
    return s;
  }
  // terms are negligible once c k^b >= 40 + c
  double lkneed = (std::log(40.0 + c) - lc) / b;
  double K0 = std::min(std::max(200.0, std::ceil(400.0 * b)), 1e5);
  double kneed = (lkneed > std::log(1e15)) ? 1e15 : std::ceil(std::exp(lkneed));
  double K = std::min(kneed, K0);
  for (double k = 1; k <= K; k += 1.0) {
    double ckb = std::exp(lc + b * std::log(k));
    double f = std::exp(-ckb);
    if (f == 0.0) break;
    s.E += f;
    s.dE += ckb / c * f;
    if (second) s.E2 += (2.0 * k - 1.0) * f;
  }
  if (kneed > K0) {
    // Euler-Maclaurin remainder beyond K:
    //   sum_{k > K} g(k) = int_K^Inf g - g(K) / 2 - g'(K) / 12
    //                      + g'''(K) / 720 - ...
    double u = c * std::pow(K, b);
    auto I = [&](double j) {   // int_K^Inf x^(j - 1) exp(-c x^b) dx
      return std::exp(std::lgamma(j / b) - (j / b) * std::log(c) -
                      std::log(b) + R::pgamma(u, j / b, 1.0, 0, 1));
    };
    double f = std::exp(-u);
    double g1 = c * b * std::pow(K, b - 1.0);
    double g2 = c * b * (b - 1.0) * std::pow(K, b - 2.0);
    double g3 = c * b * (b - 1.0) * (b - 2.0) * std::pow(K, b - 3.0);
    double f1 = -g1 * f;
    double f2 = (g1 * g1 - g2) * f;
    double f3 = (-g1 * g1 * g1 + 3.0 * g1 * g2 - g3) * f;
    s.E += I(1.0) - f / 2.0 - f1 / 12.0 + f3 / 720.0;
    double Kb = std::pow(K, b);
    double h = Kb * f, h1 = b * std::pow(K, b - 1.0) * f + Kb * f1;
    s.dE += I(b + 1.0) - h / 2.0 - h1 / 12.0;
    if (second) {
      double q = 2.0 * K - 1.0;
      double gg = q * f, gg1 = 2.0 * f + q * f1;
      double gg3 = 6.0 * f2 + q * f3;
      s.E2 += 2.0 * I(2.0) - I(1.0) - gg / 2.0 - gg1 / 12.0 + gg3 / 720.0;
    }
  }
  return s;
}

// log(c) with mean mu
double dw_solve(double mu, double b) {
  double u = b * (std::lgamma(1.0 + 1.0 / b) - std::log(mu + 0.5));
  double lo = NEG_INF, hi = R_PosInf;
  for (int it = 0; it < 100; ++it) {
    double c = std::exp(u);
    DWSums s = dw_sums(c, b, false);
    if (ISNAN(s.E)) return NA_REAL;
    double f = s.E - mu;
    if (std::fabs(f) <= 1e-13 * mu) return u;
    if (f > 0) lo = u; else hi = u;    // E decreases in c
    if (R_FINITE(lo) && R_FINITE(hi) &&
        hi - lo <= 1e-14 * std::max(1.0, std::fabs(u))) return u;
    double slope = -c * s.dE;
    double un = (slope < 0) ? u - f / slope : NA_REAL;
    if (!R_FINITE(un) || un <= lo || un >= hi || std::fabs(un - u) > 5.0) {
      if (R_FINITE(lo) && R_FINITE(hi)) un = 0.5 * (lo + hi);
      else un = u + ((f > 0) ? 5.0 : -5.0);
    }
    if (std::fabs(un - u) <= 1e-15 * std::max(1.0, std::fabs(u))) return un;
    u = un;
  }
  return NA_REAL;
}

// ---- gamma renewal count (unit window) ---------------------------------------

struct GCSums { double E, dE, E2; };

// Work is capped (5000 terms past the leading ones); beyond it the result is
// NA (only reached for extreme shapes, roughly a < 0.005, i.e. a dispersion
// index above about 200).
GCSums gc_sums(double s, double a, bool second) {
  GCSums r = {0.0, 0.0, 0.0};
  if (!R_FINITE(s) || s <= 0 || !R_FINITE(a) || a <= 0) {
    r.E = r.dE = r.E2 = NA_REAL;
    return r;
  }
  // P(Gamma(k a) <= s) is 1 to double precision while k a is more than
  // 10 sd + 40 below s: those terms are counted without evaluation
  double klow = std::floor((s - 10.0 * std::sqrt(s) - 40.0) / a);
  if (klow < 0) klow = 0;
  r.E = klow;
  if (second) r.E2 = klow * klow;   // sum_{k <= klow} (2k - 1)
  double prev = 2.0;
  bool done = false;
  for (double k = klow + 1.0; k <= klow + 5000.0; k += 1.0) {
    double p = R::pgamma(s, k * a, 1.0, 1, 0);
    r.E += p;
    r.dE += R::dgamma(s, k * a, 1.0, 0);
    if (second) r.E2 += (2.0 * k - 1.0) * p;
    // Past the bulk, the terms decrease ever faster, so the remainder is
    // below the geometric bound p r / (1 - r) with r = p / prev
    if (k * a > s && p < prev) {
      double ratio = p / prev;
      if (p == 0.0 || p / (1.0 - ratio) < 1e-16 * std::max(r.E, 1e-300)) {
        done = true;
        break;
      }
    }
    prev = p;
  }
  if (!done) r.E = r.dE = r.E2 = NA_REAL;
  return r;
}

// log(s) with mean mu
double gc_solve(double mu, double a) {
  double s0 = (mu < 1.0) ? std::pow(mu * std::exp(std::lgamma(a + 1.0)),
                                   1.0 / a) : a * mu + (a - 1.0) / 2.0;
  if (!(s0 > 0) || !R_FINITE(s0)) s0 = a * mu;
  double v = std::log(s0);
  double lo = NEG_INF, hi = R_PosInf;
  for (int it = 0; it < 100; ++it) {
    double s = std::exp(v);
    GCSums r = gc_sums(s, a, false);
    if (ISNAN(r.E)) return NA_REAL;
    double f = r.E - mu;
    if (std::fabs(f) <= 1e-13 * mu) return v;
    if (f < 0) lo = v; else hi = v;     // E increases in s
    if (R_FINITE(lo) && R_FINITE(hi) &&
        hi - lo <= 1e-14 * std::max(1.0, std::fabs(v))) return v;
    double slope = s * r.dE;
    double vn = (slope > 0) ? v - f / slope : NA_REAL;
    if (!R_FINITE(vn) || vn <= lo || vn >= hi || std::fabs(vn - v) > 5.0) {
      if (R_FINITE(lo) && R_FINITE(hi)) vn = 0.5 * (lo + hi);
      else vn = v + ((f < 0) ? 5.0 : -5.0);
    }
    if (std::fabs(vn - v) <= 1e-15 * std::max(1.0, std::fabs(v))) return vn;
    v = vn;
  }
  return NA_REAL;
}

}  // namespace

// DPO (fam = 0) or hyper-Poisson (fam = 1) with mean mu and shape par.
// what: 0 log P(Y = y), 1 log P(Y <= y), 2 log P(Y > y), 3 variance,
// 4 log(lambda). Inputs are recycled and validated in R (mu > 0, par > 0,
// y whole and >= 0); a failed solve gives NA.
// [[Rcpp::export]]
NumericVector tab_mp_cpp(NumericVector y, NumericVector mu, NumericVector par,
                         int fam, int what) {
  R_xlen_t n = mu.size();
  NumericVector out(n);
  Table t;
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i % 256 == 0) Rcpp::checkUserInterrupt();
    double v = tab_solve(t, fam, mu[i], par[i]);
    if (ISNAN(v)) {
      out[i] = NA_REAL;
      continue;
    }
    if (what == 0) out[i] = t.logw(y[i]) - t.logZ;
    else if (what == 1) out[i] = tab_lower(t, y[i]);
    else if (what == 2) out[i] = tab_upper(t, y[i]);
    else if (what == 3) out[i] = t.var;
    else out[i] = v;
  }
  return out;
}

// Discrete Weibull with mean mu and shape b. what: 0 log(c), c = -log(theta);
// 1 variance.
// [[Rcpp::export]]
NumericVector dw_mp_cpp(NumericVector mu, NumericVector b, int what) {
  R_xlen_t n = mu.size();
  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i % 256 == 0) Rcpp::checkUserInterrupt();
    double u = dw_solve(mu[i], b[i]);
    if (what == 0 || ISNAN(u)) {
      out[i] = u;
      continue;
    }
    DWSums s = dw_sums(std::exp(u), b[i], true);
    out[i] = std::max(s.E2 - s.E * s.E, 0.0);
  }
  return out;
}

// Gamma count (unit window) with mean mu and shape a. what: 0 log(rate);
// 1 variance.
// [[Rcpp::export]]
NumericVector gc_mp_cpp(NumericVector mu, NumericVector a, int what) {
  R_xlen_t n = mu.size();
  NumericVector out(n);
  for (R_xlen_t i = 0; i < n; ++i) {
    if (i % 256 == 0) Rcpp::checkUserInterrupt();
    double v = gc_solve(mu[i], a[i]);
    if (what == 0 || ISNAN(v)) {
      out[i] = v;
      continue;
    }
    GCSums r = gc_sums(std::exp(v), a[i], true);
    out[i] = std::max(r.E2 - r.E * r.E, 0.0);
  }
  return out;
}
