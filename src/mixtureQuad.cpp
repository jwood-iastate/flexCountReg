// Log probabilities of mixed count distributions by adaptive quadrature.
//
//   P(Y = y) = integral of p(y | lambda(t)) g(t) dt,  lambda(t) = c exp(s t)
//
// where t is the log of the mixing variable (or a standard normal for the
// lognormal mixtures), p is the Poisson or Poisson-Lindley kernel and g the
// density of t. Kinds:
//   0  PLN  Poisson-lognormal:        t ~ N(0, 1), c = mean, s = sigma
//   1  PGE  Poisson-generalized exp.: x = exp(t) ~ GE(shape a, scale 1),
//           c = mean / A, A = digamma(a + 1) - digamma(1), s = 1
//   2  PW   Poisson-Weibull:          x = exp(t) ~ Weibull(alpha, 1),
//           c = mean / Gamma(1 + 1 / alpha), s = 1
//   3  PLL  Poisson-Lindley-lognormal: kernel Poisson-Lindley with rate
//           lambda, t ~ N(0, 1), c = mean theta (theta + 1) / (theta + 2),
//           s = sigma
//
// The log integrand h(t) is maximized by a safeguarded Newton method with
// analytic derivatives; the integral is then computed on the log scale with
// a sinh-sinh (double exponential) rule centred at the mode and scaled by
// the curvature there. The rule is deterministic and resolves the region
// where the integrand has its mass for every count, including counts far in
// the tail where fixed simulation draws carry no information.
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <cmath>
#include <vector>
using namespace Rcpp;

namespace {

const double LOG_2PI = 1.8378770664093454836;

struct MixSpec {
  int kind;
  double y, logc, s, a, theta, lgy1;
};

// a * b with 0 * (+-Inf) = 0 (a term whose coefficient is zero is absent)
inline double mul0(double a, double b) { return a == 0.0 ? 0.0 : a * b; }

// h(t) and, when d1/d2 are given, its first two derivatives
double mix_h(const MixSpec &m, double t, double *d1 = nullptr,
             double *d2 = nullptr) {
  double L = m.logc + m.s * t;
  double lam = std::exp(L);
  double k, k1 = 0, k2 = 0;
  if (m.kind == 3) {
    // Poisson-Lindley kernel with rate lam:
    // theta^2 / (theta + 1) lam^y (theta + lam + y + 1) / (theta + lam)^(y+2)
    double th = m.theta;
    double B = th + lam + m.y + 1.0, C = th + lam;
    k = 2.0 * std::log(th) - std::log1p(th) + mul0(m.y, L) + std::log(B) -
      (m.y + 2.0) * std::log(C);
    if (d1) {
      k1 = m.s * (m.y + lam / B - (m.y + 2.0) * lam / C);
      k2 = m.s * m.s * lam * ((th + m.y + 1.0) / (B * B) -
        (m.y + 2.0) * th / (C * C));
    }
  } else {
    k = mul0(m.y, L) - lam - m.lgy1;
    if (d1) {
      k1 = m.s * (m.y - lam);
      k2 = -m.s * m.s * lam;
    }
  }
  double g, g1 = 0, g2 = 0;
  if (m.kind == 0 || m.kind == 3) {
    g = -0.5 * t * t - 0.5 * LOG_2PI;
    g1 = -t;
    g2 = -1.0;
  } else if (m.kind == 1) {
    // GE(a, 1) in t = log(x): log a + t - x + (a - 1) log(1 - exp(-x))
    double x = std::exp(t);
    double om = -std::expm1(-x);          // 1 - exp(-x)
    // log(1 - exp(-x)), = t - x / 2 + O(x^2) for tiny x (om underflows)
    double lom = (x < 1e-8) ? t - x / 2.0 : std::log(om);
    g = std::log(m.a) + t - x + mul0(m.a - 1.0, lom);
    if (d1) {
      double q = (x > 700) ? 0.0 : x / std::expm1(x);   // x / (e^x - 1)
      double r = (om > 0) ? x / om : 1.0;               // x e^x / (e^x - 1)
      g1 = 1.0 - x + (m.a - 1.0) * q;
      g2 = -x + (m.a - 1.0) * q * (1.0 - r);
    }
  } else {
    // Weibull(alpha, 1) in t = log(x): log alpha + alpha t - x^alpha
    double xa = std::exp(m.a * t);
    g = std::log(m.a) + m.a * t - xa;
    g1 = m.a - m.a * xa;
    g2 = -m.a * m.a * xa;
  }
  if (d1) {
    *d1 = k1 + g1;
    *d2 = k2 + g2;
  }
  double h = k + g;
  return ISNAN(h) ? R_NegInf : h;
}

// Mode of h by safeguarded Newton iterations with backtracking
double mix_mode(const MixSpec &m) {
  // Start where the rate matches the count (clamped to a sensible range)
  double t = (std::log(std::max(m.y, 0.5)) - m.logc) / m.s;
  if (m.kind == 0 || m.kind == 3) t = std::max(-8.0, std::min(8.0, t));
  else t = std::max(-30.0, std::min(10.0, t));
  double d1, d2;
  double ht = mix_h(m, t, &d1, &d2);
  for (int it = 0; it < 300; ++it) {
    double step;
    if (d2 < 0 && R_FINITE(d2)) step = -d1 / d2;
    else step = (d1 > 0) ? 1.0 : -1.0;
    if (!R_FINITE(step)) step = (d1 > 0) ? 1.0 : -1.0;
    step = std::max(-5.0, std::min(5.0, step));
    double tn = t + step, hn = mix_h(m, tn);
    int bt = 0;
    while (!(hn >= ht - 1e-12 * std::fabs(ht)) && bt < 60) {
      step /= 2.0;
      tn = t + step;
      hn = mix_h(m, tn);
      ++bt;
    }
    if (!(hn >= ht - 1e-12 * std::fabs(ht))) break;
    t = tn;
    ht = mix_h(m, t, &d1, &d2);
    if (std::fabs(step) < 1e-11 * std::max(1.0, std::fabs(t))) break;
  }
  return t;
}

// Nodes and log weights of the sinh-sinh rule on the real line
struct Rule {
  std::vector<double> u, logw;
  Rule() {
    const double hstep = 1.0 / 32.0, vmax = 3.25;
    const double hp = M_PI / 2.0;
    int K = (int) std::ceil(vmax / hstep);
    for (int k = -K; k <= K; ++k) {
      double v = k * hstep;
      double sv = std::sinh(v);
      u.push_back(std::sinh(hp * sv));
      logw.push_back(std::log(hstep * hp * std::cosh(v)) +
                     std::log(std::cosh(hp * sv)));
    }
  }
};

const Rule &rule() {
  static const Rule r;
  return r;
}

double mix_logpmf_one(MixSpec m) {
  double tm = mix_mode(m);
  double d1, d2;
  double hm = mix_h(m, tm, &d1, &d2);
  if (!R_FINITE(hm)) return hm;
  double sc = (d2 < 0 && R_FINITE(d2)) ? 1.0 / std::sqrt(-d2) : 1.0;
  const Rule &r = rule();
  // Sum of exp(h(t_k) - h(t*) + log w_k), from the centre node outwards in
  // each direction. h is unimodal, so once a term is below 1e-20 of the
  // peak and h is falling, the remaining terms in that direction are
  // smaller still and are skipped.
  int K = (int) r.u.size(), mid = K / 2;
  double total = std::exp(r.logw[mid]);   // u = 0: h(t*) - h(t*) = 0
  for (int dir = -1; dir <= 1; dir += 2) {
    double prev = 0.0;
    for (int k = mid + dir; k >= 0 && k < K; k += dir) {
      double hk = mix_h(m, tm + sc * r.u[k]) - hm;
      double term = hk + r.logw[k];
      if (term > -745) total += std::exp(term);
      if (term < -46.0 && hk < prev) break;
      prev = hk;
    }
  }
  return hm + std::log(sc) + std::log(total);
}

}  // namespace

// Log P(Y = x) for vectors of equal length (recycled in R). par1 is sigma
// (PLN), shape (PGE), alpha (PW) or theta (PLL); par2 is sigma (PLL).
// Invalid parameters give NaN; x that is not a nonnegative whole number
// gives -Inf; NA gives NA.
// [[Rcpp::export]]
NumericVector mix_logpmf_cpp(NumericVector x, NumericVector mean,
                             NumericVector par1, NumericVector par2,
                             int kind) {
  int n = x.size();
  NumericVector out(n);
  for (int i = 0; i < n; ++i) {
    double xi = x[i], mu = mean[i], p1 = par1[i], p2 = par2[i];
    if (ISNAN(xi) || ISNAN(mu) || ISNAN(p1) || (kind == 3 && ISNAN(p2))) {
      out[i] = NA_REAL;
      continue;
    }
    bool ok = R_FINITE(mu) && mu > 0 && R_FINITE(p1);
    if (kind == 0) ok = ok && p1 >= 0;
    else ok = ok && p1 > 0;
    if (kind == 3) ok = ok && R_FINITE(p2) && p2 >= 0;
    if (!ok) {
      out[i] = R_NaN;
      continue;
    }
    if (!R_FINITE(xi) || xi < 0 || xi != std::floor(xi)) {
      out[i] = R_NegInf;
      continue;
    }
    MixSpec m;
    m.kind = kind;
    m.y = xi;
    m.lgy1 = std::lgamma(xi + 1.0);
    m.a = p1;
    m.theta = p1;
    switch (kind) {
    case 0:
      m.logc = std::log(mu); m.s = p1; break;
    case 1:
      m.logc = std::log(mu) -
        std::log(R::digamma(p1 + 1.0) - R::digamma(1.0));
      m.s = 1.0; break;
    case 2:
      m.logc = std::log(mu) - std::lgamma(1.0 + 1.0 / p1);
      m.s = 1.0; break;
    default:
      m.logc = std::log(mu) + std::log(p1) + std::log1p(p1) -
        std::log(p1 + 2.0);
      m.s = p2;
    }
    if ((kind == 0 || kind == 3) && m.s == 0.0) {
      // No mixing: the kernel at t = 0 (normal density omitted)
      out[i] = mix_h(m, 0.0) + 0.5 * LOG_2PI;
      continue;
    }
    out[i] = mix_logpmf_one(m);
    if ((i & 0x3FF) == 0) Rcpp::checkUserInterrupt();
  }
  return out;
}
