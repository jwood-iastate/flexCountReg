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
//   4  NBGE negative binomial-generalized exponential: kernel NB2 with
//           mean lambda and size r = 1 / alpha, x = exp(t) ~ GE(a, 1),
//           c = mean / A (as PGE), s = 1
//   5  PW by scale: as PW with c = lambda sigma given directly ("mean"
//           is then lambda sigma), so that Weibull shapes whose
//           Gamma(1 + 1 / alpha) overflows can be used
//
// The log integrand h(t) is maximized by a safeguarded Newton method with
// analytic derivatives; the integral is then computed on the log scale with
// a sinh-sinh (double exponential) rule centered at the mode and scaled by
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
  double r, lgr;   // NBGE: NB size 1 / alpha and lgamma(y + r) - lgamma(r)
};

// lgamma(y + r) - lgamma(r) for whole y >= 0: summed exactly for large r
// (where the difference of two large lgamma values loses precision)
double lgamma_ratio(double y, double r) {
  if (r < 1e6 || y > 1000) return std::lgamma(y + r) - std::lgamma(r);
  double out = 0.0;
  for (double j = 0; j < y; j += 1.0) out += std::log(r + j);
  return out;
}

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
  } else if (m.kind == 4) {
    // NB2 kernel with mean lam and size r:
    // lgamma(y + r) - lgamma(r) - lgamma(y + 1) + r log(r / (r + lam))
    //   + y log(lam / (r + lam))
    double rl = m.r + lam;
    k = m.lgr - m.lgy1 - m.r * std::log1p(lam / m.r) +
      mul0(m.y, L - std::log(rl));
    if (d1) {
      k1 = m.s * (m.y - (m.y + m.r) * lam / rl);
      k2 = -m.s * m.s * (m.y + m.r) * m.r * lam / (rl * rl);
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
  } else if (m.kind == 1 || m.kind == 4) {
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

// Mode of h for the Poisson-Weibull mixture. With A = y + alpha,
//   h'(t) = A - c exp(t) - alpha exp(alpha t)
// is strictly decreasing, so the mode is its unique root. It is bracketed
// in closed form: where either term alone reaches A, h' <= 0; where both
// are at most A / 2, h' >= 0. Newton steps are taken inside the bracket,
// with bisection whenever a step leaves it. A start from the count and
// fixed step limits (as for the other kinds) can stall: for large alpha
// the steps are about 1 / alpha long, and for small alpha the mode can be
// thousands of units (about 1 / alpha) away.
double pw_mode(const MixSpec &m) {
  double lA = std::log(m.y + m.a), la = std::log(m.a);
  double hi = std::min(lA - m.logc, (lA - la) / m.a);
  double lo = std::min(lA - M_LN2 - m.logc, (lA - M_LN2 - la) / m.a);
  double t = 0.5 * (lo + hi);
  for (int it = 0; it < 400; ++it) {
    double d1, d2;
    mix_h(m, t, &d1, &d2);
    if (d1 == 0.0) return t;
    if (d1 > 0) lo = t; else hi = t;
    if (hi - lo <= 4e-16 * std::max(1.0, std::fabs(t))) break;
    double tn = (d2 < 0 && R_FINITE(d2)) ? t - d1 / d2 : R_NaN;
    if (!(tn > lo && tn < hi)) tn = 0.5 * (lo + hi);
    if (std::fabs(tn - t) <= 1e-15 * std::max(1.0, std::fabs(t))) {
      t = tn;
      break;
    }
    t = tn;
  }
  return t;
}

// Mode of h by safeguarded Newton iterations with backtracking
double mix_mode(const MixSpec &m) {
  if (m.kind == 2) return pw_mode(m);
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

// A log probability is at most 0. Rounding can give a tiny positive value;
// anything larger (or not a number) means the rule failed, which is
// reported as NaN rather than as a probability above 1.
double mix_check_logp(double out) {
  if (out > 0) out = (out < 1e-10) ? 0.0 : R_NaN;
  if (ISNAN(out)) out = R_NaN;
  return out;
}

// Sum over the sinh-sinh nodes v = +-(v0 + j hs), j = 0, 1, ..., |v| <= 3.25,
// v0 > 0, of exp(h(t(v)) - h(t*)) times dt / dv / sc, with the same early
// stopping as the fixed rule (used to halve the step of that rule)
double de_sum_offset(const MixSpec &m, double tm, double sc, double hm,
                     double v0, double hs) {
  const double hp = M_PI / 2.0, vmax = 3.25;
  double total = 0.0;
  for (int dir = -1; dir <= 1; dir += 2) {
    double prev = 0.0;
    for (int j = 0; ; ++j) {
      double v = v0 + j * hs;
      if (v > vmax) break;
      double sv = std::sinh(v);
      double hk = mix_h(m, tm + sc * dir * std::sinh(hp * sv)) - hm;
      double term = hk + std::log(hp * std::cosh(v)) +
        std::log(std::cosh(hp * sv));
      if (term > -745) total += std::exp(term);
      if (term < -46.0 && hk < prev) break;
      prev = hk;
    }
  }
  return total;
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
  double even = total;                    // nodes of the rule with step 2 h
  for (int dir = -1; dir <= 1; dir += 2) {
    double prev = 0.0;
    for (int k = mid + dir; k >= 0 && k < K; k += dir) {
      double hk = mix_h(m, tm + sc * r.u[k]) - hm;
      double term = hk + r.logw[k];
      if (term > -745) {
        double e = std::exp(term);
        total += e;
        if ((k - mid) % 2 == 0) even += e;
      }
      if (term < -46.0 && hk < prev) break;
      prev = hk;
    }
  }
  if (m.kind == 2) {
    // Poisson-Weibull: error check. The error of a double exponential rule
    // with step h is about the square of its error with step 2 h, so when
    // the sums with steps h and 2 h (the even nodes) agree to 1e-7, the
    // finer one is accurate to about 1e-14. Otherwise the step is halved
    // (reusing the earlier nodes) until two successive sums agree that
    // well. For small alpha and y = 0 the integrand has a steep edge (where
    // the Poisson rate passes 1) inside the much wider bulk of the Weibull
    // mixing density, and the fixed step left a relative error of up to
    // about 1e-4 in 1 - P(Y = 0).
    double hs = 1.0 / 32.0;   // the step of rule()
    if (!(std::fabs(total - 2.0 * even) <= 1e-7 * total)) {
      for (int lev = 0; lev < 6; ++lev) {
        double Sn = 0.5 * total +
          0.5 * hs * de_sum_offset(m, tm, sc, hm, 0.5 * hs, hs);
        hs *= 0.5;
        bool conv = std::fabs(Sn - total) <= 1e-7 * Sn;
        total = Sn;
        if (conv) break;
      }
    }
  }
  return mix_check_logp(hm + std::log(sc) + std::log(total));
}

}  // namespace

// Log P(Y = x) for vectors of equal length (recycled in R). par1 is sigma
// (PLN), shape (PGE, NBGE), alpha (PW, kind 5) or theta (PLL); par2 is
// sigma (PLL) or the NB dispersion alpha > 0 (NBGE).
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
    if (ISNAN(xi) || ISNAN(mu) || ISNAN(p1) ||
        ((kind == 3 || kind == 4) && ISNAN(p2))) {
      out[i] = NA_REAL;
      continue;
    }
    bool ok = R_FINITE(mu) && mu > 0 && R_FINITE(p1);
    if (kind == 0) ok = ok && p1 >= 0;
    else ok = ok && p1 > 0;
    if (kind == 3) ok = ok && R_FINITE(p2) && p2 >= 0;
    if (kind == 4) ok = ok && R_FINITE(p2) && p2 > 0;
    if (!ok) {
      out[i] = R_NaN;
      continue;
    }
    if (!R_FINITE(xi) || xi < 0 || xi != std::floor(xi)) {
      out[i] = R_NegInf;
      continue;
    }
    MixSpec m;
    m.kind = (kind == 5) ? 2 : kind;
    m.r = (kind == 4) ? 1.0 / p2 : 0.0;
    m.lgr = (kind == 4) ? lgamma_ratio(xi, m.r) : 0.0;
    m.y = xi;
    m.lgy1 = std::lgamma(xi + 1.0);
    m.a = p1;
    m.theta = p1;
    switch (kind) {
    case 0:
      m.logc = std::log(mu); m.s = p1; break;
    case 1:
    case 4:
      m.logc = std::log(mu) -
        std::log(R::digamma(p1 + 1.0) - R::digamma(1.0));
      m.s = 1.0; break;
    case 2:
      m.logc = std::log(mu) - std::lgamma(1.0 + 1.0 / p1);
      m.s = 1.0; break;
    case 5:
      m.logc = std::log(mu);
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
