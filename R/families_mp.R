# Helpers for the regression families added in 0.1.2 (registry entries in
# families.R). Every family is parameterized by its mean mu = exp(X beta)
# (times the reporting probability), so coefficients are effects on the
# expected count. Families whose mean has no closed form (DW, DPO,
# HYPERPOIS, GAMMACOUNT) find their natural parameter from mu numerically
# (src/meanParam.cpp).
#
# The registry passes the distribution parameters as exp(W delta). Shares
# rho in [0, 1] (DELAPORTE, HERMITE*) use a logit link: rho =
# plogis(log(alpha)), labelled "logit(rho)".

.fcr_logit_par <- function(alpha) stats::plogis(log(alpha))

# Numerical domain of the DW, DPO, HYPERPOIS and GAMMACOUNT regression
# families: shape/precision parameter in [1e-6, 1e6] and mean at most 1e9.
# Outside it the distributions are numerically degenerate (e.g. a double
# Poisson with phi = 1e8 is a point mass) or the mean equation cannot be
# solved to useful precision; the density is then treated as 0, so
# optimizers step back instead of spending their time there.
# For DPO the normalizing sums hold at most 2e6 terms (src/meanParam.cpp),
# which also limits mu / phi to about 1e10.
.fcr_mp_par_range <- c(1e-6, 1e6)
.fcr_mp_mu_max <- 1e9

# Recycle (y, mu, par) and flag the rows with valid mean-model parameters
.fcr_mp_args <- function(y, mu, par) {
  r <- .fcr_recycle(y, mu, par)
  ok <- !is.na(r[[1]]) & is.finite(r[[2]]) & r[[2]] > 0 &
    r[[2]] <= .fcr_mp_mu_max &
    is.finite(r[[3]]) & r[[3]] >= .fcr_mp_par_range[1] &
    r[[3]] <= .fcr_mp_par_range[2]
  list(y = as.numeric(r[[1]]), mu = as.numeric(r[[2]]),
       par = as.numeric(r[[3]]), ok = ok)
}

# Apply fun(y, mu, par) to valid rows; NaN elsewhere
.fcr_mp_apply <- function(y, mu, par, fun) {
  a <- .fcr_mp_args(y, mu, par)
  out <- rep(NaN, length(a$y))
  if (any(a$ok)) out[a$ok] <- fun(a$y[a$ok], a$mu[a$ok], a$par[a$ok])
  out
}

# ---- discrete Weibull (mean mu, shape b) ------------------------------------

.dw_mp_c <- function(mu, b) exp(dw_mp_cpp(mu, b, 0L))

.dw_mp_logpmf <- function(y, mu, b) {
  .fcr_mp_apply(y, mu, b, function(y, mu, b)
    .dweib_logpmf_c(y, .dw_mp_c(mu, b), b))
}

.dw_mp_logcdf <- function(q, mu, b, lower) {
  .fcr_mp_apply(floor(q), mu, b, function(y, mu, b)
    .dweib_logcdf_c(y, .dw_mp_c(mu, b), b, lower))
}

.dw_mp_var <- function(mu, b) {
  .fcr_mp_apply(0, mu, b, function(y, mu, b) dw_mp_cpp(mu, b, 1L))
}

.dw_mp_rand <- function(n, mu, b) {
  r <- .fcr_recycle(seq_len(n), mu, b)
  cc <- .dw_mp_c(r[[2]], r[[3]])
  floor((stats::rexp(n) / cc)^(1 / r[[3]]))
}

# ---- double Poisson (fam 0) and hyper-Poisson (fam 1) -------------------------

.tab_mp <- function(y, mu, par, fam, what) {
  .fcr_mp_apply(y, mu, par, function(y, mu, par)
    tab_mp_cpp(y, mu, par, fam, what))
}

.dpo_mp_rand <- function(n, mu, phi) {
  r <- .fcr_recycle(seq_len(n), mu, phi)
  lambda <- exp(.tab_mp(0, r[[2]], r[[3]], 0L, 4L))
  rdpo(n, lambda, r[[3]])
}

.hp_mp_rand <- function(n, mu, gamma) {
  r <- .fcr_recycle(seq_len(n), mu, gamma)
  lambda <- exp(.tab_mp(0, r[[2]], r[[3]], 1L, 4L))
  rhpois(n, lambda, r[[3]])
}

# ---- gamma count (mean mu over a unit window, shape a) ----------------------

# Parameters of the gamma count with mean mu: the rate is carried on the log
# scale as well (`logs`), since for small shapes it underflows (a mean of
# 0.5 at shape 0.001 needs a rate of about exp(-1100))
.gc_mp_pars <- function(mu, a) {
  v <- gc_mp_cpp(mu, a, 0L)
  list(shape = a, rate = exp(v), time = 1, logs = v)
}

.gc_mp_logpmf <- function(y, mu, a) {
  .fcr_mp_apply(y, mu, a, function(y, mu, a)
    .gcount_logpmf(y, .gc_mp_pars(mu, a)))
}

.gc_mp_logcdf <- function(q, mu, a, lower) {
  .fcr_mp_apply(floor(q), mu, a, function(y, mu, a)
    .gcount_logcdf(y, .gc_mp_pars(mu, a), lower))
}

.gc_mp_var <- function(mu, a) {
  .fcr_mp_apply(0, mu, a, function(y, mu, a) gc_mp_cpp(mu, a, 1L))
}

# Draws by inversion of the CDF (as rgcount()), with the log rate so that
# rates that underflow work
.gc_mp_rand <- function(n, mu, a) {
  r <- .fcr_recycle(seq_len(n), mu, a)
  p <- .gc_mp_pars(r[[2]], r[[3]])
  .fcr_q(stats::runif(n), p, function(p) is.finite(p$logs) &
           is.finite(p$shape) & p$shape > 0,
         .gcount_logcdf, .gcount_qstart, TRUE, FALSE)
}

# ---- Polya-Aeppli as Poisson-Tweedie with power 1.5 ---------------------------

.pa_mp <- function(mu, phi) .pt15_to_paeppli(mu, phi)
