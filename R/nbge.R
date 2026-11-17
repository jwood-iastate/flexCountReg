# Negative binomial-generalized exponential mixture (NBGE)

.nbge_valid <- function(p) {
  is.finite(p$mu) & p$mu > 0 & is.finite(p$shape) & p$shape > 0 &
    is.finite(p$alpha) & p$alpha >= 0
}

# Log PMF by adaptive quadrature (src/mixtureQuad.cpp, kind 4); alpha = 0 is
# the Poisson-generalized exponential (kind 1).
.nbge_logpmf <- function(y, p) {
  out <- numeric(length(y))
  pois <- p$alpha == 0
  if (any(pois)) {
    out[pois] <- mix_logpmf_cpp(y[pois], p$mu[pois], p$shape[pois],
                                rep(0, sum(pois)), 1L)
  }
  if (any(!pois)) {
    out[!pois] <- mix_logpmf_cpp(y[!pois], p$mu[!pois], p$shape[!pois],
                                 p$alpha[!pois], 4L)
  }
  out
}

# Tails from the PMF. The lower tail P(Y <= y) is summed over 0..y for all
# elements at once (one vectorized quadrature call). The upper tail is 1
# minus it while that is at least 1e-3 (relative error below about 1e-13);
# smaller upper tails are summed directly from y + 1 until the terms are
# negligible, so they keep their relative accuracy.
.nbge_logcdf <- function(y, p, lower) {
  n <- length(y)
  out <- numeric(n)
  small <- y <= 5000
  low <- rep(NA_real_, n)
  if (any(small)) {
    idx <- which(small)
    reps <- y[idx] + 1
    rows <- rep(idx, reps)
    k <- sequence(reps) - 1
    lp <- .nbge_logpmf(k, .fcr_subset_pars(p, rows))
    m <- tapply(lp, rows, max)
    m_r <- m[as.character(rows)]
    s <- tapply(exp(lp - m_r), rows, sum)
    low[idx] <- pmin(as.vector(m + log(s)), 0)
  }
  use_low <- !is.na(low) & (lower | low <= log1p(-1e-3))
  if (any(use_low)) {
    out[use_low] <- if (lower) low[use_low] else .fcr_log1mexp(low[use_low])
  }
  for (i in which(!use_low)) {
    pi <- .fcr_subset_pars(p, i)
    ls <- .fcr_log_upper_tail(y[i], function(x)
      .nbge_logpmf(x, .fcr_subset_pars(pi, rep(1L, length(x)))))
    out[i] <- if (lower) .fcr_log1mexp(ls) else ls
  }
  out
}

.nbge_moments <- function(mu, shape, alpha) {
  A <- digamma(shape + 1) - digamma(1)
  B <- trigamma(1) - trigamma(shape + 1)
  M2 <- 1 + B / A^2
  list(mean = mu, var = mu + ((1 + alpha) * M2 - 1) * mu^2)
}

.nbge_qstart <- function(p, llow) {
  m <- .nbge_moments(p$mu, p$shape, p$alpha)
  .fcr_qstart_normal(m$mean, m$var, llow)
}

# G = X / A with X ~ GE(shape, 1) by inversion, then Y | G ~ NB2
.nbge_rand <- function(n, p) {
  A <- digamma(p$shape + 1) - digamma(1)
  x <- -log1p(-stats::runif(n)^(1 / p$shape))
  lam <- p$mu * x / A
  out <- stats::rpois(n, lam)
  nb <- p$alpha > 0
  out[nb] <- stats::rnbinom(sum(nb), size = 1 / p$alpha[nb], mu = lam[nb])
  out
}

#' Negative Binomial-Generalized Exponential Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the negative binomial-generalized exponential (NBGE) mixture: a
#' negative binomial (NB2) count whose mean is multiplied by a generalized
#' exponential random effect with mean one.
#'
#' @inheritParams DiscreteWeibull
#' @param mu vector of (positive) means.
#' @param shape vector of (positive) shape parameters \eqn{k} of the
#'   generalized exponential mixing distribution.
#' @param alpha vector of (non-negative) NB2 dispersion parameters;
#'   \code{alpha = 0} gives the Poisson-generalized exponential distribution
#'   (\code{\link{dpge}}).
#'
#' @details
#' The distribution is defined hierarchically:
#' \deqn{X \sim \mathrm{GE}(k, 1), \quad F_X(x) = (1 - e^{-x})^k, \quad
#'   G = X / A, \quad A = \psi(k + 1) - \psi(1),}
#' \deqn{Y \mid G \sim \mathrm{NB2}(\mathrm{mean} = \mu G,
#'   \mathrm{dispersion} = \alpha),}
#' where \eqn{\psi} is the digamma function, so that \eqn{E[G] = 1} and the
#' conditional variance is \eqn{\mu G + \alpha (\mu G)^2}. With
#' \eqn{B = \psi'(1) - \psi'(k + 1)} and \eqn{M_2 = 1 + B / A^2}
#' (\eqn{= E[G^2]}),
#' \deqn{E[Y] = \mu, \qquad
#'   \mathrm{Var}[Y] = \mu + \left((1 + \alpha) M_2 - 1\right)\mu^2.}
#' The scale of the generalized exponential distribution cancels in this
#' mean parameterization and is not a parameter. \eqn{\alpha = 0} gives the
#' Poisson-generalized exponential distribution and, as the
#' generalized-exponential heterogeneity vanishes, the distribution
#' approaches NB2. The two dispersion mechanisms have confounded effects on
#' the variance, so \eqn{k} and \eqn{\alpha} are only separately identified
#' through the shape of the distribution; expect weak identification in
#' small samples.
#'
#' This construction is specified here; a distribution with a similar name
#' elsewhere may be defined differently.
#'
#' Probabilities are computed on the log scale by adaptive quadrature over
#' the logarithm of the mixing variable (the same deterministic,
#' mode-centred rule as \code{\link{dpge}}), so tail probabilities keep
#' their relative accuracy. Upper tails are summed directly from the
#' probabilities. Random deviates are generated from the definition.
#'
#' @returns \code{dnbge} gives the probability mass, \code{pnbge} the
#'   distribution function, \code{qnbge} the quantile function, and
#'   \code{rnbge} random deviates. The length of the result is determined by
#'   \code{n} for \code{rnbge}, and is the maximum of the lengths of the
#'   numerical arguments for the other functions.
#'
#' @references
#' Gupta, R. D. and Kundu, D. (1999). Generalized exponential distributions.
#' \emph{Australian & New Zealand Journal of Statistics}, 41(2), 173-188.
#'
#' Greene, W. (2008). Functional forms for the negative binomial model for
#' count data. \emph{Economics Letters}, 99(3), 585-590.
#'
#' @examples
#' dnbge(0:6, mu = 2, shape = 1.5, alpha = 0.5)
#' # alpha = 0 gives the Poisson-generalized exponential distribution
#' all.equal(dnbge(0:10, mu = 2, shape = 1.5, alpha = 0),
#'           dpge(0:10, mean = 2, shape = 1.5, scale = 1))
#' pnbge(5, mu = 2, shape = 1.5, alpha = 0.5)
#' qnbge(c(0.5, 0.95), mu = 2, shape = 1.5, alpha = 0.5)
#' set.seed(8)
#' y <- rnbge(5000, mu = 2, shape = 1.5, alpha = 0.5)
#' c(mean = mean(y), var = var(y))
#'
#' @name NBGE
#' @rdname NBGE
#' @export
dnbge <- function(x, mu, shape, alpha, log = FALSE) {
  .fcr_d(x, list(mu = mu, shape = shape, alpha = alpha), .nbge_valid,
         .nbge_logpmf, log)
}

#' @rdname NBGE
#' @export
pnbge <- function(q, mu, shape, alpha, lower.tail = TRUE, log.p = FALSE) {
  .fcr_p(q, list(mu = mu, shape = shape, alpha = alpha), .nbge_valid,
         .nbge_logcdf, lower.tail, log.p)
}

#' @rdname NBGE
#' @export
qnbge <- function(p, mu, shape, alpha, lower.tail = TRUE, log.p = FALSE) {
  .fcr_q(p, list(mu = mu, shape = shape, alpha = alpha), .nbge_valid,
         .nbge_logcdf, .nbge_qstart, lower.tail, log.p)
}

#' @rdname NBGE
#' @export
rnbge <- function(n, mu, shape, alpha) {
  .fcr_r(n, list(mu = mu, shape = shape, alpha = alpha), .nbge_valid,
         .nbge_rand)
}
