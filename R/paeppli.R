# Polya-Aeppli (Poisson-geometric) distribution

.pa_valid <- function(p) {
  is.finite(p$lambda) & p$lambda >= 0 & is.finite(p$prob) & p$prob >= 0 &
    p$prob < 1
}

.pa_logpmf <- function(y, p) pa_log_cpp(y, p$lambda, p$prob, 0L)

.pa_logcdf <- function(y, p, lower) {
  pa_log_cpp(y, p$lambda, p$prob, if (lower) 1L else 2L)
}

.pa_moments <- function(lambda, prob) {
  list(mean = lambda / (1 - prob),
       var = lambda * (1 + prob) / (1 - prob)^2)
}

.pa_qstart <- function(p, llow) {
  m <- .pa_moments(p$lambda, p$prob)
  .fcr_qstart_normal(m$mean, m$var, llow)
}

# M clusters ~ Poisson(lambda); given M, Y - M ~ NB(size = M, prob = 1 - q)
.pa_rand <- function(n, p) {
  M <- stats::rpois(n, p$lambda)
  out <- M
  pos <- M > 0   # rnbinom() with size 0 is not used: it can return NA
  out[pos] <- M[pos] + stats::rnbinom(sum(pos), size = M[pos],
                                      prob = 1 - p$prob[pos])
  out
}

# Poisson-Tweedie distribution with power 1.5 (mean mu, Var = mu + phi
# mu^1.5) as a Polya-Aeppli distribution. The Tweedie intensity with power
# 1.5 is compound Poisson with event rate 2 sqrt(mu) / phi and exponential
# jumps of mean theta = phi sqrt(mu) / 2; its Laplace transform in the
# Poisson probability generating function gives the Polya-Aeppli PGF
# exp(lambda (z - 1) / (1 - q z)) with q = theta / (1 + theta) and
# lambda = mu (1 - q). Used by the regression family in a later release.
.pt15_to_paeppli <- function(mu, phi) {
  d <- phi * sqrt(mu)
  q <- d / (2 + d)
  list(lambda = mu * (1 - q), prob = q)
}

#' Polya-Aeppli Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Polya-Aeppli (Poisson-geometric) distribution.
#'
#' @inheritParams DiscreteWeibull
#' @param lambda vector of (non-negative) Poisson rates \eqn{\lambda} of
#'   clusters.
#' @param prob vector of geometric parameters \eqn{0 \leq q < 1}: each
#'   cluster has \eqn{c \geq 1} events with probability
#'   \eqn{(1 - q) q^{c - 1}}.
#'
#' @details
#' A Polya-Aeppli random variable is a Poisson number of clusters with
#' geometric cluster sizes:
#' \deqn{Y = \sum_{i = 1}^{M} C_i, \quad M \sim \mathrm{Poisson}(\lambda),
#'   \quad P(C_i = c) = (1 - q) q^{c - 1}, \; c = 1, 2, \ldots}
#' so that \eqn{P(Y = 0) = e^{-\lambda}} and, for \eqn{y \geq 1},
#' \deqn{P(Y = y) = e^{-\lambda} \sum_{j = 1}^{y} \frac{\lambda^j}{j!}
#'   \binom{y - 1}{j - 1} q^{y - j} (1 - q)^j.}
#' The mean and variance are
#' \deqn{E[Y] = \frac{\lambda}{1 - q}, \qquad
#'   \mathrm{Var}[Y] = \frac{\lambda (1 + q)}{(1 - q)^2}.}
#' \eqn{q = 0} gives the Poisson distribution with mean \eqn{\lambda}. The
#' parameterization is that of the \pkg{polyaAeppli} package.
#'
#' The Polya-Aeppli distribution is the Poisson-Tweedie distribution with
#' power parameter 1.5: a count with mean \eqn{\mu} and variance
#' \eqn{\mu + \phi\mu^{1.5}} is Polya-Aeppli with \eqn{q = d / (2 + d)} and
#' \eqn{\lambda = \mu (1 - q)}, where \eqn{d = \phi\sqrt{\mu}}.
#'
#' Probabilities are computed by summing the decomposition over the number
#' of clusters on the log scale (all terms are positive); both tails are
#' computed directly. Random deviates are generated from the definition.
#'
#' @returns \code{dpaeppli} gives the probability mass, \code{ppaeppli} the
#'   distribution function, \code{qpaeppli} the quantile function, and
#'   \code{rpaeppli} random deviates. The length of the result is
#'   determined by \code{n} for \code{rpaeppli}, and is the maximum of the
#'   lengths of the numerical arguments for the other functions.
#'
#' @references
#' Johnson, N. L., Kemp, A. W. and Kotz, S. (2005). \emph{Univariate
#' Discrete Distributions} (3rd ed.), Section 9.7. Wiley.
#'
#' Bonat, W. H., Jorgensen, B., Kokonendji, C. C., Hinde, J. and Demetrio,
#' C. G. B. (2018). Extended Poisson-Tweedie: properties and regression
#' models for count data. \emph{Statistical Modelling}, 18(1), 24-49.
#'
#' @examples
#' dpaeppli(0:6, lambda = 2, prob = 0.3)
#' # prob = 0 gives the Poisson distribution
#' all.equal(dpaeppli(0:10, lambda = 2, prob = 0), dpois(0:10, 2))
#' ppaeppli(4, lambda = 2, prob = 0.3)
#' qpaeppli(c(0.5, 0.95), lambda = 2, prob = 0.3)
#' set.seed(7)
#' y <- rpaeppli(5000, lambda = 2, prob = 0.3)
#' c(mean = mean(y), var = var(y))   # theory: 2.857 and 5.306
#'
#' @name PolyaAeppli
#' @rdname PolyaAeppli
#' @export
dpaeppli <- function(x, lambda, prob, log = FALSE) {
  .fcr_d(x, list(lambda = lambda, prob = prob), .pa_valid, .pa_logpmf, log)
}

#' @rdname PolyaAeppli
#' @export
ppaeppli <- function(q, lambda, prob, lower.tail = TRUE, log.p = FALSE) {
  .fcr_p(q, list(lambda = lambda, prob = prob), .pa_valid, .pa_logcdf,
         lower.tail, log.p)
}

#' @rdname PolyaAeppli
#' @export
qpaeppli <- function(p, lambda, prob, lower.tail = TRUE, log.p = FALSE) {
  .fcr_q(p, list(lambda = lambda, prob = prob), .pa_valid, .pa_logcdf,
         .pa_qstart, lower.tail, log.p)
}

#' @rdname PolyaAeppli
#' @export
rpaeppli <- function(n, lambda, prob) {
  .fcr_r(n, list(lambda = lambda, prob = prob), .pa_valid, .pa_rand)
}
