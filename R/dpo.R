# Double Poisson (Efron, 1986), exactly normalized

.dpo_valid <- function(p) {
  is.finite(p$lambda) & p$lambda > 0 & is.finite(p$phi) & p$phi > 0
}

# Unnormalized log weights, 0 log 0 = 0. Terms constant in y are omitted
# (normalization removes them); phi = 1 gives y log(lambda) - lgamma(y + 1).
.dpo_logw <- function(y, lambda, phi) {
  ylogy <- ifelse(y > 0, y * log(pmax(y, 1)), 0)
  phi * y * log(lambda) + (1 - phi) * ylogy + (phi - 1) * y - lgamma(y + 1)
}

# Initial table size: about the location plus 10 approximate SDs
.dpo_k0 <- function(lambda, phi) lambda + 10 * sqrt(lambda / phi + 1) + 30

.dpo_logpmf <- function(y, p) {
  pois <- p$phi == 1
  out <- numeric(length(y))
  out[pois] <- stats::dpois(y[pois], p$lambda[pois], log = TRUE)
  if (any(!pois)) {
    out[!pois] <- .fcr_wtable_logpmf(y[!pois], .fcr_subset_pars(p, !pois),
                                     .dpo_logw, .dpo_k0)
  }
  out
}

.dpo_logcdf <- function(y, p, lower) {
  pois <- p$phi == 1
  out <- numeric(length(y))
  out[pois] <- stats::ppois(y[pois], p$lambda[pois], lower.tail = lower,
                            log.p = TRUE)
  if (any(!pois)) {
    out[!pois] <- .fcr_wtable_logcdf(y[!pois], .fcr_subset_pars(p, !pois),
                                     lower, .dpo_logw, .dpo_k0)
  }
  out
}

.dpo_qstart <- function(p, llow) {
  .fcr_qstart_normal(p$lambda, p$lambda / p$phi, llow)
}

.dpo_rand <- function(n, p) {
  .fcr_q(stats::runif(n), p, .dpo_valid, .dpo_logcdf, .dpo_qstart,
         TRUE, FALSE)
}

# Exact mean and variance of the normalized distribution
.dpo_moments <- function(lambda, phi) {
  r <- .fcr_recycle(lambda, phi)
  p <- list(lambda = r[[1]], phi = r[[2]])
  out <- list(mean = rep(NA_real_, length(r[[1]])),
              var = rep(NA_real_, length(r[[1]])))
  ok <- .dpo_valid(p)
  ok[is.na(ok)] <- FALSE
  bad <- !ok & !.fcr_any_na(p, length(ok))
  out$mean[bad] <- NaN
  out$var[bad] <- NaN
  if (any(ok)) {
    m <- .fcr_wtable_moments(.fcr_subset_pars(p, ok), .dpo_logw, .dpo_k0)
    out$mean[ok] <- m$mean
    out$var[ok] <- m$var
  }
  out
}

#' Double Poisson Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the double Poisson distribution of Efron (1986), normalized exactly.
#'
#' @inheritParams DiscreteWeibull
#' @param lambda vector of (positive) location parameters \eqn{\lambda}.
#' @param phi vector of (positive) precision parameters \eqn{\phi}.
#'   \eqn{\phi = 1} is the Poisson distribution, \eqn{\phi < 1} gives
#'   overdispersion and \eqn{\phi > 1} underdispersion.
#'
#' @details
#' The double Poisson distribution has probabilities proportional to the
#' weights (with \eqn{0 \log 0 = 0})
#' \deqn{w(y) = \exp\left\{\phi y \log\lambda + (1 - \phi) y \log y +
#'   (\phi - 1) y - \log y!\right\}, \quad y = 0, 1, 2, \ldots}
#' which is Efron's (1986) expression with the factors that do not depend on
#' \eqn{y} removed. Here the normalizing constant
#' \eqn{\sum_{j \geq 0} w(j)} is computed numerically (by summing the
#' weights around their mode until the remaining tail, bounded by a
#' geometric series, is negligible), so the probabilities sum to one. Efron's
#' closed-form approximation of the constant is not used: it is least
#' accurate for small \eqn{\lambda} and strong dispersion.
#'
#' \eqn{\lambda} is a location parameter: the mean is close to, but not
#' exactly, \eqn{\lambda}, and the variance is approximately
#' \eqn{\lambda / \phi}. The exact moments are computed from the normalized
#' probabilities when needed. Some references (e.g. the \code{DPO} family of
#' \pkg{gamlss.dist}) use the dispersion \eqn{\sigma = 1 / \phi} instead of
#' the precision.
#'
#' Upper tail probabilities (\code{lower.tail = FALSE}) are summed directly,
#' never computed as one minus the distribution function. \code{rdpo()}
#' uses inversion of the distribution function.
#'
#' @returns \code{ddpo} gives the probability mass, \code{pdpo} the
#'   distribution function, \code{qdpo} the quantile function, and
#'   \code{rdpo} random deviates. The length of the result is determined by
#'   \code{n} for \code{rdpo}, and is the maximum of the lengths of the
#'   numerical arguments for the other functions.
#'
#' @references
#' Efron, B. (1986). Double exponential families and their use in
#' generalized linear regression. \emph{Journal of the American Statistical
#' Association}, 81(395), 709-721.
#'
#' @examples
#' ddpo(0:6, lambda = 3, phi = 0.5)
#' # phi = 1 gives the Poisson distribution
#' all.equal(ddpo(0:10, lambda = 3, phi = 1), dpois(0:10, 3))
#' sum(ddpo(0:200, lambda = 3, phi = 2.5))   # exactly normalized
#' pdpo(c(2, 5), lambda = 3, phi = 2.5)
#' qdpo(c(0.1, 0.5, 0.9), lambda = 3, phi = 0.5)
#' set.seed(2)
#' y <- rdpo(2000, lambda = 5, phi = 3)
#' c(mean = mean(y), var = var(y))   # underdispersed
#'
#' @name DoublePoisson
#' @rdname DoublePoisson
#' @export
ddpo <- function(x, lambda, phi, log = FALSE) {
  .fcr_d(x, list(lambda = lambda, phi = phi), .dpo_valid, .dpo_logpmf, log)
}

#' @rdname DoublePoisson
#' @export
pdpo <- function(q, lambda, phi, lower.tail = TRUE, log.p = FALSE) {
  .fcr_p(q, list(lambda = lambda, phi = phi), .dpo_valid, .dpo_logcdf,
         lower.tail, log.p)
}

#' @rdname DoublePoisson
#' @export
qdpo <- function(p, lambda, phi, lower.tail = TRUE, log.p = FALSE) {
  .fcr_q(p, list(lambda = lambda, phi = phi), .dpo_valid, .dpo_logcdf,
         .dpo_qstart, lower.tail, log.p)
}

#' @rdname DoublePoisson
#' @export
rdpo <- function(n, lambda, phi) {
  .fcr_r(n, list(lambda = lambda, phi = phi), .dpo_valid, .dpo_rand)
}
