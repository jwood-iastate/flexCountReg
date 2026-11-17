# Hyper-Poisson distribution (Bardwell and Crow, 1964)

.hpois_valid <- function(p) {
  is.finite(p$lambda) & p$lambda > 0 & is.finite(p$gamma) & p$gamma > 0
}

# Unnormalized log weights: lambda^y / (gamma)_y, (gamma)_y the rising
# factorial Gamma(gamma + y) / Gamma(gamma)
.hpois_logw <- function(y, lambda, gamma) {
  y * log(lambda) - (lgamma(gamma + y) - lgamma(gamma))
}

.hpois_k0 <- function(lambda, gamma) {
  max(lambda - gamma, 0) + 10 * sqrt(lambda + 1) + 30
}

.hpois_logpmf <- function(y, p) {
  pois <- p$gamma == 1
  out <- numeric(length(y))
  out[pois] <- stats::dpois(y[pois], p$lambda[pois], log = TRUE)
  if (any(!pois)) {
    out[!pois] <- .fcr_wtable_logpmf(y[!pois], .fcr_subset_pars(p, !pois),
                                     .hpois_logw, .hpois_k0)
  }
  out
}

.hpois_logcdf <- function(y, p, lower) {
  pois <- p$gamma == 1
  out <- numeric(length(y))
  out[pois] <- stats::ppois(y[pois], p$lambda[pois], lower.tail = lower,
                            log.p = TRUE)
  if (any(!pois)) {
    out[!pois] <- .fcr_wtable_logcdf(y[!pois], .fcr_subset_pars(p, !pois),
                                     lower, .hpois_logw, .hpois_k0)
  }
  out
}

.hpois_qstart <- function(p, llow) {
  .fcr_qstart_normal(pmax(p$lambda - p$gamma + 1, 0), p$lambda, llow)
}

.hpois_rand <- function(n, p) {
  .fcr_q(stats::runif(n), p, .hpois_valid, .hpois_logcdf, .hpois_qstart,
         TRUE, FALSE)
}

.hpois_moments <- function(lambda, gamma) {
  r <- .fcr_recycle(lambda, gamma)
  p <- list(lambda = r[[1]], gamma = r[[2]])
  out <- list(mean = rep(NA_real_, length(r[[1]])),
              var = rep(NA_real_, length(r[[1]])))
  ok <- .hpois_valid(p)
  ok[is.na(ok)] <- FALSE
  bad <- !ok & !.fcr_any_na(p, length(ok))
  out$mean[bad] <- NaN
  out$var[bad] <- NaN
  if (any(ok)) {
    m <- .fcr_wtable_moments(.fcr_subset_pars(p, ok), .hpois_logw,
                             .hpois_k0)
    out$mean[ok] <- m$mean
    out$var[ok] <- m$var
  }
  out
}

#' Hyper-Poisson Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the hyper-Poisson distribution of Bardwell and Crow (1964).
#'
#' @inheritParams DiscreteWeibull
#' @param lambda vector of (positive) rate parameters \eqn{\lambda}.
#' @param gamma vector of (positive) dispersion parameters \eqn{\gamma}.
#'   \eqn{\gamma = 1} is the Poisson distribution, \eqn{\gamma > 1} gives
#'   overdispersion and \eqn{\gamma < 1} underdispersion.
#'
#' @details
#' The hyper-Poisson distribution has probabilities
#' \deqn{P(Y = y) = \frac{\lambda^y}{(\gamma)_y \,
#'   {}_1F_1(1; \gamma; \lambda)}, \quad y = 0, 1, 2, \ldots}
#' where \eqn{(\gamma)_y = \Gamma(\gamma + y) / \Gamma(\gamma)} is the
#' rising factorial and \eqn{{}_1F_1(1; \gamma; \lambda) =
#' \sum_{j \geq 0} \lambda^j / (\gamma)_j} the confluent hypergeometric
#' function, the normalizing constant. Successive probabilities satisfy
#' \eqn{P(Y = y + 1) / P(Y = y) = \lambda / (\gamma + y)}.
#'
#' \eqn{\lambda} is not the mean. The mean is
#' \deqn{E[Y] = \lambda - (\gamma - 1)(1 - P(Y = 0)),}
#' and moments are computed from the normalized probabilities when needed.
#' Regression models usually reparameterize the distribution by its mean
#' (Saez-Castillo and Conde-Sanchez, 2013); symbols in other packages can
#' differ from those used here.
#'
#' The normalizing constant is computed by summing the series around its
#' mode until the remaining tail (bounded by a geometric series, because the
#' ratio of successive terms decreases) is negligible, so the computations
#' do not overflow for large \eqn{\lambda}. Upper tail probabilities are
#' summed directly. \code{rhpois()} uses inversion of the distribution
#' function.
#'
#' @returns \code{dhpois} gives the probability mass, \code{phpois} the
#'   distribution function, \code{qhpois} the quantile function, and
#'   \code{rhpois} random deviates. The length of the result is determined
#'   by \code{n} for \code{rhpois}, and is the maximum of the lengths of the
#'   numerical arguments for the other functions.
#'
#' @references
#' Bardwell, G. E. and Crow, E. L. (1964). A two-parameter family of
#' hyper-Poisson distributions. \emph{Journal of the American Statistical
#' Association}, 59(305), 133-141.
#'
#' Saez-Castillo, A. J. and Conde-Sanchez, A. (2013). A hyper-Poisson
#' regression model for overdispersed and underdispersed count data.
#' \emph{Computational Statistics & Data Analysis}, 61, 148-157.
#'
#' @examples
#' dhpois(0:6, lambda = 3, gamma = 2)
#' # gamma = 1 gives the Poisson distribution
#' all.equal(dhpois(0:10, lambda = 3, gamma = 1), dpois(0:10, 3))
#' phpois(4, lambda = 3, gamma = 0.5)
#' qhpois(c(0.1, 0.5, 0.9), lambda = 3, gamma = 2)
#' set.seed(5)
#' y <- rhpois(2000, lambda = 4, gamma = 0.3)
#' c(mean = mean(y), var = var(y))   # underdispersed
#'
#' @name HyperPoisson
#' @rdname HyperPoisson
#' @export
dhpois <- function(x, lambda, gamma, log = FALSE) {
  .fcr_d(x, list(lambda = lambda, gamma = gamma), .hpois_valid,
         .hpois_logpmf, log)
}

#' @rdname HyperPoisson
#' @export
phpois <- function(q, lambda, gamma, lower.tail = TRUE, log.p = FALSE) {
  .fcr_p(q, list(lambda = lambda, gamma = gamma), .hpois_valid,
         .hpois_logcdf, lower.tail, log.p)
}

#' @rdname HyperPoisson
#' @export
qhpois <- function(p, lambda, gamma, lower.tail = TRUE, log.p = FALSE) {
  .fcr_q(p, list(lambda = lambda, gamma = gamma), .hpois_valid,
         .hpois_logcdf, .hpois_qstart, lower.tail, log.p)
}

#' @rdname HyperPoisson
#' @export
rhpois <- function(n, lambda, gamma) {
  .fcr_r(n, list(lambda = lambda, gamma = gamma), .hpois_valid,
         .hpois_rand)
}
