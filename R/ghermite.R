# Generalized Hermite distribution of order m, mean parameterized

.ghermite_valid <- function(p) {
  is.finite(p$mu) & p$mu >= 0 & is.finite(p$rho) & p$rho >= 0 &
    p$rho <= 1 & is.finite(p$m) & p$m >= 1 & p$m == round(p$m)
}

# Y = N1 + m Nm, N1 ~ Poisson((1 - rho) mu), Nm ~ Poisson(rho mu / m)
.ghermite_logpmf <- function(y, p) {
  herm_log_cpp(y, (1 - p$rho) * p$mu, p$rho * p$mu / p$m, p$m, 0L)
}

.ghermite_logcdf <- function(y, p, lower) {
  herm_log_cpp(y, (1 - p$rho) * p$mu, p$rho * p$mu / p$m, p$m,
               if (lower) 1L else 2L)
}

.ghermite_moments <- function(mu, rho, m = 2) {
  list(mean = mu, var = mu * (1 + (m - 1) * rho))
}

.ghermite_qstart <- function(p, llow) {
  mo <- .ghermite_moments(p$mu, p$rho, p$m)
  .fcr_qstart_normal(mo$mean, mo$var, llow)
}

.ghermite_rand <- function(n, p) {
  stats::rpois(n, (1 - p$rho) * p$mu) +
    p$m * stats::rpois(n, p$rho * p$mu / p$m)
}

#' Generalized Hermite Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the (generalized) Hermite distribution of order \eqn{m}, in a mean
#' parameterization.
#'
#' @inheritParams DiscreteWeibull
#' @param mu vector of means \eqn{\mu \geq 0}.
#' @param rho vector of shares \eqn{0 \leq \rho \leq 1} of the mean that come
#'   from events occurring in batches of size \eqn{m}.
#' @param m order: the (fixed, whole) batch size, \eqn{m \geq 1}. The
#'   default \code{m = 2} is the classical Hermite distribution.
#'
#' @details
#' A generalized Hermite random variable counts single events and events
#' that occur in batches of \eqn{m}:
#' \deqn{Y = N_1 + m N_m, \quad N_1 \sim \mathrm{Poisson}((1 - \rho)\mu),
#'   \quad N_m \sim \mathrm{Poisson}(\rho\mu / m),}
#' independently, so that
#' \deqn{E[Y] = \mu, \qquad \mathrm{Var}[Y] = \mu(1 + (m - 1)\rho).}
#' The variance-to-mean ratio (dispersion index) \eqn{d = 1 + (m - 1)\rho}
#' lies between 1 and \eqn{m}: the distribution is overdispersed but with
#' bounded dispersion, and it can be multimodal. \eqn{\rho = 0} (or
#' \eqn{m = 1}) gives the Poisson distribution; with \eqn{\rho = 1} all of
#' the mass is on multiples of \eqn{m}, so other counts have probability
#' zero. In the \pkg{hermite} package parameterization (order \eqn{m},
#' mean, dispersion index) the dispersion index is \eqn{d}.
#'
#' \eqn{m} is a fixed model choice, not a continuous parameter. Probabilities
#' are computed by summing the convolution of the two Poisson components on
#' the log scale; both tails are computed directly. Random deviates are
#' generated from the definition.
#'
#' @returns \code{dghermite} gives the probability mass, \code{pghermite}
#'   the distribution function, \code{qghermite} the quantile function, and
#'   \code{rghermite} random deviates. The length of the result is
#'   determined by \code{n} for \code{rghermite}, and is the maximum of the
#'   lengths of the numerical arguments for the other functions.
#'
#' @references
#' Kemp, C. D. and Kemp, A. W. (1965). Some properties of the "Hermite"
#' distribution. \emph{Biometrika}, 52(3/4), 381-394.
#'
#' Puig, P. (2003). Characterizing additively closed discrete models by a
#' property of their maximum likelihood estimators, with an application to
#' generalized Hermite distributions. \emph{Journal of the American
#' Statistical Association}, 98(463), 687-692.
#'
#' Morina, D., Higueras, M., Puig, P. and Oliveira, M. (2015). Generalized
#' Hermite distribution modelling with the R package hermite. \emph{The R
#' Journal}, 7(2), 263-274.
#'
#' @examples
#' dghermite(0:6, mu = 2, rho = 0.4)
#' # rho = 0 is Poisson
#' all.equal(dghermite(0:10, mu = 2, rho = 0), dpois(0:10, 2))
#' # rho = 1: only multiples of m
#' dghermite(0:6, mu = 2, rho = 1, m = 3)
#' pghermite(4, mu = 2, rho = 0.4, m = 3)
#' qghermite(c(0.5, 0.95), mu = 2, rho = 0.4)
#' set.seed(4)
#' y <- rghermite(5000, mu = 3, rho = 0.5, m = 3)
#' c(mean = mean(y), var = var(y))   # theory: 3 and 3 * (1 + 2 * 0.5) = 6
#'
#' @name GeneralizedHermite
#' @rdname GeneralizedHermite
#' @export
dghermite <- function(x, mu, rho, m = 2, log = FALSE) {
  .fcr_d(x, list(mu = mu, rho = rho, m = m), .ghermite_valid,
         .ghermite_logpmf, log)
}

#' @rdname GeneralizedHermite
#' @export
pghermite <- function(q, mu, rho, m = 2, lower.tail = TRUE, log.p = FALSE) {
  .fcr_p(q, list(mu = mu, rho = rho, m = m), .ghermite_valid,
         .ghermite_logcdf, lower.tail, log.p)
}

#' @rdname GeneralizedHermite
#' @export
qghermite <- function(p, mu, rho, m = 2, lower.tail = TRUE, log.p = FALSE) {
  .fcr_q(p, list(mu = mu, rho = rho, m = m), .ghermite_valid,
         .ghermite_logcdf, .ghermite_qstart, lower.tail, log.p)
}

#' @rdname GeneralizedHermite
#' @export
rghermite <- function(n, mu, rho, m = 2) {
  .fcr_r(n, list(mu = mu, rho = rho, m = m), .ghermite_valid,
         .ghermite_rand)
}
