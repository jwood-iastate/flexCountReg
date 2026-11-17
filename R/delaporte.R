# Delaporte distribution: Poisson + negative binomial, mean parameterized

.delap_valid <- function(p) {
  is.finite(p$mu) & p$mu >= 0 & is.finite(p$rho) & p$rho >= 0 &
    p$rho <= 1 & is.finite(p$k) & p$k > 0
}

# Components: A ~ Poisson((1 - rho) mu), B ~ NB(size = k, mean = rho mu)
.delap_logpmf <- function(y, p) {
  delap_log_cpp(y, (1 - p$rho) * p$mu, p$k, p$rho * p$mu, 0L)
}

.delap_logcdf <- function(y, p, lower) {
  delap_log_cpp(y, (1 - p$rho) * p$mu, p$k, p$rho * p$mu,
                if (lower) 1L else 2L)
}

.delap_moments <- function(mu, rho, k) {
  list(mean = mu, var = mu + (rho * mu)^2 / k)
}

.delap_qstart <- function(p, llow) {
  m <- .delap_moments(p$mu, p$rho, p$k)
  .fcr_qstart_normal(m$mean, m$var, llow)
}

.delap_rand <- function(n, p) {
  stats::rpois(n, (1 - p$rho) * p$mu) +
    stats::rnbinom(n, size = p$k, mu = p$rho * p$mu)
}

#' Delaporte Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the Delaporte distribution in a mean parameterization.
#'
#' @inheritParams DiscreteWeibull
#' @param mu vector of means \eqn{\mu \geq 0}.
#' @param rho vector of shares \eqn{0 \leq \rho \leq 1} of the mean that come
#'   from the negative binomial (gamma-mixed) component.
#' @param k vector of (positive) negative binomial size (shape) parameters.
#'
#' @details
#' A Delaporte random variable is the sum of two independent counts,
#' \deqn{Y = A + B, \quad A \sim \mathrm{Poisson}((1 - \rho)\mu), \quad
#'   B \sim \mathrm{NB}(\mathrm{size} = k, \mathrm{mean} = \rho\mu),}
#' so that
#' \deqn{E[Y] = \mu, \qquad \mathrm{Var}[Y] = \mu + \frac{(\rho\mu)^2}{k}.}
#' Equivalently it is a Poisson distribution whose rate is a shifted gamma
#' variable. \eqn{\rho = 0} gives the Poisson distribution with mean
#' \eqn{\mu} and \eqn{\rho = 1} the negative binomial distribution with mean
#' \eqn{\mu} and size \eqn{k} (NB2 with \eqn{\alpha = 1/k}). Both endpoints
#' are accepted here.
#'
#' In the \eqn{(\alpha, \beta, \lambda)} parameterization of the
#' \pkg{Delaporte} package, \eqn{\alpha = k}, \eqn{\beta = \rho\mu / k} and
#' \eqn{\lambda = (1 - \rho)\mu}.
#'
#' Probabilities are computed by summing the convolution of the two
#' components on the log scale (all terms are positive). Both tails are
#' computed directly from the same decomposition, so
#' \code{lower.tail = FALSE} is accurate far into the upper tail. Random
#' deviates are generated as the sum of independent Poisson and negative
#' binomial draws. The work per probability grows linearly with the count.
#'
#' @returns \code{ddelaporte} gives the probability mass, \code{pdelaporte}
#'   the distribution function, \code{qdelaporte} the quantile function, and
#'   \code{rdelaporte} random deviates. The length of the result is
#'   determined by \code{n} for \code{rdelaporte}, and is the maximum of the
#'   lengths of the numerical arguments for the other functions.
#'
#' @references
#' Ruohonen, M. (1988). On a model for the claim number process.
#' \emph{ASTIN Bulletin}, 18(1), 57-68.
#'
#' Vose, D. (2008). \emph{Risk Analysis: A Quantitative Guide} (3rd ed.).
#' Wiley.
#'
#' @examples
#' ddelaporte(0:6, mu = 2, rho = 0.5, k = 1.5)
#' # rho = 0 is Poisson, rho = 1 is negative binomial
#' all.equal(ddelaporte(0:10, mu = 2, rho = 0, k = 1.5), dpois(0:10, 2))
#' all.equal(ddelaporte(0:10, mu = 2, rho = 1, k = 1.5),
#'           dnbinom(0:10, size = 1.5, mu = 2))
#' pdelaporte(5, mu = 2, rho = 0.5, k = 1.5)
#' qdelaporte(c(0.5, 0.95), mu = 2, rho = 0.5, k = 1.5)
#' set.seed(3)
#' y <- rdelaporte(5000, mu = 4, rho = 0.6, k = 2)
#' c(mean = mean(y), var = var(y))   # theory: 4 and 4 + 2.4^2 / 2 = 6.88
#'
#' @name Delaporte
#' @rdname Delaporte
#' @export
ddelaporte <- function(x, mu, rho, k, log = FALSE) {
  .fcr_d(x, list(mu = mu, rho = rho, k = k), .delap_valid, .delap_logpmf,
         log)
}

#' @rdname Delaporte
#' @export
pdelaporte <- function(q, mu, rho, k, lower.tail = TRUE, log.p = FALSE) {
  .fcr_p(q, list(mu = mu, rho = rho, k = k), .delap_valid, .delap_logcdf,
         lower.tail, log.p)
}

#' @rdname Delaporte
#' @export
qdelaporte <- function(p, mu, rho, k, lower.tail = TRUE, log.p = FALSE) {
  .fcr_q(p, list(mu = mu, rho = rho, k = k), .delap_valid, .delap_logcdf,
         .delap_qstart, lower.tail, log.p)
}

#' @rdname Delaporte
#' @export
rdelaporte <- function(n, mu, rho, k) {
  .fcr_r(n, list(mu = mu, rho = rho, k = k), .delap_valid, .delap_rand)
}
