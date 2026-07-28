#' Raw Moment of a Lognormal Distribution
#'
#' Computes a raw moment of a lognormal distribution using its closed-form
#' expression. The lognormal distribution is specified by its log-mean
#' (\eqn{\mu}) and log-standard deviation (\eqn{\sigma}).
#'
#' @param mu Numeric vector of means of the log-transformed variable,
#'   corresponding to \eqn{\mu}.
#' @param sigma Non-negative numeric vector of standard deviations of the
#'   log-transformed variable, corresponding to \eqn{\sigma}.
#' @param n Numeric vector specifying the order of the raw moment.
#'
#' @returns A numeric vector containing the raw moments of the specified
#'   lognormal distribution.
#'
#' @details
#' Let \eqn{X \sim N(\mu, \sigma^2)} and \eqn{Y = \exp(X)}. Then \eqn{Y}
#' follows a lognormal distribution, and its raw moment of order \eqn{n} is:
#' \deqn{
#' E[Y^n] =
#' \exp\left(n\mu + \frac{n^2\sigma^2}{2}\right).
#' }
#'
#' In particular, the mean is:
#' \deqn{
#' E[Y] =
#' \exp\left(\mu + \frac{\sigma^2}{2}\right),
#' }
#' and the variance can be calculated as:
#' \deqn{
#' \mathrm{Var}(Y) =
#' E[Y^2] - E[Y]^2.
#' }
#'
#' This function can be used to adjust predictions from generalized linear
#' mixed models with normally distributed random parameters and a log link.
#'
#' @examples
#' mu <- 0
#' sigma <- 1
#'
#' # Mean of the lognormal distribution
#' moment_lognormal(mu, sigma, n = 1)
#'
#' # Variance of the lognormal distribution
#' moment_1 <- moment_lognormal(mu, sigma, n = 1)
#' moment_2 <- moment_lognormal(mu, sigma, n = 2)
#' moment_2 - moment_1^2
#'
#' @export
#' @name moment_lognormal
moment_lognormal <- function(mu, sigma, n) {
  if (!is.numeric(mu) || !is.numeric(sigma) || !is.numeric(n)) {
    stop("`mu`, `sigma`, and `n` must be numeric.", call. = FALSE)
  }
  
  if (any(sigma < 0, na.rm = TRUE)) {
    stop("`sigma` must be non-negative.", call. = FALSE)
  }
  
  exp(n * mu + n^2 * sigma^2 / 2)
}