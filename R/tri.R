#' Triangle Distribution
#'
#' These functions provide density, distribution function, quantile function,
#' and random number generation for the Triangle Distribution, specified by its
#' mean, standard deviation, and optional lower and upper bounds.
#'
#' The Triangle Distribution is defined by three points: a (minimum), b
#' (maximum), and c (mode), where the density is zero outside the interval [a,
#' b], increases linearly from a to c, and decreases linearly from c to b.
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mode numeric value or vector of mode values for the distribution.
#' @param sigma single value or vector indicating both the positive and negative
#'   max differences from the mean (if the difference is the same).
#' @param upper single value or vector for the upper limit of the distribution
#'   (must be used with `lower`).
#' @param lower single value or vector for the lower limit of the distribution
#'   (must be used with `upper`).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dtri} computes the density (PDF) of the Triangle Distribution.
#'
#' \code{ptri} computes the CDF of the Triangle Distribution.
#'
#' \code{qtri} computes the quantile function of the Triangle Distribution.
#'
#' \code{rtri} generates random numbers from the Triangle Distribution.
#'
#' The mode and standard deviation parameters define the distribution's location
#' and scale, respectively, while the lower and upper bounds explicitly set the
#' minimum and maximum values of the distribution.
#'
#' @examples
#' dtri(4, mode=8, upper=13, lower=1)
#' ptri(c(0, 1, 2, 3, 5, 7, 9, 10), mode = 3, upper=9, lower = 1)
#' qtri(c(0.1, 0.3, 0.5, 0.9, 0.95), mode = 3, upper = 9, lower = 1)
#' rtri(30, mode = 5, sigma = 3)
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @export
#' @name Triangular
#' 
#' @rdname Triangular
#' @export
dtri <- Vectorize(function(
    x, mode = 0, sigma = 1, upper = NA, lower = NA, log = FALSE) {
  dtri_cpp(x, mode, sigma, upper, lower, log)
})

#' @rdname Triangular
#' @export
ptri <- Vectorize(function(
    q, mode = 0, sigma = 1, upper = NA, 
    lower = NA, lower.tail = TRUE, log.p = FALSE) {
  ptri_cpp(q, mode, sigma, upper, lower, lower.tail, log.p)
})

#' @rdname Triangular
#' @export
qtri <- Vectorize(function(p, mode = 0, sigma = 1, upper = NA, lower = NA) {
  qtri_cpp(p, mode, sigma, upper, lower)
})

#' @rdname Triangular
#' @export
rtri <- function(n, mode = 0, sigma = 1, upper = NA, lower = NA) {
  rtri_cpp(n, mode, sigma, upper, lower)
}
