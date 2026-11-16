#' Inverse Gamma Distribution
#'
#' These functions provide the density function, distribution function, quantile
#' function, and random number generation for the Inverse-Gamma (IG)
#' Distribution
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param shape numeric value or vector of shape values for the distribution
#'   (the values have to be greater than 0).
#' @param scale single value or vector of values for the scale parameter of the
#'   distribution (the values have to be greater than 0).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dinvgamma} computes the density (PDF) of the Inverse-Gamma
#' Distribution.
#'
#' \code{pinvgamma} computes the CDF of the Inverse-Gamma Distribution.
#'
#' \code{qinvgamma} computes the quantile function of the Inverse-Gamma
#' Distribution.
#'
#' \code{rinvgamma} generates random numbers from the Inverse-Gamma
#' Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Inverse-Gamma
#' distribution:
#' \deqn{f(x | \alpha, \beta) = 
#'     \frac{\beta^\alpha}{\Gamma(\alpha)}
#'     \left(\frac{1}{x}\right)^{\alpha+1} e^{-\frac{\beta}{x}}}
#' 
#' Where \eqn{\alpha} is the shape parameter and \eqn{\beta} is a scale
#' parameter with the restrictions that \eqn{\alpha > 0} and \eqn{\eta > 0}, and
#' \eqn{x > 0}.
#' 
#' The CDF of the Inverse-Gamma distribution is:
#' \deqn{F(x | \alpha, \beta) = 
#'    \frac{\alpha. \Gamma \left(\frac{\beta}{x}\right)}{\Gamma(\alpha)} = 
#'    Q\left(\alpha, \frac{\beta}{x} \right)}
#' 
#' Where the numerator is the incomplete gamma function and \eqn{Q(\cdot)} is
#' the regularized gamma function.
#'
#' The mean of the distribution is (provided \eqn{\alpha>1}):
#' \deqn{\mu=\frac{\beta}{\alpha-1}}
#' 
#' The variance of the distribution is (for \eqn{\alpha>2}):
#' \deqn{\sigma^2=\frac{\beta^2}{(\alpha-1)^2(\alpha-2)}}
#' 
#' @returns dinvgamma gives the density, pinvgamma gives the distribution 
#'  function, qinvgamma gives the quantile function, and rinvgamma generates
#'  random  deviates.
#' 
#'  The length of the result is determined by n for rinvgamma, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#' 
#' @examples
#' dinvgamma(1, shape = 3, scale = 2)
#' pinvgamma(c(0.1, 0.5, 1, 3, 5, 10, 30), shape = 3, scale = 2)
#' qinvgamma(c(0.1, 0.3, 0.5, 0.9, 0.95), shape = 3, scale = 2)
#' rinvgamma(30, shape = 3, scale = 2)
#'
#' @importFrom stats dgamma pgamma qgamma runif
#' @export
#' @name invgamma
#' @rdname invgamma
#' @export
dinvgamma <- function(x, shape = 2.5, scale = 1, log = FALSE) {
  
  # --- Input Validation ---
  if (any(shape <= 0, na.rm = TRUE)) {
    warning("'shape' must be positive")
  }
  if (any(scale <= 0, na.rm = TRUE)) {
    warning("'scale' must be positive")
  }
  
  # --- Vectorization Setup ---
  n <- max(length(x), length(shape), length(scale))
  x <- rep_len(x, n)
  shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)
  
  # --- Compute Log-Density ---
  # f(x) = (scale^shape / Gamma(shape)) * x^(-shape-1) * exp(-scale/x)
  # log f(x) = shape*log(scale) - lgamma(shape) - (shape+1)*log(x) - scale/x
  
  log_p <- 
    shape * log(scale) - lgamma(shape) - (shape + 1) * log(x) - scale / x
  
  # Handle invalid x (must be positive), invalid parameters and NAs
  log_p[!is.na(x) & x <= 0] <- -Inf
  log_p[.invgamma_bad(shape, scale)] <- NaN
  log_p[is.na(x) | is.na(shape) | is.na(scale)] <- NA_real_
  
  # --- Return ---
  if (log) {
    return(log_p)
  } else {
    return(exp(log_p))
  }
}


#' @rdname invgamma
#' @export
pinvgamma <- function(q, shape = 2.5, scale = 1, 
                      lower.tail = TRUE, log.p = FALSE) {
  
  # --- Input Validation ---
  if (any(shape <= 0, na.rm = TRUE)) warning("'shape' must be positive")
  if (any(scale <= 0, na.rm = TRUE)) warning("'scale' must be positive")
  
  # --- Vectorization Setup ---
  n <- max(length(q), length(shape), length(scale))
  q <- rep_len(q, n)
  shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)
  
  # --- Compute CDF ---
  # If X ~ InvGamma(shape, scale), then 1/X ~ Gamma(shape, rate = scale)
  # P(X <= q) = P(1/X >= 1/q) 
  ##         = 1 - P(1/X < 1/q) = 1 - pgamma(1/q, shape, rate = scale)
  # Or equivalently: 
  # P(X <= q) = pgamma(1/q, shape, rate = scale, lower.tail = FALSE)
  
  bad <- .invgamma_bad(shape, scale)
  cdf <- rep(NaN, n)
  cdf[!bad] <- pgamma(1/q[!bad], shape = shape[!bad], rate = scale[!bad],
                      lower.tail = FALSE)
  
  # Handle boundaries
  cdf[!bad & !is.na(q) & q <= 0] <- 0
  cdf[!bad & !is.na(q) & is.infinite(q) & q > 0] <- 1
  cdf[is.na(q) | is.na(shape) | is.na(scale)] <- NA_real_
  
  # --- Apply Tail and Log Options ---
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  return(cdf)
}


#' @rdname invgamma
#' @export
qinvgamma <- function(p, shape = 2.5, scale = 1, 
                      lower.tail = TRUE, log.p = FALSE) {
  
  # --- Input Handling ---
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  
  # --- Input Validation ---
  if (any(shape <= 0, na.rm = TRUE)) warning("'shape' must be positive")
  if (any(scale <= 0, na.rm = TRUE)) warning("'scale' must be positive")
  if (any(p < 0 | p > 1, na.rm = TRUE)) warning("'p' must be in [0, 1]")
  
  # --- Vectorization Setup ---
  n <- max(length(p), length(shape), length(scale))
  p <- rep_len(p, n)
  shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)
  
  # --- Compute Quantile ---
  # Q_InvGamma(p) = 1 / Q_Gamma(1-p)
  # where Q_Gamma uses the same shape and rate = scale
  
  ok <- !.invgamma_bad(shape, scale) & !is.na(p) & p >= 0 & p <= 1
  q_val <- rep(NaN, n)
  q_val[ok] <- 1 / qgamma(1 - p[ok], shape = shape[ok], rate = scale[ok])
  
  # Handle boundaries
  q_val[ok & p == 0] <- 0
  q_val[ok & p == 1] <- Inf
  q_val[is.na(p) | is.na(shape) | is.na(scale)] <- NA_real_
  
  return(q_val)
}


#' @rdname invgamma
#' @export
rinvgamma <- function(n, shape = 2.5, scale = 1) {
  n <- .fcr_rn(n)
  if (n == 0) return(numeric(0))
  
  shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)
  bad <- .invgamma_bad(shape, scale) | is.na(shape) | is.na(scale)
  .fcr_warn_na(bad, "'shape' and 'scale' must be positive.")
  # If Y ~ Gamma(shape, rate = scale), then 1/Y ~ InvGamma(shape, scale)
  out <- rep(NA_real_, n)
  out[!bad] <- 1 / rgamma(sum(!bad), shape = shape[!bad], rate = scale[!bad])
  out
}

# Rows with invalid inverse-gamma parameters (NA parameters are not flagged).
.invgamma_bad <- function(shape, scale) {
  (!is.na(shape) & shape <= 0) | (!is.na(scale) & scale <= 0)
}