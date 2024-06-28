#' Inverse Gamma Distribution
#'
#' These functions provide the density function, distribution function, quantile function, and random number generation for the Inverse-Gamma (IG) Distribution
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param shape numeric value or vector of shape values for the distribution (the values have to be greater than 0).
#' @param scale single value or vector of values for the scale parameter of the distribution (the values have to be greater than 0).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dinvgamma} computes the density (PDF) of the Inverse-Gamma Distribution.
#'
#' \code{pinvgamma} computes the CDF of the Inverse-Gamma Distribution.
#'
#' \code{qinvgamma} computes the quantile function of the Inverse-Gamma Distribution.
#'
#' \code{rinvgamma} generates random numbers from the Inverse-Gamma Distribution.
#'
#' The Probability Density Function (PDF) for the Inverse-Gamma distribution:
#' \deqn{f(x|\alpha,\beta)=\frac{\beta^\alpha}{\Gamma(\alpha)}\left(\frac{1}{x}\right)^{\alpha+1}e^{-\frac{\beta}{x}}}
#' 
#' Where \eqn{\alpha} is the shape parameter and \eqn{\beta} is a scale parameter with the restrictions that \eqn{\alpha>0} and \eqn{\eta>0}, and \eqn{x>0}.
#' 
#' The CDF of the Inverse-Gamma distribution is:
#' \deqn{F(x|\alpha,\beta)=\frac{\alpha. \Gamma\left(\frac{\beta}{x}\right)}{\Gamma(\alpha)}=Q\left(\alpha, \frac{\beta}{x} \right)}
#' 
#' Where the numerator is the incomplete gamma function and \eqn{Q(\cdot)} is the regularized gamma function.
#'
#' The mean of the distribution is (provided \eqn{\alpha>1}):
#' \deqn{\mu=\frac{\beta}{\alpha-1}}
#' 
#' The variance of the distribution is (for \eqn{\alpha>2}):
#' \deqn{\sigma^2=\frac{\beta^2}{(\alpha-1)^2(\alpha-2)}}
#' 
#' @examples
#' dinvgamma(1, shape=3, scale=2)
#' pinvgamma(c(0.1, 0.5, 1, 3, 5, 10, 30), shape=3, scale=2)
#' qinvgamma(c(0.1,0.3,0.5,0.9,0.95), shape=3, scale=2)
#' rinvgamma(30, shape=3, scale=2)
#'
#' @import stats
#' @export
#' @name invgamma

#' @rdname invgamma
#' @export
dinvgamma <- Vectorize(function(x, shape = 2.5, scale = 1, log = FALSE) {
  if (x <= 0) { # Ensure the input x is strictly positive
    warning("x must be greater than 0")
    return(NA)  # Return NA for invalid input values
  }
  
  # Calculate the log-density using the gamma density function
  # Convert scale to rate for the gamma function
  rate <- 1 / scale
  # Logarithm of the gamma density evaluated at 1/x
  log.p <- stats::dgamma(1/x, shape, rate = rate, log = TRUE) - 2 * log(x)
  
  # Return log density or convert log density to density
  if (log) {
    return(log.p)
  } else {
    return(exp(log.p))
  }
})


#' @rdname invgamma
#' @export
pinvgamma <- Vectorize(function(q, shape=2.5, scale = 1, lower.tail=TRUE, log.p=FALSE){
  p <- stats::pgamma(1/q, shape, rate=1/scale, lower.tail = !lower.tail, log.p = log.p)
  return(p)
})

#' @rdname invgamma
#' @export
qinvgamma <- Vectorize(function(p,shape=2.5, scale = 1, lower.tail = TRUE,  log.p = FALSE) {
  q <- stats::qgamma(1-p, shape, rate=1/scale, lower.tail = lower.tail, log.p = log.p)^(-1)
  return(q)
})


#' @rdname invgamma
#' @export
rinvgamma <- function(n, shape=2.5, scale = 1) {
  u <- runif(n)
  y <- sapply(u, function(p) qinvgamma(p, shape, scale))
  return(y)
}
