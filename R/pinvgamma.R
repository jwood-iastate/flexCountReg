#' Poisson-Inverse-Gamma Distribution
#'
#' These functions provide the density function, distribution function,
#' quantile function, and random number generation for the
#' Poisson-Inverse-Gamma (PInvGamma) Distribution
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mu numeric value or vector of mean values for the distribution (the
#' values have to be greater than 0).
#' @param eta single value or vector of values for the scale parameter of the
#' distribution (the values have to be greater than 0).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#' otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dpinvgamma} computes the density (PDF) of the Poisson-Inverse-Gamma
#' Distribution.
#'
#' \code{ppinvgamma} computes the CDF of the Poisson-Inverse-Gama Distribution.
#'
#' \code{qpinvgamma} computes the quantile function of the
#' Poisson-Inverse-Gamma Distribution.
#'
#' \code{rpinvgamma} generates random numbers from the Poisson-Inverse-Gamma
#' Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Poisson-Inverse-Gamma
#' distribution is:
#' \deqn{
#' f(x|\eta,\mu)=\frac{2\left(\mu\left(\frac{1}{\eta}+1\right)\right)^{
#' \frac{x+\frac{1}{eta}+2}{2}}}{x!\Gamma\left(\frac{1}{\eta}+2\right)}
#' K_{x-\frac{1}{\eta}-2}\left(2\sqrt{\mu\left(\frac{1}{\eta}+1\right)}\right)
#' }
#'
#' Where \eqn{\eta} is a shape parameter with the restriction that
#' \eqn{\eta>0}, \eqn{\mu>0} is the mean value,  \eqn{y} is a non-negative
#' integer, and \eqn{K_i(z)} is the modified Bessel function of the second
#' kind. This formulation uses the mean directly.
#'
#' The variance of the distribution is:
#' \deqn{\sigma^2=\mu+\eta\mu^2}
#' 
#' @returns dpinvgamma gives the density, ppinvgamma gives the distribution 
#'  function, qpinvgamma gives the quantile function, and rcom generates random 
#'  deviates.
#' 
#' The length of the result is determined by n for rpinvgamma, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' dpinvgamma(1, mu=0.75, eta=1)
#' ppinvgamma(c(0,1,2,3,5,7,9,10), mu=0.75, eta=3)
#' qpinvgamma(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, eta=0.5)
#' rpinvgamma(30, mu=0.75, eta=1.5)
#'
#' @importFrom stats runif
#' @export
#' @name PoissonInverseGamma
#' @rdname PoissonInverseGamma
#' @export
dpinvgamma <- function(x, mu = 1, eta = 1, log = FALSE) {
  r <- .fcr_recycle(x, mu, eta)
  x <- r[[1]]; mu <- r[[2]]; eta <- r[[3]]
  
  is_count <- .fcr_is_count(x)
  if (any(!is_count & !is.na(x))) {
    warning("The value of `x` must be a non-negative whole number.")
  }
  bad_par <- !is.na(mu) & !is.na(eta) & (mu <= 0 | eta <= 0)
  if (any(bad_par)) {
    warning("The values of `mu` and `eta` must be greater than 0.")
  }
  
  xx <- ifelse(is_count, round(x), 0)
  lp <- suppressWarnings(.ldpinvgamma(xx, mu, eta))
  lp[!is_count] <- -Inf
  lp[bad_par] <- NaN
  lp[is.na(x) | is.na(mu) | is.na(eta)] <- NA
  if (log) lp else exp(lp)
}

# Log PMF of the Poisson-inverse-gamma:
#   f(x) = 2 B^((x + c)/2) K_{x - c}(2 sqrt(B)) / (x! Gamma(c)),
# with c = 1/eta + 2 and B = mu (1/eta + 1). Uses lgamma() and a log Bessel
# function that stays finite for large orders.
.ldpinvgamma <- function(x, mu, eta) {
  inv_eta <- 1 / eta
  cc <- inv_eta + 2
  B <- mu * (inv_eta + 1)
  log(2) + (x + cc) / 2 * log(B) +
    .fcr_log_besselK(2 * sqrt(B), x - cc) -
    lgamma(x + 1) - lgamma(cc)
}

#' @rdname PoissonInverseGamma
#' @export
ppinvgamma <- function(q, mu = 1, eta = 1, lower.tail = TRUE,
                       log.p = FALSE) {
  .fcr_p_contract(.ppinvgamma_core, q, list(mu = mu, eta = eta),
                  lower.tail, log.p,
                  logdens = function(x, ...) dpinvgamma(x, ..., log = TRUE))
}

.ppinvgamma_core <- function(q, mu = 1, eta = 1, lower.tail = TRUE,
                       log.p = FALSE) {
  if (any(eta <= 0, na.rm = TRUE)) warning("'eta' must be positive")
  if (any(mu <= 0, na.rm = TRUE)) warning("'mu' must be positive")
  
  r <- .fcr_recycle(q, mu, eta)
  q <- floor(r[[1]]); mu <- r[[2]]; eta <- r[[3]]
  cdf <- rep(NA_real_, length(q))
  cdf[!is.na(q) & q < 0] <- 0
  bad_par <- !(is.finite(mu) & mu > 0 & is.finite(eta) & eta > 0)
  cdf[bad_par & !is.na(q) & q >= 0] <- NaN
  for (i in which(!is.na(q) & q >= 0 & !bad_par)) {
    lp <- suppressWarnings(.ldpinvgamma(0:q[i], mu[i], eta[i]))
    m <- max(lp)
    cdf[i] <- if (is.finite(m)) min(1, exp(m + log(sum(exp(lp - m))))) else NA
  }
  
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  cdf
}

#' @rdname PoissonInverseGamma
#' @export
qpinvgamma <- function(p, mu = 1, eta = 1) {
  .fcr_qcount(p, function(y, mu, eta) {
    cumsum(exp(.ldpinvgamma(y, mu, eta)))
  }, mu, eta)
}

#' @rdname PoissonInverseGamma
#' @export
rpinvgamma <- function(n, mu = 1, eta = 1) {
  
  if(length(n) != 1L || is.na(n) || n < 0 || n != floor(n)) {
    stop("'n' must be a non-negative integer.")
  }
  
  n <- as.integer(n)
  
  if(n == 0L) {
    return(numeric(0))
  }
  
  if(any(mu <= 0, na.rm = TRUE)) {
    stop("'mu' must be positive.")
  }
  
  if(any(eta <= 0, na.rm = TRUE)) {
    stop("'eta' must be positive.")
  }
  
  mu <- rep_len(mu, n)
  eta <- rep_len(eta, n)
  
  alpha <- 1 / eta + 2
  beta  <- mu * (1 / eta + 1)
  
  lambda <- rinvgamma(
    n = n,
    shape = alpha,
    scale = beta
  )
  
  rpois(n, lambda)
}