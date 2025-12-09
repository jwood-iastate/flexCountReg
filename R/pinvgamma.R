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
dpinvgamma <- Vectorize(function(x, mu=1, eta = 1,  log=FALSE){
  
  # Test to make sure the value of x is an integer
  # Improved integer check
  if(abs(x - round(x)) > .Machine$double.eps^0.5 || x < 0)
    warning("The value of `x` must be a non-negative whole number")
  
  if(eta <= 0) warning("The value of `eta` must be greater than 0.")
  
  # Calculation
  term1 <- 2 * (mu * (1/eta + 1))^((x + 1/eta + 2)/2)
  term2 <- besselK(2 * sqrt(mu * (1/eta + 1)), x - 1/eta - 2)
  denom <- factorial(x) * gamma(1/eta + 2)
  
  p <- (term1 * term2) / denom
  
  if (log) return(log(p))
  else return(p)
})

#' @rdname PoissonInverseGamma
#' @export
ppinvgamma <- function(q, mu = 1, eta = 1, lower.tail = TRUE,
                       log.p = FALSE) {
  
  # --- Input Validation ---
  
  if (any(eta <= 0, na.rm = TRUE)) warning("'eta' must be positive")
  if (any(mu <= 0, na.rm = TRUE)) warning("'mu' must be positive")
  
  
  # --- Vectorization Setup ---
  
  n <- max(length(q), length(mu), length(eta))
  q <- rep_len(as.integer(floor(q)), n)
  mu <- rep_len(mu, n)
  eta <- rep_len(eta, n)
  
  
  # --- Initialize Result ---
  
  cdf <- rep(NA_real_, n)
  
  
  # --- Handle Special Cases ---
  
  invalid_q <- is.na(q) | q < 0
  cdf[invalid_q] <- 0
  
  
  # --- Helper: Safe log-Bessel-K ---
  
  log_besselK_safe <- function(x, nu) {
    if (!is.finite(x) || x <= 0) return(-Inf)
    
    if (x > 700) {
      # Use scaled version: besselK(x, nu, expon.scaled=TRUE) = exp(x) *
      # K_nu(x)
      k_scaled <- tryCatch(
        besselK(x, nu, expon.scaled = TRUE),
        warning = function(w) NA_real_,
        error = function(e) NA_real_
      )
      if (is.na(k_scaled) || k_scaled <= 0) return(-Inf)
      return(log(k_scaled) - x)
    } else {
      k_val <- tryCatch(
        besselK(x, nu),
        warning = function(w) NA_real_,
        error = function(e) NA_real_
      )
      if (is.na(k_val) || k_val <= 0) return(-Inf)
      return(log(k_val))
    }
  }
  
  
  # --- Compute CDF for Valid Cases ---
  
  valid <- !invalid_q
  
  if (any(valid)) {
    for (i in which(valid)) {
      mu_i <- mu[i]
      eta_i <- eta[i]
      q_i <- q[i]
      
      # Pre-compute constants
      # PMF: p(x) = 2 * (mu*(1/eta + 1))^((x + 1/eta + 2)/2) *
      #             K_{x - 1/eta - 2}(2*sqrt(mu*(1/eta + 1))) /
      #             (x! * Gamma(1/eta + 2))
      
      inv_eta <- 1 / eta_i
      base_term <- mu_i * (inv_eta + 1)
      sqrt_base <- 2 * sqrt(base_term)
      log_denom_const <- lgamma(inv_eta + 2)
      
      # Compute log-PMF for 0:q_i
      log_pmf <- numeric(q_i + 1)
      
      for (y in 0:q_i) {
        # Bessel order
        nu_bessel <- y - inv_eta - 2
        
        # log(K_{nu}(sqrt_base))
        log_K <- log_besselK_safe(sqrt_base, nu_bessel)
        
        # Exponent: (y + 1/eta + 2) / 2
        exponent <- (y + inv_eta + 2) / 2
        
        # log(PMF) = log(2) + exponent * log(base_term) + log_K -
        # lgamma(y+1) - log_denom_const
        log_pmf[y + 1] <- log(2) + exponent * log(base_term) + log_K -
          lgamma(y + 1) - log_denom_const
      }
      
      # Log-sum-exp for CDF
      finite_mask <- is.finite(log_pmf)
      if (!any(finite_mask)) {
        cdf[i] <- NA_real_
      } else {
        max_log <- max(log_pmf[finite_mask])
        log_cdf <- max_log + log(sum(exp(log_pmf[finite_mask] -
                                           max_log)))
        cdf[i] <- exp(log_cdf)
      }
    }
  }
  
  
  # --- Apply Tail and Log Options ---
  
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  return(cdf)
}

#' @rdname PoissonInverseGamma
#' @export
qpinvgamma <- Vectorize(function(p, mu=1, eta = 1) {
  y <- 0
  p_value <- ppinvgamma(y, mu, eta)
  while(p_value < p){
    y <- y + 1
    p_value <- ppinvgamma(y, mu, eta)
  }
  return(y)
})

#' @rdname PoissonInverseGamma
#' @export
rpinvgamma <- function(n, mu=1, eta = 1) {
  u <- runif(n)
  y <- sapply(u, function(p) qpinvgamma(p, mu, eta))
  return(y)
}
