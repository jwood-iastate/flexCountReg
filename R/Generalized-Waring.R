#' Generalized Waring Distribution
#'
#' These functions provide the probability mass function, distribution
#' function, quantile function, and random number generation for the
#' Generalized Waring distribution.
#'
#' The Generalized Waring distribution is a three-parameter count distribution
#' used to model overdispersed count data.
#'
#' @param y non-negative integer vector of count outcomes.
#' @param q non-negative integer vector of quantiles.
#' @param p numeric vector of probabilities.
#' @param n integer number of random numbers to generate.
#' @param mu positive numeric vector of distribution means.
#' @param k positive numeric vector of shape parameters.
#' @param rho numeric vector of shape parameters greater than 1. A finite
#'   variance requires \eqn{\rho > 2}.
#' @param log logical; if TRUE, log probabilities are returned.
#' @param log.p logical; if TRUE, probabilities are returned as log
#'   probabilities.
#' @param lower.tail logical; if TRUE, probabilities are
#'   \eqn{P[X \leq q]}; otherwise, \eqn{P[X > q]}.
#'
#' @details
#' \code{dgwar} computes the probability mass function (PMF) of the
#' Generalized Waring distribution.
#'
#' \code{pgwar} computes the CDF of the Generalized Waring distribution.
#'
#' \code{qgwar} computes the quantile function of the Generalized Waring
#' distribution.
#'
#' \code{rgwar} generates random numbers from the Generalized Waring
#' distribution.
#'
#' The probability mass function for the Generalized Waring distribution is:
#' \deqn{
#' f(y \mid a_x, k, \rho) =
#' \frac{
#'   \Gamma(a_x + \rho)\Gamma(k + \rho)
#'   (a_x)_y(k)_y
#' }{
#'   y!\Gamma(\rho)\Gamma(a_x + k + \rho)
#'   (a_x + k + \rho)_y
#' }
#' }
#' where
#' \eqn{(\alpha)_r = \frac{\Gamma(\alpha+r)}{\Gamma(\alpha)}},
#' and \eqn{a_x > 0}, \eqn{k > 0}, and \eqn{\rho > 0}.
#'
#' When \eqn{\rho > 1}, the mean is:
#' \deqn{
#' E[Y] = \frac{a_x k}{\rho - 1}.
#' }
#'
#' Therefore, the distribution can be parameterized in terms of its mean using:
#' \deqn{
#' a_x = \frac{\mu(\rho - 1)}{k}.
#' }
#'
#' For a regression model:
#' \deqn{
#' \mu = \exp(X\beta).
#' }
#'
#' When \eqn{\rho > 2}, the variance under this mean parameterization is:
#' \deqn{
#' \mathrm{Var}(Y) =
#' \frac{\mu(k+\mu)(k+\rho-1)}{k(\rho-2)}.
#' }
#'
#' @returns
#' \code{dgwar} gives the probability mass function, \code{pgwar} gives the
#' distribution function, \code{qgwar} gives the quantile function, and
#' \code{rgwar} generates random deviates.
#'
#' The length of the result is determined by \code{n} for \code{rgwar} and is
#' the maximum of the lengths of the numerical arguments for the other
#' functions.
#'
#' @examples
#' dgwar(0, mu = 1, k = 2, rho = 3)
#' pgwar(c(0, 1, 2, 3), mu = 1, k = 2, rho = 3)
#' qgwar(0.8, mu = 1, k = 2, rho = 3)
#' rgwar(10, mu = 1, k = 2, rho = 3)
#'
#' @importFrom stats runif
#' @export
#' @name Generalized-Waring
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg

#' @rdname Generalized-Waring
#' @export
dgwar <- Vectorize(function(y, mu, k, rho, log = FALSE) {

  pmf <- genWaring_cpp(y, mu, k, rho)
  
  if (log) pmf <- log(pmf)
  return(pmf)
})


#' @rdname Generalized-Waring
#' @export
pgwar <- function(q, mu, k, rho, lower.tail = TRUE, log.p = FALSE) {
  
  # --- 1. Vectorization Setup ---
  # Ensure all inputs are the same length before checking values
  n <- max(length(q), length(mu), length(k), length(rho))
  q   <- rep_len(q, n)
  mu  <- rep_len(mu, n)
  k   <- rep_len(k, n)
  rho <- rep_len(rho, n)
  
  # --- 2. Initialize Result ---
  # Default to NaN. If parameters are invalid, this remains NaN.
  cdf <- rep(NaN, n)
  
  # --- 3. Identify Cases ---
  
  # Valid parameters: mu > 0, k > 0, rho > 1 (rho must be >1 for the mean
  # formula used)
  valid_params <- 
    (mu > 0) & (k > 0) & (rho > 1) & !is.na(mu) & !is.na(k) & !is.na(rho)
  
  # Case A: Parameters are valid, but q < 0. CDF is 0.
  # (We treat NA in q as resulting in NA, so we check !is.na(q))
  is_neg_q <- valid_params & !is.na(q) & (q < 0)
  cdf[is_neg_q] <- 0
  
  # Case B: Parameters are valid and q >= 0. Calculate CDF.
  calc_idx <- which(valid_params & !is.na(q) & (q >= 0))
  
  # --- 4. Compute CDF for Valid Cases ---
  if (length(calc_idx) > 0) {
    for (i in calc_idx) {
      mu_i  <- mu[i]
      k_i   <- k[i]
      rho_i <- rho[i]
      q_i   <- floor(q[i]) # Integers only for summation
      
      # Calculate parameter 'a' derived from the mean formula
      a_i <- mu_i * (rho_i - 1) / k_i
      
      # Double check a_i valid just in case
      if (a_i <= 0) {
        cdf[i] <- NaN
        next
      }
      
      # Compute log-PMF for y in 0..q_i
      # We calculate vectorially for 0:q_i to avoid inner loops overhead
      y_seq <- 0:q_i
      
      # Prepare terms for lgamma
      # Numerator terms: Gamma(a+y) + Gamma(k+y) + Gamma(k+rho) + Gamma(a+rho)
      log_num <- lgamma(a_i + y_seq) + lgamma(k_i + y_seq) + 
        lgamma(k_i + rho_i) + lgamma(a_i + rho_i)
      
      # Denominator terms: 
      # Gamma(a) + Gamma(k) + Gamma(rho) + Gamma(a+k+rho+y) + y!
      # Note: 
      # (a+k+rho)_y in denominator combines with 
      # Gamma(a+k+rho) to form Gamma(a+k+rho+y)
      log_den <- lgamma(a_i) + lgamma(k_i) + lgamma(rho_i) + 
        lgamma(a_i + k_i + rho_i + y_seq) + lgamma(y_seq + 1)
      
      log_pmf <- log_num - log_den
      
      # Log-sum-exp for stability
      max_log <- max(log_pmf)
      
      # Check for numerical issues
      if (!is.finite(max_log)) {
        cdf[i] <- NaN 
      } else {
        log_cdf <- max_log + log(sum(exp(log_pmf - max_log)))
        cdf[i]  <- exp(log_cdf)
      }
    }
  }
  
  # --- 5. Apply Tail and Log Options ---
  
  # Clamp result to [0,1] to handle tiny floating point errors
  cdf[valid_params & (q >= 0) & !is.nan(cdf)] <- 
    pmin(pmax(cdf[valid_params & (q >= 0) & !is.nan(cdf)], 0), 1)
  
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  return(cdf)
}


#' @rdname Generalized-Waring
#' @export
qgwar <- Vectorize(function(p, mu, k, rho) {
  if (any(p < 0 | p > 1)) {
    warning("All p values must be in the interval [0, 1].")
  }
  
  y <- 0
  p_value <- pgwar(y, mu, k, rho)
  while (p_value < p) {
    y <- y + 1
    p_value <- pgwar(y, mu, k, rho)
  }
  return(y)
})


#' @rdname Generalized-Waring
#' @export
rgwar <- function(n, mu, k, rho) {
  p <- runif(n)
  random_counts <- qgwar(p, mu, k, rho)
  return(random_counts)
}
