#' Generalized Waring Distribution
#'
#' These functions provide density, distribution function, quantile function,
#' and random number generation for the Generalized Waring Distribution.
#'
#' The Generalized Waring distribution is a 3-parameter count distribution that
#' is used to model overdispersed count data.
#'
#' @param y non-negative integer vector of count outcomes.
#' @param q non-negative integer vector of quantiles.
#' @param p numeric vector of probabilities.
#' @param n integer number of random numbers to generate.
#' @param mu numeric vector of means of the distribution.
#' @param k non-negative numeric parameter of the distribution.
#' @param rho non-negative numeric parameter of the distribution.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dgwar} computes the density (PMF) of the Generalized Waring
#' Distribution.
#'
#' \code{pgwar} computes the CDF of the Generalized Waring Distribution.
#'
#' \code{qwaring} computes the quantile function of the Generalized Waring
#' Distribution.
#'
#' \code{rwaring} generates random numbers from the Generalized Waring
#' Distribution.
#'
#' The Probability Mass Function (PMF) for the Generalized Waring (GW)
#' distribution is:
#' \deqn{f(y|a_x,k,\rho) = 
#'    \frac{\Gamma(a_x+\rho)\Gamma(k+\rho)\left(a_x\right)_y(k)_y}
#'    {y!\Gamma(\rho)\Gamma(a_x+k+\rho)(a_x+k+\rho)_y}}
#' Where \eqn{(\alpha)_r=\frac{\Gamma(\alpha+r)}{\Gamma(\alpha)}}, 
#' and \eqn{a_x, \ k, \ \rho)>0}.
#'
#' The mean value is:
#' \deqn{E[Y]=\frac{a_x K}{\rho-1}}
#'
#' Thus, we can use:
#' \deqn{a_x=\frac{\mu(\rho-1)}{k}}
#' 
#' This results in a regression model where:
#' \deqn{\mu=e^{X\beta}}
#' \deqn{\sigma^2 =
#'   \mu \left(1-\frac{1}{\alpha+\rho+1} \right) + 
#'   \mu^2\frac{(\alpha+\rho)^2}{\alpha\rho(\alpha+\rho+1)}}
#'   
#' @returns dgwar gives the density, pgwar gives the distribution function, 
#'  qgwar gives the quantile function, and rgwar generates random  deviates.
#' 
#'  The length of the result is determined by n for rgwar, and is the maximum of 
#'  the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' dgwar(0, mu=1, k=2, rho=3)
#' pgwar(c(0,1,2,3), mu=1, k=2, rho=3)
#' qgwar(0.8, mu=1, k=2, rho=3)
#' rgwar(10, mu=1, k=2, rho=3)
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
      a_i <- mu_i * k_i / (rho_i - 1)
      
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
