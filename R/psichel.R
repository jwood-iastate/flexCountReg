#' Sichel Distribution
#'
#' @description
#' Density, distribution function, quantile function, and random generation
#' for the Sichel distribution.
#'
#' @param x numeric value or vector of non-negative integer values.
#' @param q quantile or vector of quantiles.
#' @param p probability or vector of probabilities.
#' @param n number of random values to generate.
#' @param mu numeric; mean of the distribution (mu > 0).
#' @param sigma numeric; scale parameter (sigma > 0).
#' @param gamma numeric; shape parameter (can be any real number).
#' @param log,log.p logical; if TRUE, probabilities are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities are P[X <= x].
#'
#' @details
#' The Sichel distribution is a three-parameter discrete distribution that
#' generalizes the Poisson-inverse Gaussian distribution. It is useful for
#' modeling overdispersed count data.
#'
#' The PMF is:
#' \deqn{f(y|\mu, \sigma, \gamma) = 
#' \frac{(\mu/c)^y K_{y+\gamma}(\alpha)}{K_\gamma(1/\sigma) y! 
#' (\alpha\sigma)^{y+\gamma}}}
#'
#' @references
#' Rigby, R. A., Stasinopoulos, D. M., & Akantziliotou, C. (2008).
#' A framework for modelling overdispersed count data, including the
#' Poisson-shifted generalized inverse Gaussian distribution.
#' Computational Statistics & Data Analysis, 53(2), 381-393.
#' 
#' @details dsichel gives the density, psichel gives the distribution 
#'  function, qsichel gives the quantile function, and rsichel 
#'  generates random  deviates.
#' 
#'  The length of the result is determined by n for rsichel, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' # Basic usage
#' dsichel(0:10, mu = 5, sigma = 1, gamma = -0.5)
#'
#' # Log-probabilities for numerical stability
#' dsichel(0:10, mu = 5, sigma = 1, gamma = -0.5, log = TRUE)
#'
#' # CDF
#' psichel(5, mu = 5, sigma = 1, gamma = -0.5)
#'
#' @name SichelDistribution
#' @export
dsichel <- function(x, mu = 1, sigma = 1, gamma = 1, log = FALSE) {
  
  
  # === Input Validation ===
  
  
  # Check parameter constraints
  if (any(mu <= 0, na.rm = TRUE)) warning("'mu' must be positive")
  if (any(sigma <= 0, na.rm = TRUE)) warning("'sigma' must be positive")
  # gamma can be any real number, no constraint needed
  
  
  # === Vectorization Setup ===
  
  
  # Determine output length and recycle inputs
  n <- max(length(x), length(mu), length(sigma), length(gamma))
  x <- rep_len(x, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  gamma <- rep_len(gamma, n)
  
  
  # === Initialize Output ===
  
  
  log_p <- rep(-Inf, n)
  
  
  # === Identify Valid Cases ===
  
  
  # Valid: non-NA, non-negative integers
  valid <- !is.na(x) & !is.na(mu) & !is.na(sigma) & !is.na(gamma) &
    is.finite(x) & x >= 0 & x == floor(x)
  
  if (!any(valid)) {
    if (log) return(log_p) else return(exp(log_p))
  }
  
  
  # === Core Computation (for valid cases only) ===
  
  
  x_v <- x[valid]
  mu_v <- mu[valid]
  sigma_v <- sigma[valid]
  gamma_v <- gamma[valid]
  
  # Pre-compute common terms
  inv_sigma <- 1 / sigma_v
  
  # Compute Bessel functions with error handling
  # K_gamma(1/sigma) - appears in denominator
  log_K_g <- log_besselK_safe(inv_sigma, gamma_v)
  
  # K_{gamma+1}(1/sigma) - needed for c
  log_K_g1 <- log_besselK_safe(inv_sigma, gamma_v + 1)
  
  # Check for invalid Bessel computations
  bessel_valid <- is.finite(log_K_g) & is.finite(log_K_g1)
  
  if (!any(bessel_valid)) {
    warning("Bessel function computation failed for all inputs. ",
            "Consider different parameter values.")
    if (log) return(log_p) else return(exp(log_p))
  }
  
  # Compute c = K_{gamma+1}(1/sigma) / K_gamma(1/sigma)
  # In log space: log(c) = log_K_g1 - log_K_g
  log_c <- log_K_g1 - log_K_g
  c_val <- exp(log_c)
  
  # Compute alpha = sqrt(sigma^{-2} + 2*mu/(c*sigma))
  # Need to handle this carefully to avoid overflow
  alpha_sq <- inv_sigma^2 + 2 * mu_v / (c_val * sigma_v)
  
  # Check for valid alpha
  alpha_valid <- bessel_valid & is.finite(alpha_sq) & alpha_sq > 0
  
  if (!any(alpha_valid)) {
    warning("Alpha computation failed. Consider different parameter values.")
    if (log) return(log_p) else return(exp(log_p))
  }
  
  alpha <- sqrt(alpha_sq)
  
  # Compute K_{x+gamma}(alpha)
  log_K_xg_a <- log_besselK_safe(alpha, x_v + gamma_v)
  
  # Final validity check
  final_valid <- alpha_valid & is.finite(log_K_xg_a)
  
  # === Compute Log-Probability ===
  
  # log f(y) = y*log(mu/c) + log(K_{y+gamma}(alpha)) 
  #            - log(K_gamma(1/sigma)) - log(y!) - (y+gamma)*log(alpha*sigma)
  
  log_p_calc <- x_v * (log(mu_v) - log_c) +
    log_K_xg_a -
    log_K_g -
    lgamma(x_v + 1) -  # log(y!) = lgamma(y+1)
    (x_v + gamma_v) * log(alpha * sigma_v)
  
  # Only assign valid results
  log_p_calc[!final_valid] <- -Inf
  
  # Handle any remaining numerical issues
  log_p_calc[!is.finite(log_p_calc)] <- -Inf
  
  # Assign back to full result vector
  log_p[valid] <- log_p_calc
  
  
  # === Return Result ===
  
  
  if (log) {
    return(log_p)
  } else {
    return(exp(log_p))
  }
}


#' Safe computation of log(besselK)
#' 
#' Computes log of the modified Bessel function of the second kind,
#' handling edge cases that cause besselK to return Inf, 0, or NaN.
#' 
#' @param x numeric; argument of the Bessel function (x > 0)
#' @param nu numeric; order of the Bessel function
#' @return log(besselK(x, nu)), or -Inf/Inf for edge cases
#' @noRd
log_besselK_safe <- function(x, nu) {
  n <- length(x)
  result <- rep(-Inf, n)
  
  # besselK is only defined for x > 0
  valid <- is.finite(x) & x > 0 & is.finite(nu)
  
  if (!any(valid)) {
    return(result)
  }
  
  x_v <- x[valid]
  nu_v <- nu[valid]
  
  # Compute besselK
  # Use expon.scaled = TRUE for large x to avoid underflow
  # besselK(x, nu, expon.scaled = TRUE) returns exp(x) * K_nu(x)
  
  # Determine which values need scaling
  needs_scaling <- x_v > 700  # exp(709) overflows, so be conservative
  
  result_v <- numeric(length(x_v))
  
  # For moderate x, use direct computation
  if (any(!needs_scaling)) {
    K_direct <- besselK(x_v[!needs_scaling], nu_v[!needs_scaling])
    
    # Handle cases where besselK returns 0 (underflow) or Inf (overflow)
    K_direct[K_direct == 0] <- .Machine$double.xmin
    K_direct[!is.finite(K_direct)] <- NA
    
    result_v[!needs_scaling] <- log(K_direct)
  }
  
  # For large x, use scaled version
  if (any(needs_scaling)) {
    K_scaled <- besselK(
      x_v[needs_scaling], nu_v[needs_scaling], expon.scaled = TRUE)
    
    # log(K_nu(x)) = log(exp(-x) * K_scaled) = -x + log(K_scaled)
    K_scaled[K_scaled == 0] <- .Machine$double.xmin
    K_scaled[!is.finite(K_scaled)] <- NA
    
    result_v[needs_scaling] <- -x_v[needs_scaling] + log(K_scaled)
  }
  
  result[valid] <- result_v
  
  return(result)
}

#' @rdname SichelDistribution
#' @export
psichel <- function(q, mu = 1, sigma = 1, gamma = 1,
                    lower.tail = TRUE, log.p = FALSE) {
  
  # --- Input Validation ---
  if (any(sigma <= 0, na.rm = TRUE)) warning("'sigma' must be positive")
  if (any(mu <= 0, na.rm = TRUE)) warning("'mu' must be positive")
  
  # --- Vectorization Setup ---
  n <- max(length(q), length(mu), length(sigma), length(gamma))
  q <- rep_len(as.integer(floor(q)), n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  gamma <- rep_len(gamma, n)
  
  # --- Initialize Result ---
  cdf <- rep(NA_real_, n)
  
  # --- Handle Special Cases ---
  invalid_q <- is.na(q) | q < 0
  cdf[invalid_q] <- 0
  
  # --- Compute CDF for Valid Cases ---
  valid <- !invalid_q
  
  if (any(valid)) {
    for (i in which(valid)) {
      mu_i <- mu[i]
      sigma_i <- sigma[i]
      gamma_i <- gamma[i]
      q_i <- q[i]
      
      # Pre-compute constants
      inv_sigma <- 1 / sigma_i
      
      # Note: Using the global log_besselK_safe function defined previously
      log_K_gamma_inv_sigma <- log_besselK_safe(inv_sigma, gamma_i)
      log_K_gamma1_inv_sigma <- log_besselK_safe(inv_sigma, gamma_i + 1)
      
      # c = K_{gamma+1}(1/sigma) / K_gamma(1/sigma)
      log_c <- log_K_gamma1_inv_sigma - log_K_gamma_inv_sigma
      c_val <- exp(log_c)
      
      # Check validity
      if (!is.finite(c_val) || c_val <= 0) {
        cdf[i] <- NA_real_
        next
      }
      
      # alpha^2 = sigma^{-2} + 2*mu/(c*sigma)
      alpha_sq <- inv_sigma^2 + 2 * mu_i / (c_val * sigma_i)
      if (alpha_sq <= 0) {
        cdf[i] <- NA_real_
        next
      }
      alpha_val <- sqrt(alpha_sq)
      
      # Compute log-PMF for 0:q_i
      # If q_i is 0, this creates a vector of length 1, which is safe.
      log_pmf <- numeric(q_i + 1)
      
      for (y in 0:q_i) {
        log_mu_over_c_y <- y * (log(mu_i) - log_c)
        log_K_y_gamma_alpha <- log_besselK_safe(alpha_val, y + gamma_i)
        log_factorial_y <- lgamma(y + 1)
        log_alpha_sigma_pow <- (y + gamma_i) * log(alpha_val * sigma_i)
        
        log_pmf[y + 1] <- log_mu_over_c_y + log_K_y_gamma_alpha - 
          log_K_gamma_inv_sigma - log_factorial_y - log_alpha_sigma_pow
      }
      
      # Log-sum-exp for cumulative probability
      max_log <- max(log_pmf[is.finite(log_pmf)])
      if (!is.finite(max_log)) {
        cdf[i] <- NA_real_
        next
      }
      
      # sum() works safely on length 1 vectors
      log_cdf <- max_log + log(sum(exp(log_pmf - max_log)))
      cdf[i] <- exp(log_cdf)
    }
  }
  
  # --- Apply Tail and Log Options ---
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  return(cdf)
}

#' @rdname SichelDistribution
#' @export
qsichel <- function(p, mu = 1, sigma = 1, gamma = 1,
                    lower.tail = TRUE, log.p = FALSE) {
  
  # Input validation
  if (any(mu <= 0, na.rm = TRUE)) warning("'mu' must be positive")
  if (any(sigma <= 0, na.rm = TRUE)) warning("'sigma' must be positive")
  
  # Handle log.p
  if (log.p) p <- exp(p)
  
  # Handle lower.tail
  if (!lower.tail) p <- 1 - p
  
  # Validate probabilities
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    warning("'p' must be in [0, 1]")
  }
  
  # Vectorize
  n <- max(length(p), length(mu), length(sigma), length(gamma))
  p <- rep_len(p, n)
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  gamma <- rep_len(gamma, n)
  
  # Compute quantiles
  q <- mapply(function(pi, mui, sigmai, gammai) {
    if (is.na(pi)) return(NA_integer_)
    if (pi == 0) return(0L)
    if (pi == 1) return(Inf)
    
    y <- 0L
    cum_p <- dsichel(0, mu = mui, sigma = sigmai, gamma = gammai)
    
    # Safety check: if parameters are invalid, dsichel might return NA
    if (is.na(cum_p)) return(NA_integer_)
    
    max_iter <- 10000 
    
    # Added !is.na check to prevent "missing value where TRUE/FALSE needed"
    # if cum_p becomes NA during iteration
    while (!is.na(cum_p) && cum_p < pi && y < max_iter) {
      y <- y + 1L
      prob_y <- dsichel(y, mu = mui, sigma = sigmai, gamma = gammai)
      if(is.na(prob_y)) break 
      cum_p <- cum_p + prob_y
    }
    
    if (y >= max_iter) {
      warning("Maximum iterations reached in qsichel")
    }
    
    return(y)
    
  }, p, mu, sigma, gamma, SIMPLIFY = TRUE)
  
  return(as.integer(q))
}

#' @rdname SichelDistribution
#' @export
rsichel <- function(n, mu = 1, sigma = 1, gamma = 1) {
  
  if (length(n) != 1 || n < 0 || n != floor(n)) 
    warning("'n' must be a non-negative integer")
  
  if (n == 0) return(integer(0))
  
  # Input validation
  if (any(mu <= 0)) warning("'mu' must be positive")
  if (any(sigma <= 0)) warning("'sigma' must be positive")
  
  # Generate uniform random numbers and use inverse CDF
  u <- runif(n)
  
  # Recycle parameters if needed
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  gamma <- rep_len(gamma, n)
  
  # Use qsichel for each uniform variate
  qsichel(u, mu = mu, sigma = sigma, gamma = gamma)
}

