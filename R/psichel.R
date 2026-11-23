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
#' @param lower.tail logical; if TRUE, probabilities are \eqn{P[X <= x]}.
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
#' @returns dsichel gives the density, psichel gives the distribution 
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
  bad_par <- !is.na(mu) & !is.na(sigma) & (mu <= 0 | sigma <= 0)
  log_p[bad_par] <- NaN
  log_p[is.na(x) | is.na(mu) | is.na(sigma) | is.na(gamma)] <- NA_real_
  valid <- !is.na(x) & !is.na(mu) & !is.na(sigma) & !is.na(gamma) &
    is.finite(x) & x >= 0 & x == floor(x) & !bad_par
  
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
  # Order recurrence on the log scale (.fcr_log_besselK): besselK() itself
  # underflows to 0 for large arguments and overflows for large orders, and
  # clamping those values gave -Inf log-probabilities where the PMF is finite.
  out <- .fcr_log_besselK(x, nu)
  out[is.na(out)] <- -Inf
  out
}

#' @rdname SichelDistribution
#' @export
psichel <- function(q, mu = 1, sigma = 1, gamma = 1,
                    lower.tail = TRUE, log.p = FALSE) {
  .fcr_p_contract(.psichel_core, q, list(mu = mu, sigma = sigma, gamma = gamma),
                  lower.tail, log.p,
                  logdens = function(x, ...) dsichel(x, ..., log = TRUE))
}

.psichel_core <- function(q, mu = 1, sigma = 1, gamma = 1,
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
  invalid_q <- !is.na(q) & q < 0
  cdf[invalid_q] <- 0
  bad_par <- !is.na(mu) & !is.na(sigma) & (mu <= 0 | sigma <= 0)
  cdf[bad_par] <- NaN
  
  # --- Compute CDF for Valid Cases ---
  valid <- !is.na(q) & !invalid_q & !bad_par & !is.na(mu) & !is.na(sigma) &
    !is.na(gamma)
  
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
      
      # Log-PMF for 0:q_i in one sweep of the Bessel order recurrence
      log_pmf <- .fcr_sichel_logpmf_seq(q_i, mu_i, sigma_i, gamma_i)
      
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
  
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  if (any(p < 0 | p > 1, na.rm = TRUE)) {
    warning("'p' must be in [0, 1]")
  }
  
  r <- .fcr_recycle(p, mu, sigma, gamma)
  p <- r[[1]]; mu <- r[[2]]; sigma <- r[[3]]; gamma <- r[[4]]
  
  # The support is searched in blocks that double in size (.fcr_qsearch), and
  # each block's PMF comes from one sweep of the Bessel order recurrence, so
  # the cost is linear in the quantile. Stepping one count at a time with an
  # iteration cap stopped short for heavy-tailed parameters.
  vapply(seq_along(p), function(i) {
    if (is.na(p[i]) || is.na(gamma[i])) return(NA_real_)
    if (is.na(mu[i]) || is.na(sigma[i]) || mu[i] <= 0 || sigma[i] <= 0) {
      return(NaN)
    }
    .fcr_qsearch(p[i], function(y) {
      cumsum(exp(.fcr_sichel_logpmf_seq(max(y), mu[i], sigma[i], gamma[i])))
    })
  }, numeric(1))
}

# Log-PMF of the Sichel distribution at 0, 1, ..., K for scalar parameters.
# log K_{y + gamma}(alpha) comes from the upward order recurrence
#   K_{v+1}(z) = K_{v-1}(z) + (2 v / z) K_v(z),
# carried out on the log scale. It is stable for v > 0, so orders up to 1
# (and all orders below 0, where the 2 v / z term is negative) are computed
# directly. The cost is O(K) instead of O(K^2) for K separate evaluations.
.fcr_sichel_logpmf_seq <- function(K, mu, sigma, gamma) {
  y <- 0:K
  lKg <- .fcr_log_besselK(1 / sigma, gamma)
  log_c <- .fcr_log_besselK(1 / sigma, gamma + 1) - lKg
  alpha <- sqrt(1 / sigma^2 + 2 * mu / (exp(log_c) * sigma))
  if (!is.finite(lKg) || !is.finite(log_c) || !is.finite(alpha)) {
    return(rep(NA_real_, K + 1))
  }
  v <- y + gamma
  lK <- numeric(K + 1)
  n_direct <- min(K + 1, max(2, which(v > 1)[1], na.rm = TRUE))
  lK[seq_len(n_direct)] <- .fcr_log_besselK(alpha, v[seq_len(n_direct)])
  if (n_direct < K + 1) {
    for (k in (n_direct + 1):(K + 1)) {
      # log K_{v_k} from log K_{v_k - 2} and log K_{v_k - 1}
      lK[k] <- lK[k - 1] + log(exp(lK[k - 2] - lK[k - 1]) +
                                 2 * v[k - 1] / alpha)
    }
  }
  y * (log(mu) - log_c) + lK - lKg - lgamma(y + 1) -
    (y + gamma) * log(alpha * sigma)
}

#' @rdname SichelDistribution
#' @export
rsichel <- function(n, mu = 1, sigma = 1, gamma = 1) {
  
  if (length(n) > 1) n <- length(n)
  if (n < 0 || n != floor(n)) stop("'n' must be a non-negative integer")
  
  if (n == 0) return(integer(0))
  
  # Parameter checks (and their warnings) are done by qsichel()
  # Generate uniform random numbers and use inverse CDF
  u <- runif(n)
  
  # Recycle parameters if needed
  mu <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  gamma <- rep_len(gamma, n)
  
  # Use qsichel for each uniform variate
  qsichel(u, mu = mu, sigma = sigma, gamma = gamma)
}

