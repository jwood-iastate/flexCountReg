#' @title One-Parameter Lindley Distribution
#' 
#' @description Distribution function for the one-parameter Lindley distribution with parameter theta.
#'
#' @param x a single value or vector of positive values.
#' @param p a single value or vector of probabilities.
#' @param q a single value or vector of quantiles.
#' @param n number of random values to generate.
#' @param theta distribution parameter value. Default is 1.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p). If
#'   FALSE, probabilities p are given directly. Default is FALSE.
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are
#'   returned, otherwise \eqn{P(X > x)} is returned. Default is TRUE.
#'
#' @details
#' Probability density function (PDF)
#' \deqn{f(x\mid \theta )=\frac{\theta ^{2}}{(1+\theta )}(1+x)e^{-\theta x}}
#'
#' Cumulative distribution function (CDF)
#' \deqn{F(x\mid \theta ) =
#'   1 - \left(1+ \frac{\theta x}{1+\theta }\right)e^{-\theta x}}
#'
#' Quantile function (Inverse CDF)
#' \deqn{
#' Q(p\mid\theta) = -1 - \frac{1}{\theta}
#'   - \frac{1}{\theta} W_{-1}\!\left((1+\theta)(p-1)e^{-(1+\theta)}\right)
#' }
#'
#' where \eqn{W_{-1}()} is the negative branch of the Lambert W function.
#' 
#' The moment generating function (MGF) is:
#' \deqn{M_X(t)=\frac{\theta^2(\theta-t+1)}{(\theta+1)(\theta-t)^2}}
#' 
#' The distribution mean and variance are:
#' \deqn{\mu=\frac{\theta+2}{\theta(1+\theta)}}
#' \deqn{\sigma^2=\frac{\mu}{\theta+2}\left(\frac{6}{\theta}-4\right)-\mu^2}
#' 
#' @importFrom lamW lambertWm1
#'
#' @examples
#' x <- seq(0, 5, by = 0.1)
#' p <- seq(0.1, 0.9, by = 0.1)
#' q <- c(0.2, 3, 0.2)
#' dlindley(x, theta = 1.5)
#' dlindley(x, theta=0.5, log=TRUE)
#' plindley(q, theta = 1.5)
#' plindley(q, theta = 0.5, lower.tail = FALSE)
#' qlindley(p, theta = 1.5)
#' qlindley(p, theta = 0.5)
#' 
#' set.seed(123154)
#' rlindley(5, theta = 1.5)
#' rlindley(5, theta = 0.5)
#'
#' @rdname Lindley
#' @export
dlindley <- function(x, theta = 1, log = FALSE) {
  
  
  # --- Input Validation ---
  
  if (any(theta <= 0, na.rm = TRUE)) {
    warning("'theta' must be positive")
  }
  
  
  # --- Vectorization Setup ---
  
  n <- max(length(x), length(theta))
  x <- rep_len(x, n)
  theta <- rep_len(theta, n)
  
  
  # --- Compute in Log-Space ---
  
  # log(p) = 2*log(theta) - log(1+theta) + log(1+x) - theta*x
  log_p <- 2 * log(theta) - log(1 + theta) + log(1 + x) - theta * x
  
  # Handle invalid x (x must be > 0 for Lindley)
  invalid <- is.na(x) | x <= 0
  log_p[invalid] <- -Inf
  
  
  # --- Return ---
  
  if (log) {
    return(log_p)
  } else {
    return(exp(log_p))
  }
}


#' @rdname Lindley
#' @export
plindley <- function(q, theta = 1, lower.tail = TRUE, log.p = FALSE) {
  

  # --- Vectorization Setup ---
  n <- max(length(q), length(theta))
  q <- rep_len(q, n)
  theta <- rep_len(theta, n)
  
  # --- Compute CDF ---
  # F(x) = 1 - (1 + theta*x/(1+theta)) * exp(-theta*x)
  # This is already a closed-form expression - fully vectorizable
  
  cdf <- 1 - (1 + theta * q / (1 + theta)) * exp(-theta * q)
  
  # Handle boundaries
  cdf[q <= 0] <- 0
  cdf[!is.finite(q) & q > 0] <- 1
  
  # --- Apply Tail and Log Options ---
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  return(cdf)
}


#' @rdname Lindley
#' @export
qlindley <- function(p, theta = 1, lower.tail = TRUE, log.p = FALSE) {
  
  # Calculate the probability value (P[X <= q])
  if (log.p) {
    # 1. Handle log.p and lower.tail simultaneously using log(1 - exp(x))
    if (!lower.tail) {
      # p is log(P[X > q]). We need log(P[X <= q]) = log(1 - P[X > q])
      # Use log1p(x) = log(1+x) for stability. log(1 - exp(p)) is stable.
      p_val <- log(1 - exp(p)) 
    } else {
      # p is already log(P[X <= q])
      p_val <- p
    }
    
    # 2. Check for valid log-probabilities (must be <= 0)
    if (any(p_val > 0, na.rm = TRUE)) {
      warning("'p' must be <= 0 when log.p = TRUE and lower.tail = TRUE")
    }
    
    # 3. Convert to probability (0 < P < 1)
    p <- exp(p_val) 
    
  } else {
    # If not log.p, handle lower.tail directly
    if (!lower.tail) p <- 1 - p
    
    # Check for valid probabilities (0 <= P <= 1)
    if (any(p < 0 | p > 1, na.rm = TRUE)) {
      warning("'p' must be in [0, 1]")
    }
  }
  
  # --- Input Validation (Existing Checks) ---
  if (any(theta <= 0, na.rm = TRUE)) {
    warning("'theta' must be positive")
  }
  
  # --- Vectorization Setup ---
  n <- max(length(p), length(theta))
  p <- rep_len(p, n)
  theta <- rep_len(theta, n)
  
  # --- Compute Quantile (Rest of code is correct) ---
  # Argument to Lambert W
  w_arg <- (1 + theta) * (p - 1) * exp(-(1 + theta))
  
  # Lambert W (negative branch)
  w_val <- lamW::lambertWm1(w_arg)
  
  # Quantile
  q_val <- -1 - 1/theta - w_val/theta
  
  # Handle boundaries
  q_val[p == 0] <- 0
  q_val[p == 1] <- Inf
  
  return(q_val)
}

#' @rdname Lindley
#' @export
rlindley <- function(n, theta = 1) {
  
  if (length(n) != 1 || n < 0 || n != floor(n) || n==round(n)) {
    warning("'n' must be a non-negative integer")
  }
  
  # Recycle theta if needed
  theta <- rep_len(theta, n)
  
  # Generate via inverse CDF (fully vectorized)
  u <- runif(n)
  qlindley(u, theta = theta)
}
