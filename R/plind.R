#' Poisson-Lindley Distribution
#'
#' These functions provide density, distribution function, quantile function,
#' and random number generation for the Poisson-Lindley (PL) Distribution
#'
#' The Poisson-Lindley is a 2-parameter count distribution that captures high
#' densities for small integer values. This makes it ideal for data that are
#' zero-inflated.
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mean numeric value or vector of mean values for the distribution (the
#'   values have to be greater than 0).
#' @param theta single value or vector of values for the theta parameter of the
#'   distribution (the values have to be greater than 0).
#' @param lambda alternative parameterization (use instead of the mean); numeric
#'   value or vector of values for lambda parameter of the distribution (the
#'   values have to be greater than 0).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dplind} computes the density (PDF) of the Poisson-Lindley Distribution.
#'
#' \code{pplind} computes the CDF of the Poisson-Lindley Distribution.
#'
#' \code{qplind} computes the quantile function of the Poisson-Lindley
#' Distribution.
#'
#' \code{rplind} generates random numbers from the Poisson-Lindley Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Poisson-Lindley (PL)
#' distribution is:
#' \deqn{f(y| \theta, \lambda) = 
#'   \frac{\theta^2 \lambda^y (\theta + \lambda + y + 1)}
#'     {(\theta + 1)(\theta + \lambda)^{y + 2}}}
#'
#' Where \eqn{\theta} and \eqn{\lambda} are distribution parameters with the
#' restrictions that \eqn{\theta > 0} and \eqn{\lambda > 0}, and \eqn{y} is a
#' non-negative integer.
#'
#' The expected value of the distribution is:
#' \deqn{\mu=\frac{\lambda(\theta + 2)}{\theta(\theta + 1)}}
#'
#' The default is to use the input mean value for the distribution. However, the
#' lambda parameter can be used as an alternative to the mean value.
#' 
#' @returns dplind gives the density, pplind gives the distribution 
#'  function, qplind gives the quantile function, and rplind generates
#'  random  deviates.
#' 
#'  The length of the result is determined by n for rplind, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' dplind(0, mean = 0.75, theta = 7)
#' pplind(c(0, 1, 2, 3, 5, 7, 9, 10), mean = 0.75, theta = 7)
#' qplind(c(0.1, 0.3, 0.5, 0.9, 0.95), lambda = 4.67, theta = 7)
#' rplind(30, mean = 0.75, theta = 7)
#'
#' @importFrom stats runif
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @export
#' @name PoissonLindley



msg1 <- ("The value of `x` must be a non-negative whole number")
msg2 <- 
  ('The values of `mean` and `theta` both have to have values greater than 0.')

#' @rdname PoissonLindley
#' @export
dplind <- Vectorize(function(
    x, mean = 1, theta = 1, lambda = NULL, log = FALSE){
  
  # Test to make sure the value of x is an integer
  tst <- !(is.na(nchar(strsplit(as.character(x), "\\.")[[1]][2]) > 0))
  if(tst | x < 0) warning(msg1)
  if(is.null(lambda)){
    if(mean <= 0 | theta <= 0){
      warning(msg2)
    }
    else{
      lambda <- mean * theta * (theta + 1) / (theta + 2)
      p <- (theta^2 * lambda^x * (theta + lambda + x + 1))/
        ((theta + 1) * (theta + lambda)^(x + 2))
    }
  }
  else{
    if(lambda <= 0 | theta <= 0){
      warning(msg2)
    }
    else{
      p <- dplind_cpp(x, mean, theta)
    }
  }
  if (log) return(log(p))
  else return(p)
})

#' @rdname PoissonLindley
#' @export
pplind <- function(q, mean = 1, theta = 1, lambda = NULL,
                   lower.tail = TRUE, log.p = FALSE) {
  
  # --- Input Validation ---
  if (is.null(lambda)) {
    if (any(mean <= 0, na.rm = TRUE) || any(theta <= 0, na.rm = TRUE)) 
      warning("'mean' and 'theta' must be positive")
  } else {
    if (any(lambda <= 0, na.rm = TRUE) || any(theta <= 0, na.rm = TRUE)) 
      warning("'lambda' and 'theta' must be positive")
    # Convert lambda to mean
    mean <- lambda * (theta + 2) / (theta * (theta + 1))
  }
  
  # --- Vectorization Setup ---
  n <- max(length(q), length(mean), length(theta))
  q <- rep_len(as.integer(floor(q)), n)
  mean <- rep_len(mean, n)
  theta <- rep_len(theta, n)
  
  # --- Initialize Result ---
  cdf <- rep(NA_real_, n)
  
  # --- Handle Special Cases ---
  # q < 0 -> CDF = 0
  invalid_q <- is.na(q) | q < 0
  cdf[invalid_q] <- 0
  
  # --- Compute CDF for Valid Cases ---
  valid <- !invalid_q
  
  if (any(valid)) {
    # Get unique parameter combinations
    param_df <- data.frame(
      idx = which(valid),
      q = q[valid],
      mean = mean[valid],
      theta = theta[valid]
    )
    
    unique_params <- unique(param_df[, c("mean", "theta")])
    
    for (i in seq_len(nrow(unique_params))) {
      m_i <- unique_params$mean[i]
      t_i <- unique_params$theta[i]
      
      mask <- param_df$mean == m_i & param_df$theta == t_i
      indices <- param_df$idx[mask]
      q_vals <- param_df$q[mask]
      max_q <- max(q_vals)
      
      # Pre-compute lambda
      lambda_i <- m_i * t_i * (t_i + 1) / (t_i + 2)
      
      # Compute log-PMF for 0:max_q
      y_seq <- 0:max_q
      log_pmf <- 2 * log(t_i) + 
        y_seq * log(lambda_i) + 
        log(t_i + lambda_i + y_seq + 1) -
        log(t_i + 1) - 
        (y_seq + 2) * log(t_i + lambda_i)
      
      # Compute cumulative sum using log-sum-exp for stability
      cum_log_prob <- numeric(length(y_seq))
      cum_log_prob[1] <- log_pmf[1]
      
      # --- FIX START: Check length before looping ---
      if (length(y_seq) > 1) {
        for (k in 2:length(y_seq)) {
          a <- cum_log_prob[k - 1]
          b <- log_pmf[k]
          if (a >= b) {
            cum_log_prob[k] <- a + log1p(exp(b - a))
          } else {
            cum_log_prob[k] <- b + log1p(exp(a - b))
          }
        }
      }
      # --- FIX END ---
      
      # Assign CDF values to appropriate indices
      for (j in seq_along(indices)) {
        idx_j <- indices[j]
        q_j <- q_vals[j]
        cdf[idx_j] <- exp(cum_log_prob[q_j + 1])  # +1 because R is 1-indexed
      }
    }
  }
  
  # --- Apply Tail and Log Options ---
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  return(cdf)
}

#' @rdname PoissonLindley
#' @export
qplind <- Vectorize(function(p, mean=1, theta=1, lambda=NULL) {
  if (is.null(lambda)){
    if(mean<=0 || theta<=0){
      stop(msg2)
    }
  }
  else{
    if (lambda<=0 | theta <= 0) {
      warning(msg2)
    }
    else{
      mean <- lambda * (theta + 2) / (theta * (theta + 1))
    }
  }

  y <- 0
  p_value <- pplind(q=y, mean=mean, theta=theta)
  while(p_value < p){
    y <- y + 1
    p_value <- pplind(q=y, mean=mean, theta=theta)
  }
  return(y)
})


#' @rdname PoissonLindley
#' @export
rplind <- function(n, mean=1, theta=1, lambda=NULL) {
  
  if(is.null(lambda)){
    if(mean <= 0 | theta <= 0){
      warning(msg2)
    }
  }
  else{
    if(lambda <= 0 | theta <= 0){
      warning(msg2)
    }
    else{
      mean <- lambda * (theta + 2) / (theta * (theta + 1))
    }
  }

  u <- runif(n)
  y <- vapply(X = u, FUN = \(p) qplind(p, mean, theta), FUN.VALUE = numeric(1))
  return(y)
}

