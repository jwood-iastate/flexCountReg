#' Poisson-Lindley-Gamma (Negative Binomial-Lindley) Distribution
#'
#' These functions provide density, distribution function, quantile function, and random number generation for the Poisson-Lindley-Gamma (PLG) Distribution
#'
#' The Poisson-Lindley-Gamma is a count distribution that captures high densities for small integer values and provides flexibility for heavier tails.
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mean numeric value or vector of mean values for the distribution (the values have to be greater than 0).
#' @param theta single value or vector of values for the theta parameter of the distribution (the values have to be greater than 0).
#' @param alpha single value or vector of values for the `alpha` parameter of the gamma distribution in the special case that the mean = 1 and the variance = `alpha` (the values for `alpha` have to be greater than 0).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dplindGamma} computes the density (PDF) of the Poisson-Lindley-Gamma Distribution.
#'
#' \code{pplindGamma} computes the CDF of the Poisson-Lindley-Gamma Distribution.
#'
#' \code{qplindGamma} computes the quantile function of the Poisson-Lindley-Gamma Distribution.
#'
#' \code{rplindGamma} generates random numbers from the Poisson-Lindley-Gamma Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Poisson-Lindley-Gamma (PLG) distribution is:
#' \deqn{f(x|\mu,\theta,\alpha)=\frac{\alpha  (\theta+2) ^2 \Gamma (x+\alpha ) }{(\mu)^2(\theta +1)^3 \Gamma (\alpha )}\left(\frac{\mu\theta(\theta+1)}{\theta+2} U\left(x+1,2-\alpha ,\frac{\alpha  (\theta+2) }{\mu(\theta+1)}\right)+\alpha(x+1) U\left(x+2,3-\alpha ,\frac{\alpha  (\theta+2) }{\mu(\theta+1)}\right)\right)}
#' 
#' Where \eqn{\theta} is a distribution parameter from the Poisson-Lindley distribution with the restrictions that \eqn{\theta>0} a, \eqn{\alpha} is a parameter for the gamma distribution with the restriction \eqn{\alpha>0}, \eqn{mu} is the mean value, and \eqn{x} is a non-negative integer, and \deqn{U(a,b,z)} is the Tricomi's solution to the confluent hypergeometric function - also known as the confluent hypergeometric function of the second kind
#'
#' The expected value of the distribution is:
#' \deqn{E[x]=\mu}
#'
#' The variance is:
#' \deqn{\sigma^2=\mu+\left(\alpha+1-\frac{2}{(\theta+2)^2}\right)\mu^2}
#' 
#' While the distribution can be computed using the confluent hypergeometric function, that function has limitations in value it can be computed at (along with accuracy, in come cases). For this reason, the function uses Halton draws to perform simulation over the gamma distribution to solve the integral. This is sometimes more computationally efficient as well.
#'
#' @examples
#' dplindGamma(0, mean=0.75, theta=7, alpha=2)
#' pplindGamma(c(0,1,2,3,5,7,9,10), mean=0.75, theta=3, alpha=0.5)
#' qplindGamma(c(0.1,0.3,0.5,0.9,0.95), mean=1.67, theta=0.5, alpha=0.5)
#' rplindGamma(30, mean=0.5, theta=0.5, alpha=2)
#'
#' @importFrom stats runif
#' @importFrom gsl hyperg_U
#' @useDynLib flexCountReg
#' @name NegativeBinomialLindley
#' 
#' @rdname NegativeBinomialLindley
#' @export
dplindGamma <- Vectorize(function(x, mean=1, theta = 1, alpha=1, log=FALSE){
  
  if(mean <= 0 || theta <= 0 || alpha <= 0) warning('The values of `mean`, `theta`, and `alpha` all have to have values greater than 0.')
  
  U1 <- gsl::hyperg_U(x + 1, 2 - alpha, (alpha * (theta + 2)) / (mean * (theta + 1)))
  U2 <- gsl::hyperg_U(x + 2, 3 - alpha, (alpha * (theta + 2)) / (mean * (theta + 1)))
  co1 <- alpha*(theta+2)^2 * gamma(x+alpha) / (mean^2 * (theta+1)^3 * gamma(alpha))
  co2 <- mean*theta*(theta+1)/(theta+2)
  co3 <- alpha*(x+1)
  
  
  p <- co1*(co2*U1 + co3*U2)
  
  if (log) return(log(p))
  else return(p)
})

##' @rdname NegativeBinomialLindley
#' @export
pplindGamma <- function(q, mean = 1, theta = 1, alpha = 1,
                        lower.tail = TRUE, log.p = FALSE) {
  
  
  # --- Input Validation ---
  
  if (any(mean <= 0, na.rm = TRUE)) warning("'mean' must be positive")
  if (any(theta <= 0, na.rm = TRUE)) warning("'theta' must be positive")
  if (any(alpha <= 0, na.rm = TRUE)) warning("'alpha' must be positive")
  
  
  # --- Vectorization Setup ---
  
  n <- max(length(q), length(mean), length(theta), length(alpha))
  q <- rep_len(as.integer(floor(q)), n)
  mean <- rep_len(mean, n)
  theta <- rep_len(theta, n)
  alpha <- rep_len(alpha, n)
  
  
  # --- Initialize Result ---
  
  cdf <- rep(NA_real_, n)
  
  
  # --- Handle Special Cases ---
  
  invalid_q <- is.na(q) | q < 0
  cdf[invalid_q] <- 0
  
  
  # --- Helper: Safe log of Tricomi U function ---
  # gsl::hyperg_U can return 0 or Inf for extreme parameters
  
  log_hyperg_U_safe <- function(a, b, z) {
    u_val <- tryCatch(
      gsl::hyperg_U(a, b, z),
      warning = function(w) NA_real_,
      error = function(e) NA_real_
    )
    
    if (is.na(u_val) || !is.finite(u_val) || u_val <= 0) {
      return(-Inf)
    }
    return(log(u_val))
  }
  
  
  # --- Compute CDF for Valid Cases ---
  
  valid <- !invalid_q
  
  if (any(valid)) {
    for (i in which(valid)) {
      mean_i <- mean[i]
      theta_i <- theta[i]
      alpha_i <- alpha[i]
      q_i <- q[i]
      
      # Pre-compute constants
      # PMF uses Tricomi's confluent hypergeometric function U(a, b, z)
      # z = alpha * (theta + 2) / (mean * (theta + 1))
      z_arg <- alpha_i * (theta_i + 2) / (mean_i * (theta_i + 1))
      
      # Common coefficient parts (in log-space)
      # co1 = alpha * (theta+2)^2 * Gamma(x+alpha) / (mean^2 * (theta+1)^3 * Gamma(alpha))
      log_co1_base <- log(alpha_i) + 2 * log(theta_i + 2) -
        2 * log(mean_i) - 3 * log(theta_i + 1) - lgamma(alpha_i)
      
      # co2 = mean * theta * (theta+1) / (theta+2)
      co2 <- mean_i * theta_i * (theta_i + 1) / (theta_i + 2)
      log_co2 <- log(co2)
      
      # Compute log-PMF for 0:q_i
      log_pmf <- numeric(q_i + 1)
      
      for (y in 0:q_i) {
        # log(Gamma(y + alpha))
        log_gamma_y_alpha <- lgamma(y + alpha_i)
        
        # U1 = U(y + 1, 2 - alpha, z)
        # U2 = U(y + 2, 3 - alpha, z)
        log_U1 <- log_hyperg_U_safe(y + 1, 2 - alpha_i, z_arg)
        log_U2 <- log_hyperg_U_safe(y + 2, 3 - alpha_i, z_arg)
        
        # co3 = alpha * (y + 1)
        log_co3 <- log(alpha_i) + log(y + 1)
        
        # PMF = co1 * Gamma(y+alpha) * (co2 * U1 + co3_val * U2)
        # In log-space: need to compute log(co2 * U1 + alpha*(y+1) * U2)
        
        # This is tricky - we need log(a + b) where a = co2 * U1, b = alpha*(y+1) * U2
        log_term1 <- log_co2 + log_U1
        log_term2 <- log_co3 + log_U2
        
        # Log-sum-exp for (term1 + term2)
        if (!is.finite(log_term1) && !is.finite(log_term2)) {
          log_sum_terms <- -Inf
        } else if (!is.finite(log_term1)) {
          log_sum_terms <- log_term2
        } else if (!is.finite(log_term2)) {
          log_sum_terms <- log_term1
        } else {
          max_term <- max(log_term1, log_term2)
          log_sum_terms <- max_term + log(exp(log_term1 - max_term) + exp(log_term2 - max_term))
        }
        
        log_pmf[y + 1] <- log_co1_base + log_gamma_y_alpha + log_sum_terms
      }
      
      # Log-sum-exp for CDF
      finite_mask <- is.finite(log_pmf)
      if (!any(finite_mask)) {
        cdf[i] <- NA_real_
      } else {
        max_log <- max(log_pmf[finite_mask])
        log_cdf <- max_log + log(sum(exp(log_pmf[finite_mask] - max_log)))
        cdf[i] <- exp(log_cdf)
      }
    }
  }
  
  
  # --- Apply Tail and Log Options ---
  
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  return(cdf)
}

#' @rdname NegativeBinomialLindley
#' @export
qplindGamma <- Vectorize(function(p, mean=1, theta=1, alpha=1) {
  if(p < 0) warning("The value of `p` must be a value greater than 0 and less than 1.")
  if(is.na(p)) warning("The value of `p` cannot be an `NA` value")
  
  if(mean<=0 || theta<=0 || alpha<=0) warning('The values of `mean`, `theta`, and `alpha` all have to have values greater than 0.')
  
  
  y <- 0
  p_value <- max(pplindGamma(y, mean, theta, alpha=alpha), .Machine$double.xmin)
  while(p_value < p){
    y <- y + 1
    p_value_new <- max(pplindGamma(y, mean, theta, alpha=alpha), .Machine$double.xmin)
    if (!is.na(p_value_new)) p_value <- p_value_new else break
  }
  return(y)
})


#' @rdname NegativeBinomialLindley
#' @export
rplindGamma <- function(n, mean=1, theta=1, alpha=1) {
  
  if(mean<=0 || theta<=0  || alpha<=0) warning('The values of `mean`, `theta`, and `alpha` all have to have values greater than 0.')
  
  u <- runif(n)
  y <- lapply(u, function(p) qplindGamma(p, mean, theta, alpha=alpha))
  return(unlist(y))
}
