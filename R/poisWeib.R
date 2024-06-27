#' Poisson-Weibull Distribution Functions
#'
#' These functions provide density, distribution function, quantile function, and random number 
#' generation for the Poisson-Weibull Distribution, which is specified either by its shape and scale 
#' parameters or by its mean and standard deviation.
#'
#' The Poisson-Weibull distribution uses the Weibull distribution as a mixing distribution for a 
#' Poisson process. It is useful for modeling overdispersed count data. The density function 
#' (probability mass function) for the Poisson-Weibull distribution is given by:
#' \deqn{P(X = k) = \int_0^\infty e^{-\theta} \theta^k / k! f(\theta; \alpha, \beta) d\theta,}
#' where \eqn{f(\theta; \alpha, \beta)} is the PDF of the Weibull distribution.
#'
#' @param x A numeric value or vector of values for which the PDF or CDF is calculated.
#' @param p A numeric value or vector of probabilities for the quantile function.
#' @param n The number of random samples to generate.
#' @param alpha Shape parameter of the Weibull distribution (optional if mean and sd are provided).
#' @param beta Scale parameter of the Weibull distribution (optional if mean and sd are provided).
#' @param mean_value Mean of the Weibull distribution (optional if alpha and beta are provided).
#' @param sd_value Standard deviation of the Weibull distribution (optional if alpha and beta are provided).
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE, probabilities are P[X <= x], otherwise P[X > x].
#'
#' @details
#' \code{dpoisweibull} computes the density of the Poisson-Weibull distribution.
#'
#' \code{ppoisweibull} computes the distribution function of the Poisson-Weibull distribution.
#'
#' \code{qpoisweibull} computes the quantile function of the Poisson-Weibull distribution.
#'
#' \code{rpoisweibull} generates random numbers following the Poisson-Weibull distribution.
#'
#' The shape and scale parameters directly define the Weibull distribution, whereas the mean and 
#' standard deviation are used to compute these parameters indirectly.
#'
#' @examples
#' dpoisweibull(4, alpha = 1.5, beta = 0.5)
#' ppoisweibull(4, alpha = 1.5, beta = 0.5)
#' qpoisweibull(0.95, alpha = 1.5, beta = 0.5)
#' rpoisweibull(10, alpha = 1.5, beta = 0.5)
#'
#' @importFrom stats runif optimize
#' @export
#' @name PoissonWeibull

#' @rdname PoissonWeibull
#' @export
dpoisweibull <- function(x, alpha = NULL, beta = NULL, mean_value = NULL, sd_value = NULL, log = FALSE) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    beta <- params[2]
  }
  
  if (is.null(alpha) || is.null(beta)) stop("Parameters alpha and beta are required")
  
  result <- sapply(x, function(x_val) {
    integrand <- function(theta) {
      # Use log to calculate the expression for better numerical stability
      log_pmf <- -theta + x_val * log(theta) - lgamma(x_val + 1) +
        log(alpha) - alpha * log(beta) + (alpha - 1) * log(theta) - (theta/beta)^alpha
      exp(log_pmf)  # Return to original scale
    }
    integral <- integrate(integrand, lower = 0, upper = Inf, rel.tol = .Machine$double.eps^0.5)$value
    if (log) log(integral) else integral
  })
  
  return(result)
}

#' @rdname PoissonWeibull
#' @export
ppoisweibull <- function(q, alpha = NULL, beta = NULL, mean_value = NULL, sd_value = NULL, lower.tail = TRUE, log.p = FALSE) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    beta <- params[2]
  }
  
  if (is.null(alpha) || is.null(beta)) stop("Parameters alpha and beta are required")
  
  result <- sapply(q, function(q_val) {
    cdf_val <- sum(sapply(0:q_val, function(k) dpoisweibull(k, alpha = alpha, beta = beta)))
    if (!lower.tail) cdf_val <- 1 - cdf_val
    if (log.p) log(cdf_val) else cdf_val
  })
  
  return(result)
}

#' @rdname PoissonWeibull
#' @export
qpoisweibull <- function(p, alpha = NULL, beta = NULL, mean_value = NULL, sd_value = NULL) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    beta <- params[2]
  }
  
  if (is.null(alpha) || is.null(beta)) stop("Parameters alpha and beta are required")
  
  uniroot(function(x) ppoisweibull(x, alpha = alpha, beta = beta) - p, lower = 0, upper = 10^5)$root
}

#' @rdname PoissonWeibull
#' @export
rpoisweibull <- function(n, alpha = NULL, beta = NULL, mean_value = NULL, sd_value = NULL) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    beta <- params[2]
  }
  
  if (is.null(alpha) || is.null(beta)) stop("Parameters alpha and beta are required")
  
  theta <- rweibull(n, shape = alpha, scale = beta)
  rpois(n, lambda = theta)
}


# Helper function to calculate alpha and beta from mean and standard deviation
calculate_params <- function(mean_value = NULL, sd_value = NULL, alpha = NULL, beta = NULL) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    # Function to minimize
    objective_function <- function(alpha, mean_value, sd_value) {
      beta <- mean_value / gamma(1 + 1/alpha)
      calculated_mean <- beta * gamma(1 + 1/alpha)
      calculated_var <- beta^2 * (gamma(1 + 2/alpha) - (gamma(1 + 1/alpha))^2)
      (calculated_mean - mean_value)^2 + (sqrt(calculated_var) - sd_value)^2
    }
    
    # Initial guess for alpha
    alpha_guess <- 1.5
    res <- optimize(objective_function, c(0.1, 5), mean_value = mean_value, sd_value = sd_value)
    alpha <- res$minimum
    beta <- mean_value / gamma(1 + 1/alpha)
    c(alpha, beta)
  } else if (!is.null(alpha) && !is.null(beta)) {
    c(alpha, beta)
  } else {
    stop("Insufficient parameters: Provide either mean and standard deviation, or alpha and beta.")
  }
}
