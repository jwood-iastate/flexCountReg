#' Poisson-Weibull Distribution Functions
#'
#' These functions provide density, distribution function, quantile function, and random number 
#' generation for the Poisson-Weibull Distribution, which is specified either by its shape and scale 
#' parameters or by its mean and standard deviation.
#'
#' The Poisson-Weibull distribution uses the Weibull distribution as a mixing distribution for a 
#' Poisson process. It is useful for modeling overdispersed count data. The density function 
#' (probability mass function) for the Poisson-Weibull distribution is given by:
#' \deqn{P(y|\lambda,\alpha,\beta) = \int_0^\infty \frac{e^{-\lambda} \lambda^y }{y!} f(\theta; \alpha, \beta) d\theta}
#' where \eqn{f(\theta; \alpha, \beta)} is the PDF of the Weibull distribution and \eqn{\lambda} is the mean of the Poisson distribution.
#'
#' @param x A numeric value or vector of values for which the PDF or CDF is calculated.
#' @param q Quantile or a vector of quantiles.
#' @param p A numeric value or vector of probabilities for the quantile function.
#' @param n The number of random samples to generate.
#' @param alpha Shape parameter of the Weibull distribution (optional if mean and sd are provided).
#' @param beta Scale parameter of the Weibull distribution (optional if mean and sd are provided).
#' @param mean_value Mean of the Weibull distribution (optional if alpha and beta are provided).
#' @param sd_value Standard deviation of the Weibull distribution (optional if alpha and beta are provided).
#' @param lambda Mean value of the Poisson distribution.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE, probabilities are P[X <= x], otherwise P[X > x].
#' @param ndraws the number of Halton draws to use for the integration.
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
#' dpoisweibull(4, lambda=1.5, mean_value=1.5, sd_value=0.5, ndraws=10)
#' ppoisweibull(4, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#' qpoisweibull(0.95, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#' rpoisweibull(10, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#'
#' @import stats randtoolbox
#' @export
#' @name PoissonWeibull

#' @rdname PoissonWeibull
#' @export
dpoisweibull <- Vectorize(function(x, lambda, alpha = NULL, beta = NULL, 
                                   mean_value = NULL, sd_value = NULL, 
                                   ndraws=1500, log = FALSE) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    beta <- params[2]
  }
  
  if (is.null(alpha) || is.null(beta)) stop("Parameters alpha and beta are required")
  
  # Generate Halton draws to use as quantile values
  h <- randtoolbox::halton(ndraws)
  
  # Evaluate the density of the normal distribution at those quantiles and use the exponent to transform to lognormal values
  weibullDist <- stats::qweibull(h, shape=alpha, scale=beta)
  mu_i <- lambda * weibullDist
  
  p_pweib.i <- sapply(mu_i, stats::dpois, x=x)
  
  p <- mean(p_pweib.i)
  
  if (log) return(log(p))
  else return(p)
})

#' @rdname PoissonWeibull
#' @export
ppoisweibull <- Vectorize(function(q, lambda, alpha = NULL, beta = NULL, 
                                   mean_value = NULL, sd_value = NULL, 
                                   ndraws=1500, lower.tail = TRUE, log.p = FALSE) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    beta <- params[2]
  }
  
  if (is.null(alpha) || is.null(beta)) stop("Parameters alpha and beta are required")
  
  y <- seq(0,q,1)
  probs <- dpoisweibull(y, lambda=lambda, alpha = alpha, beta = beta, ndraws=ndraws)
  p <- sum(probs)
  
  if(!lower.tail) p <- 1-p
  
  if (log.p) return(log(p))
  else return(p)
  
  return(result)
})

#' @rdname PoissonWeibull
#' @export
qpoisweibull <- Vectorize(function(p, lambda, alpha = NULL, beta = NULL, 
                                   mean_value = NULL, sd_value = NULL, 
                                   ndraws=1500) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    beta <- params[2]
  }
  y <- 0
  p_value <- ppoisweibull(y, lambda=lambda, alpha = alpha, beta = beta, ndraws=ndraws)
  while (p_value < p){
    y <- y + 1
    p_value <- ppoisweibull(y, lambda=lambda, alpha = alpha, beta = beta, ndraws=ndraws)
  }
  return(y)
})

#' @rdname PoissonWeibull
#' @export
rpoisweibull <- function(n, lambda, alpha = NULL, beta = NULL, 
                         mean_value = NULL, sd_value = NULL, ndraws=1500) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    beta <- params[2]
  }
  
  if (is.null(alpha) || is.null(beta)) stop("Parameters alpha and beta are required")
  
  u <- stats::runif(n)
  y <- qpoisweibull(u, lambda=lambda, alpha = alpha, beta = beta, ndraws=ndraws)
  return(y)
}


# Helper function to calculate alpha and beta from mean and standard deviation
calculate_params <- function(mean_value = NULL, sd_value = NULL) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    # Function to minimize
    objective_function <- function(alpha, mean_value, sd_value) {
      beta <- mean_value / gamma(1 + 1/alpha)
      calculated_mean <- beta * gamma(1 + 1/alpha)
      calculated_var <- beta^2 * (gamma(1 + 2/alpha) - (gamma(1 + 1/alpha))^2)
      (calculated_mean - mean_value)^2 + (sqrt(calculated_var) - sd_value)^2
    }
    
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