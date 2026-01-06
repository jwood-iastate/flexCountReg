#' Poisson-Weibull Distribution Functions
#'
#' These functions provide density, distribution function, quantile function,
#' and random number generation for the Poisson-Weibull Distribution, which is
#' specified either by its shape and scale parameters or by its mean and
#' standard deviation.
#'
#' The Poisson-Weibull distribution uses the Weibull distribution as a mixing
#' distribution for a Poisson process. It is useful for modeling overdispersed
#' count data. The density function (probability mass function) for the
#' Poisson-Weibull distribution is given by:
#' \deqn{P(y|\lambda,\alpha,\sigma) = 
#' \int_0^\infty \frac{e^{-\lambda x} \lambda^y x^y }{y!} 
#'   \left(\frac{\alpha}{\sigma}\right) 
#'   \left(\frac{x}{\sigma}\right)^{\alpha-1}
#'   e^{-\left(\frac{x}{\sigma}\right)^\alpha} dx}
#' where \eqn{f(x| \alpha, \sigma)} is the PDF of the Weibull distribution and
#' \eqn{\lambda} is the mean of the Poisson distribution.
#'
#' @param x A numeric value or vector of values for which the PDF or CDF is
#'   calculated.
#' @param q Quantile or a vector of quantiles.
#' @param p A numeric value or vector of probabilities for the quantile
#'   function.
#' @param n The number of random samples to generate.
#' @param alpha Shape parameter of the Weibull distribution (optional if mean
#'   and sd are provided).
#' @param sigma Scale parameter of the Weibull distribution (optional if mean
#'   and sd are provided).
#' @param mean_value Mean of the Weibull distribution (optional if alpha and
#'   sigma are provided).
#' @param sd_value Standard deviation of the Weibull distribution (optional if
#'   alpha and sigma are provided).
#' @param lambda Mean value of the Poisson distribution.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @param log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE, probabilities are P[X <= x], otherwise
#'   P[X > x].
#' @param ndraws the number of Halton draws to use for the integration.
#'
#' @details
#' \code{dpoisweibull} computes the density of the Poisson-Weibull distribution.
#'
#' \code{ppoisweibull} computes the distribution function of the Poisson-Weibull
#' distribution.
#'
#' \code{qpoisweibull} computes the quantile function of the Poisson-Weibull
#' distribution.
#'
#' \code{rpoisweibull} generates random numbers following the Poisson-Weibull
#' distribution.
#'
#' The shape and scale parameters directly define the Weibull distribution,
#' whereas the mean and standard deviation are used to compute these parameters
#' indirectly.
#' 
#' @details dpoisweibull gives the density, ppoisweibull gives the distribution 
#'  function, qpoisweibull gives the quantile function, and rpoisweibull 
#'  generates random  deviates.
#' 
#'  The length of the result is determined by n for rpoisweibull, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' dpoisweibull(4, lambda=1.5, mean_value=1.5, sd_value=0.5, ndraws=10)
#' ppoisweibull(4, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#' qpoisweibull(0.95, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#' rpoisweibull(10, lambda=1.5, mean_value=1.5, sd_value=2, ndraws=10)
#'
#' @importFrom stats runif
#' @importFrom randtoolbox halton
#' @export
#' @name PoissonWeibull
#' 
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @rdname PoissonWeibull
#' @export
dpoisweibull <- function(x, lambda = NULL, alpha = NULL, sigma = NULL, 
                         mean_value = NULL, sd_value = NULL, 
                         ndraws = 1500, log = FALSE) {
  
  # --- Input Parameter Logic ---
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    sigma <- params[2]
  }
  
  if (is.null(alpha) || is.null(sigma)) 
    warning("Parameters alpha and sigma are required")
  
  if (is.null(lambda)) {
    if (!is.null(mean_value)) {
      lambda <- mean_value / (sigma * gamma(1 + 1 / alpha))
    } else {
      warning("Parameter lambda or mean_value is required")
    }
  }
  
  # --- Vectorization Setup ---
  # Ensure parameters are vectors of appropriate length matches x
  n <- length(x)
  if(length(lambda) == 1) lambda <- rep(lambda, n)
  if(length(alpha) == 1) alpha <- rep(alpha, n)
  if(length(sigma) == 1) sigma <- rep(sigma, n)
  
  # --- Optimization ---
  # Generate Halton draws ONCE for the whole batch
  h <- randtoolbox::halton(ndraws)
  
  # Assuming dpWeib_cpp is vectorized or we map over it. 
  # Ideally, dpWeib_cpp should accept vectors. If strictly scalar C++, 
  # we must use mapply here, but passing 'h' prevents regenerating it.
  
  # If dpWeib_cpp is vectorized in C++:
  p <- dpWeib_cpp(x, lambda, alpha, sigma, h)
  
  if (log) return(log(p))
  else return(p)
}

#' @rdname PoissonWeibull
#' @export
ppoisweibull <- Vectorize(function(
    q, lambda=NULL, alpha = NULL, sigma = NULL, 
    mean_value = NULL, sd_value = NULL, 
    ndraws=1500, lower.tail = TRUE, log.p = FALSE) {
  
  # Parameter logic repeated ensures safety inside Vectorize
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    sigma <- params[2]
  }
  if (is.null(lambda) && !is.null(mean_value)){
    lambda <- mean_value/(sigma*gamma(1+1/alpha))
  }
  
  if (q < 0) return(if(log.p) -Inf else 0)
  
  # Efficiency: dpoisweibull is now faster
  y <- 0:floor(q)
  
  # Note: dpoisweibull now handles the heavy lifting
  probs <- 
    dpoisweibull(y, lambda=lambda, alpha=alpha, sigma=sigma, ndraws=ndraws)
  p <- sum(probs)
  
  if(!lower.tail) p <- 1-p
  if (log.p) return(log(p))
  else return(p)
})

#' @rdname PoissonWeibull
#' @export
qpoisweibull <- Vectorize(function(p, lambda=NULL, alpha = NULL, sigma = NULL, 
                                   mean_value = NULL, sd_value = NULL, 
                                   ndraws=1500) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    sigma <- params[2]
  }
  if (is.null(alpha) || is.null(sigma)) 
    warning("Parameters alpha and sigma are required")
  if (is.null(lambda) && !is.null(mean_value)){
    lambda <- mean_value/(sigma*gamma(1+1/alpha))
  } else if (is.null(lambda) && is.null(mean_value)) {
    warning("Parameter lambda or mean_value is required")
  }
  y <- 0
  p_value <- 
    ppoisweibull(y, lambda=lambda, alpha = alpha, sigma = sigma, ndraws=ndraws)
  while (p_value < p){
    y <- y + 1
    p_value <- ppoisweibull(
      y, lambda = lambda, alpha = alpha, sigma = sigma, ndraws = ndraws)
  }
  return(y)
})

#' @rdname PoissonWeibull
#' @export
rpoisweibull <- function(n, lambda=NULL, alpha = NULL, sigma = NULL, 
                         mean_value = NULL, sd_value = NULL, ndraws=1500) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    params <- calculate_params(mean_value = mean_value, sd_value = sd_value)
    alpha <- params[1]
    sigma <- params[2]
  }
  if (is.null(alpha) || is.null(sigma)) 
    warning("Parameters alpha and sigma are required")
  if (is.null(lambda) && !is.null(mean_value)){
    lambda <- mean_value/(sigma*gamma(1+1/alpha))
  } else if (is.null(lambda) && is.null(mean_value)) {
    warning("Parameter lambda or mean_value is required")
  }
  
  if (is.null(alpha) || is.null(sigma)) 
    warning("Parameters alpha and sigma are required")
  
  u <- stats::runif(n)
  y <- 
    qpoisweibull(u, lambda=lambda, alpha = alpha, sigma = sigma, ndraws=ndraws)
  return(y)
}


# Helper function to calculate alpha and sigma from mean and standard deviation
calculate_params <- function(mean_value = NULL, sd_value = NULL) {
  if (!is.null(mean_value) && !is.null(sd_value)) {
    # Function to minimize
    objective_function <- function(alpha, mean_value, sd_value) {
      sigma <- mean_value / gamma(1 + 1/alpha)
      calculated_mean <- sigma * gamma(1 + 1/alpha)
      calculated_var <- sigma^2 * (gamma(1 + 2/alpha) - (gamma(1 + 1/alpha))^2)
      (calculated_mean - mean_value)^2 + (sqrt(calculated_var) - sd_value)^2
    }
    
    res <- 
      optimize(objective_function, 
               c(0.1, 5), 
               mean_value = mean_value, 
               sd_value = sd_value)
    alpha <- res$minimum
    sigma <- mean_value / gamma(1 + 1/alpha)
    c(alpha, sigma)
  } else if (!is.null(alpha) && !is.null(sigma)) {
    c(alpha, sigma)
  } else {
    msg <- paste(
      "Insufficient parameters: Provide either mean and standard deviation,",
      "or alpha and sigma."
    )
    warning(msg)
  }
}