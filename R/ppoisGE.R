#' Poisson-Generalized-Exponential Distribution
#'
#' These functions provide density, distribution function, quantile function,
#' and random number generation for the Poisson-Generalized-Exponential (PGE)
#' Distribution
#'
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mean numeric value or vector of mean values for the distribution (the
#'   values have to be greater than 0). This is NOT the value of \eqn{\lambda}.
#' @param shape numeric value or vector of shape values for the shape parameter
#'   of the generalized exponential distribution (the values have to be greater
#'   than 0).
#' @param scale single value or vector of values for the scale parameter of the
#'   generalized exponential distribution (the values have to be greater than
#'   0).
#' @param ndraws the number of Halton draws to use for the integration.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#' @param haltons an optional vector of Halton draws to use instead of ndraws.
#'
#' @details
#' \code{dpge} computes the density (PDF) of the PGE Distribution.
#'
#' \code{ppge} computes the CDF of the PGE Distribution.
#'
#' \code{qpge} computes the quantile function of the PGE Distribution.
#'
#' \code{rpge} generates random numbers from the PGE Distribution.
#' 
#' The Generalized Exponential distribution can be written as a function with a
#' shape parameter \eqn{\alpha>0} and scale parameter \eqn{\gamma>0}. The
#' distribution has strictly positive continuous values. The PDF of the
#' distribution is:
#' \deqn{f(x|\alpha,\gamma)=\frac{\alpha}{\gamma}\left(1-e^{-\frac{x}{\gamma}}\right)^{\alpha-1}e^{-\frac{x}{\gamma}}} 
#' 
#' Thus, the compound Probability Mass Function(PMF) for the PGE distribution
#' is:
#' \deqn{f(y|\lambda,\alpha,\beta)=\int_0^\infty \frac{\lambda^y x^y e^{-\lambda x}}{y!}\frac{\alpha}{\gamma}\left(1-e^{-\frac{x}{\gamma}}\right)^{\alpha-1}e^{-\frac{x}{\gamma}} dx}
#' 
#' The expected value of the distribution is:
#' \deqn{E[y]=\mu=\lambda \left(\frac{\psi(\alpha+1)-\psi(1)}{\gamma}\right)}
#' 
#' Where \eqn{\psi(\cdot)} is the digamma function.
#' 
#' The variance is:
#' \deqn{\sigma^2=\lambda \left(\frac{\psi(\alpha+1)-\psi(1)}{\gamma}\right) + \left(\frac{-\psi'(\alpha+1)+\psi'(1)}{\gamma^2}\right)\lambda^2}
#' 
#' Where \eqn{\psi'(\cdot)} is the trigamma function.
#' 
#' To ensure that \eqn{\mu=e^{X\beta}}, \eqn{\lambda} is replaced with:
#' \deqn{\lambda=\frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)}}
#' 
#' This results in:
#' \deqn{f(y|\mu,\alpha,\beta)=\int_0^\infty \frac{\left(\frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)}\right)^y x^y e^{-\left(\frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)}\right) x}}{y!}\frac{\alpha}{\gamma}\left(1-e^{-\frac{x}{\gamma}}\right)^{\alpha-1}e^{-\frac{x}{\gamma}} dx}
#' 
#' Halton draws are used to perform simulation over the lognormal distribution to solve the integral.
#'
#' @references
#' Gupta, R. D., & Kundu, D. (2007). Generalized exponential distribution: Existing results and some recent developments. Journal of Statistical planning and inference, 137(11), 3537-3547.
#'
#' @examples
#' dpge(0, mean=0.75, shape=2, scale=1, ndraws=2000)
#' ppge(c(0,1,2,3,4,5,6), mean=0.75, shape=2, scale=1, ndraws=500)
#' qpge(c(0.1,0.3,0.5,0.9,0.95), mean=0.75, shape=2, scale=1, ndraws=500)
#' rpge(30, mean=0.75,  shape=2, scale=1, ndraws=500)
#'
#' @importFrom stats runif
#' @importFrom randtoolbox halton
#' @export
#' @name PoissonGeneralizedExponential

#' @rdname PoissonGeneralizedExponential
#' @export
dpge <- Vectorize(function(x, mean=1, shape=1, scale=1, ndraws=1500, log=FALSE, haltons=NULL){
  
  if(mean<=0 || scale<=0 || shape <=0) warning('The values of `mean`, `shape`, and `scale` have to have values greater than 0.')
  
  # qge <- function(p, shape, scale){
  #   q <- log((p*shape+1)^(1/shape)+1)/scale
  #   return(q)
  # }
  # 
  lambda <- mean * scale /(digamma(shape+1)-digamma(1))

  # Generate Halton draws to use as quantile values
  if (!is.null(haltons)) h <- haltons else h <- randtoolbox::halton(ndraws)

  # Evaluate the density of the normal distribution at those quantiles and use the exponent to transform to lognormal values
  gedist <- log(1-h^(1/shape))/(-scale) # Quantile function of generalized exponential distribution applied to halton draws
  mu_i <- lambda * gedist

  p_pge.i <- sapply(mu_i, stats::dpois, x=x)

  p <- mean(p_pge.i)

  if (log) return(log(p))
  else return(p)
})

#' @rdname PoissonGeneralizedExponential
#' @export
ppge <- function(q, mean = 1, shape = 1, scale = 1, ndraws = 1500,
                 lower.tail = TRUE, log.p = FALSE, haltons = NULL) {
  
  
  # --- Input Validation ---
  
  if (any(mean <= 0, na.rm = TRUE)) warning("'mean' must be positive")
  if (any(shape <= 0, na.rm = TRUE)) warning("'shape' must be positive")
  if (any(scale <= 0, na.rm = TRUE))  warning("'scale' must be positive")
  
  
  # --- Vectorization Setup ---
  
  n <- max(length(q), length(mean), length(shape), length(scale))
  q <- rep_len(as.integer(floor(q)), n)
  mean <- rep_len(mean, n)
  shape <- rep_len(shape, n)
  scale <- rep_len(scale, n)
  
  
  # --- Initialize Result ---
  
  cdf <- rep(NA_real_, n)
  
  
  # --- Handle Special Cases ---
  
  invalid_q <- is.na(q) | q < 0
  cdf[invalid_q] <- 0
  
  
  # --- Pre-compute Halton draws (shared across all evaluations) ---
  
  if (is.null(haltons)) {
    haltons <- randtoolbox::halton(ndraws)
  }
  
  
  # --- Compute CDF for Valid Cases ---
  
  valid <- !invalid_q
  
  if (any(valid)) {
    # Group by unique parameter combinations
    param_df <- data.frame(
      idx = which(valid),
      q = q[valid],
      mean = mean[valid],
      shape = shape[valid],
      scale = scale[valid]
    )
    
    unique_params <- unique(param_df[, c("mean", "shape", "scale")])
    
    for (j in seq_len(nrow(unique_params))) {
      m_j <- unique_params$mean[j]
      sh_j <- unique_params$shape[j]
      sc_j <- unique_params$scale[j]
      
      # Find all indices with these parameters
      mask <- param_df$mean == m_j & param_df$shape == sh_j & param_df$scale == sc_j
      indices <- param_df$idx[mask]
      q_vals <- param_df$q[mask]
      max_q <- max(q_vals)
      
      # Compute lambda for this parameter set
      lambda_j <- m_j * sc_j / (digamma(sh_j + 1) - digamma(1))
      
      # Generate GE quantiles from Halton draws
      # Quantile function of GE: Q(p) = -log(1 - p^(1/shape)) / scale
      ge_quantiles <- -log(1 - haltons^(1 / sh_j)) / sc_j
      
      # Compute mu values for each Halton draw
      mu_draws <- lambda_j * ge_quantiles
      
      # Compute PMF for 0:max_q by averaging over Halton draws
      pmf <- numeric(max_q + 1)
      
      for (y in 0:max_q) {
        # Average Poisson PMF over all draws
        pmf[y + 1] <- mean(stats::dpois(y, mu_draws))
      }
      
      # Compute cumulative sums
      cum_pmf <- cumsum(pmf)
      
      # Assign CDF values
      for (k in seq_along(indices)) {
        idx_k <- indices[k]
        q_k <- q_vals[k]
        cdf[idx_k] <- cum_pmf[q_k + 1]
      }
    }
  }
  
  
  # --- Apply Tail and Log Options ---
  
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  return(cdf)
}
#' @rdname PoissonGeneralizedExponential
#' @export
qpge <- Vectorize(function(p, mean=1, shape=1, scale=1, ndraws=1500) {
  if(mean<=0 || scale<=0 || shape <=0){
    warning('The values of `mean`, `shape`, and `scale` have to have values greater than 0.')
  }

  y <- 0
  p_value <- ppge(y, mean, shape, scale, ndraws=ndraws)
  while(p_value < p){
    y <- y + 1
    p_value <- ppge(y, mean, shape, scale, ndraws=ndraws)
  }
  return(y)
})


#' @rdname PoissonGeneralizedExponential
#' @export
rpge <- function(n, mean=1, shape=1, scale=1, ndraws=1500) {
  if(mean<=0 || scale<=0 || shape <=0){
    warning('The values of `mean`, `shape`, and `scale` have to have values greater than 0.')
  }

  u <- runif(n)
  y <- sapply(u, function(p) qpge(p, mean, shape, scale, ndraws))
  return(y)
}
