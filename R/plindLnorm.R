#' Poisson-Lindley-Lognormal Distribution
#'
#' These functions provide density, distribution, quantile, and random
#' generation for the Poisson-Lindley-Lognormal (PLL) Distribution.
#'
#' The PLL is a 3-parameter count distribution that captures high mass at
#' small y and allows flexible heavy tails.
#'
#' @param x numeric value or vector of values.
#' @param q quantile or vector of quantiles.
#' @param p probability or vector of probabilities.
#' @param n number of random draws.
#' @param mean the Poisson-Lindley mean before the lognormal factor (>0); the
#'   mean of the distribution is \code{mean * exp(sigma^2 / 2)}.
#' @param theta Poisson-Lindley theta parameter (>0).
#' @param sigma lognormal sigma parameter (>0).
#' @param lambda optional Poisson scale parameter (>0), an alternative to
#'   \code{mean}: \code{mean = lambda * (theta + 2) / (theta * (theta + 1))}.
#' @param ndraws number of Halton draws.
#' @param log return log-density.
#' @param log.p return log-CDF.
#' @param lower.tail TRUE returns \eqn{P[X \leq x]}.
#' @param hdraws optional Halton draws.
#'
#' @details
#' \code{dplindLnorm} computes the PLL density.
#' \code{pplindLnorm} computes the PLL CDF.
#' \code{qplindLnorm} computes quantiles.
#' \code{rplindLnorm} generates random draws.
#'
#' The PMF is:
#' \deqn{
#' f(y|\mu,\theta,\sigma)=\int_0^\infty
#'   \frac{\theta^2\mu^y x^y(\theta+\mu x+y+1)}
#'        {(\theta+1)(\theta+\mu x)^{y+2}}
#'   \frac{\exp\left(-\frac{\ln^2(x)}{2\sigma^2}\right)}
#'        {x\sigma\sqrt{2\pi}}dx
#' }
#'
#' Where \eqn{\mu} is the \code{mean} argument (the Poisson-Lindley mean). The
#' mean of the distribution is:
#' \deqn{
#' E[y]=\mu e^{\sigma^2/2}=\frac{\lambda(\theta+2)e^{\sigma^2/2}}
#'             {\theta(\theta+1)}
#' }
#'
#' Halton draws are used to evaluate the integral.
#' 
#' @returns dplindLnorm gives the density, pplindLnorm gives the distribution 
#'  function, qplindLnorm gives the quantile function, and rplindLnorm generates
#'  random  deviates.
#' 
#'  The length of the result is determined by n for rplindLnorm, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' dplindLnorm(0, mean=0.75, theta=7, sigma=2, ndraws=10)
#' pplindLnorm(0:10, mean=0.75, theta=7, sigma=2, ndraws=10)
#' qplindLnorm(c(0.1,0.5,0.9), lambda=4.67, theta=7, sigma=2)
#' rplindLnorm(5, mean=0.75, theta=7, sigma=2)
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @name PoissonLindleyLognormal
#'
#' @rdname PoissonLindleyLognormal
#' @export
dplindLnorm <- function(x, mean=1, theta=1, sigma=1,
                        ndraws=1500, log=FALSE, hdraws=NULL) {
  
  if(!is.null(hdraws)) {
    h <- qnorm(hdraws)
  } else {
    h <- randtoolbox::halton(ndraws, normal=TRUE)
  }
  
  p <- dplindlogn_cpp(x, mean, theta, sigma, h)
  
  if(log) return(log(p))
  return(p)
}

#' @rdname PoissonLindleyLognormal
#' @export
pplindLnorm <- function(q, mean=1, theta=1, lambda=NULL,
                        sigma=1, ndraws=1500,
                        lower.tail=TRUE, log.p=FALSE) {
  
  if(!is.null(lambda)) {
    mean <- lambda * (theta + 2) / (theta * (theta + 1))
  }
  
  # FIX: Set normal=FALSE.
  # dplindLnorm applies qnorm() internally.
  h <- randtoolbox::halton(ndraws, normal=FALSE)
  
  n <- max(length(q), length(mean), length(theta), length(sigma))
  q <- rep_len(q, n)
  mean <- rep_len(mean, n)
  theta <- rep_len(theta, n)
  sigma <- rep_len(sigma, n)
  
  cdf <- numeric(n)
  
  for(i in seq_len(n)) {
    if(q[i] < 0) {
      cdf[i] <- 0
      next
    }
    y_seq <- 0:floor(q[i])
    probs <- dplindLnorm(
      y_seq, mean[i], theta[i], sigma[i],
      ndraws = ndraws, hdraws = h
    )
    cdf[i] <- sum(probs)
  }
  
  if(!lower.tail) cdf <- 1 - cdf
  if(log.p) return(log(cdf))
  return(cdf)
}

#' @rdname PoissonLindleyLognormal
#' @export
qplindLnorm <- function(p, mean = 1, theta = 1, sigma = 1, ndraws = 1500,
                        lambda = NULL) {
  if (!is.null(lambda)) {
    mean <- lambda * (theta + 2) / (theta * (theta + 1))
  }
  if (any(mean <= 0 | theta <= 0 | sigma <= 0, na.rm = TRUE)) {
    warning(paste("The values of `mean` (or `lambda`), `theta`, and",
                  "`sigma` must all be greater than 0."))
  }
  h <- randtoolbox::halton(ndraws)
  .fcr_qcount(p, function(y, mean, theta, sigma) {
    cumsum(dplindLnorm(y, mean, theta, sigma, hdraws = h))
  }, mean, theta, sigma)
}

#' @rdname PoissonLindleyLognormal
#' @export
rplindLnorm <- function(n, mean = 1, theta = 1, sigma = 1,
                        ndraws = 1500, lambda = NULL) {
  n <- .fcr_rn(n)
  theta <- rep_len(theta, n); sigma <- rep_len(sigma, n)
  if (is.null(lambda)) {
    mean <- rep_len(mean, n)
    lambda <- mean * theta * (theta + 1) / (theta + 2)
  } else {
    lambda <- rep_len(lambda, n)
  }
  bad <- !(is.finite(lambda) & lambda > 0 & is.finite(theta) & theta > 0 &
             is.finite(sigma) & sigma > 0)
  .fcr_warn_na(bad, "`mean` (or `lambda`), `theta` and `sigma` must be > 0.")
  out <- rep(NA_real_, n); ok <- !bad; k <- sum(ok)
  out[ok] <- stats::rpois(k, lambda[ok] * .fcr_rlindley_mix(k, theta[ok]) *
                            exp(stats::rnorm(k, 0, sigma[ok])))
  out
}
