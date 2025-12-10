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
#' @param mean mean (>0).
#' @param theta Poisson-Lindley theta parameter (>0).
#' @param sigma lognormal sigma parameter (>0).
#' @param lambda optional lambda parameter (>0).
#' @param ndraws number of Halton draws.
#' @param log return log-density.
#' @param log.p return log-CDF.
#' @param lower.tail TRUE returns P[X ≤ x].
#' @param hdraws optional Halton draws.
#'
#' @details
#' \code{dplind} computes the PLL density.
#' \code{pplind} computes the PLL CDF.
#' \code{qplind} computes quantiles.
#' \code{rplind} generates random draws.
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
#' Mean:
#' \deqn{
#' E[y]=\mu=\frac{\lambda(\theta+2)e^{\sigma^2/2}}
#'             {\theta(\theta+1)}
#' }
#'
#' Halton draws are used to evaluate the integral.
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
    mean <- lambda * (theta + 2) /
      (theta * (theta + 1)) * exp(sigma^2 / 2)
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
qplindLnorm <- Vectorize(function(p, mean=1, theta=1,
                                  sigma=1, ndraws=1500,
                                  lambda=NULL) {
  
  if(is.null(lambda)) {
    if(mean <= 0 || theta <= 0 || sigma <= 0) {
      warning(paste(
        "The values of `mean`, `theta`, and `sigma` must",
        "all be greater than 0."
      ))
    }
  } else {
    if(lambda <= 0 || theta <= 0 || sigma <= 0) {
      warning(paste(
        "The values of `lambda`, `theta`, and `sigma` must",
        "all be greater than 0."
      ))
    } else {
      mean <- lambda * (theta + 2) /
        (theta * (theta + 1)) * exp(sigma^2 / 2)
    }
  }
  
  y <- 0
  p_value <- pplindLnorm(y, mean, theta,
                         sigma = sigma,
                         ndraws = ndraws)
  
  while(p_value < p) {
    y <- y + 1
    p_value <- pplindLnorm(
      y, mean, theta, sigma = sigma,
      ndraws = ndraws
    )
  }
  
  return(y)
})

#' @rdname PoissonLindleyLognormal
#' @export
rplindLnorm <- function(n, mean=1, theta=1, sigma=1,
                        ndraws=1500, lambda=NULL) {
  
  if(is.null(lambda)) {
    if(mean <= 0 || theta <= 0 || sigma <= 0) {
      warning(paste(
        "The values of `mean`, `theta`, and `sigma` must",
        "all be greater than 0."
      ))
    }
  } else {
    if(lambda <= 0 || theta <= 0 || sigma <= 0) {
      warning(paste(
        "The values of `lambda`, `theta`, and `sigma` must",
        "all be greater than 0."
      ))
    } else {
      mean <- lambda * (theta + 2) /
        (theta * (theta + 1)) * exp(sigma^2 / 2)
    }
  }
  
  u <- runif(n)
  y <- sapply(
    u,
    function(p) qplindLnorm(
      p, mean, theta, sigma = sigma,
      ndraws = ndraws
    )
  )
  
  return(y)
}
