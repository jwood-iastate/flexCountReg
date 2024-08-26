#' Poisson-Lognormal Distribution
#'
#' These functions provide density, distribution function, quantile function, and random number generation for the Poisson-Lognormal (PLogN) Distribution
#'
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mean numeric value or vector of mean values for the distribution (the values have to be greater than 0).
#' @param sigma single value or vector of values for the sigma parameter of the lognormal distribution (the values have to be greater than 0).
#' @param ndraws the number of Halton draws to use for the integration.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dpLnorm} computes the density (PDF) of the Poisson-Lognormal Distribution.
#'
#' \code{ppLnorm} computes the CDF of the Poisson-Lognormal Distribution.
#'
#' \code{qpLnorm} computes the quantile function of the Poisson-Lognormal Distribution.
#'
#' \code{rpLnorm} generates random numbers from the Poisson-Lognormal Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Poisson-Lognormal distribution is:
#' \deqn{f(y|\mu,\theta,\alpha)=\int_0^\infty \frac{\mu^y x^y e^{-\mu x}}{y!}\frac{exp\left(-\frac{ln^2(x)}{2\sigma^2} \right)}{x\sigma\sqrt{2\pi}}dx}
#'
#' Where \eqn{\sigma} is a parameter for the lognormal distribution with the restriction \eqn{\sigma>0}, and \eqn{y} is a non-negative integer.
#'
#' The expected value of the distribution is:
#' \deqn{E[y]=e^{X\beta+\sigma^2/2} = \mu e^{\sigma^2/2}}
#' Halton draws are used to perform simulation over the lognormal distribution to solve the integral.
#'
#' @examples
#' dpLnorm(0, mean=0.75, sigma=2, ndraws=10)
#' ppLnorm(c(0,1,2,3,5,7,9,10), mean=0.75, sigma=2, ndraws=10)
#' qpLnorm(c(0.1,0.3,0.5,0.9,0.95), mean=0.75, sigma=2, ndraws=10)
#' rpLnorm(30, mean=0.75,  sigma=2, ndraws=10)
#'
#' @import stats randtoolbox
#' @export
#' @name Poisson-Lognormal
#' 
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @rdname Poisson-Lognormal
#' @export
dpLnorm_cpp <- Vectorize(function(x, mean=1, sigma=1, ndraws=1500, log=FALSE){
  
  if(mean<=0 || sigma<=0){
    print('The values of `mean` and `sigma` have to have values greater than 0.')
    stop()
  }
  
  # Generate Halton draws to use as quantile values
  h <- randtoolbox::halton(ndraws)
  
  p <- calculate_plogn_prob(x, mean, sigma, h)
  
  if (log) return(log(p))
  else return(p)
})

#' @rdname Poisson-Lognormal
#' @export
dpLnorm <- Vectorize(function(x, mean=1, sigma=1, ndraws=1500, log=FALSE){
  
  if(mean<=0 || sigma<=0){
    print('The values of `mean` and `sigma` have to have values greater than 0.')
    stop()
  }
  
  # Generate Halton draws to use as quantile values
  h <- randtoolbox::halton(ndraws)
  
  # Evaluate the density of the normal distribution at those quantiles and use the exponent to transform to lognormal values
  lnormdist <- exp(stats::qnorm(h, 0, sigma))
  mu_i <- outer(mean, lnormdist)
  
  p_plogn.i <- sapply(mu_i, stats::dpois, x=x)
  
  p <- mean(p_plogn.i)
  
  if (log) return(log(p))
  else return(p)
})

#' @rdname Poisson-Lognormal
#' @export
ppLnorm <- Vectorize(function(q, mean=1, sigma=1, ndraws=1500, lower.tail=TRUE, log.p=FALSE){
  if(mean<=0 || sigma<=0){
    print('The values of `mean` and `sigma` have to have values greater than 0.')
    stop()
  }
  
  y <- seq(0,q,1)
  probs <- dpLnorm(y, mean, sigma=sigma, ndraws=ndraws)
  p <- sum(probs)
  
  if(!lower.tail) p <- 1-p
  
  if (log.p) return(log(p))
  else return(p)
})

#' @rdname Poisson-Lognormal
#' @export
qpLnorm <- Vectorize(function(p, mean=1, sigma=1, ndraws=1500) {
  if(mean<=0 || sigma<=0){
    print('The values of `mean`  and `sigma` have to have values greater than 0.')
    stop()
  }
  
  y <- 0
  p_value <- ppLnorm(y, mean, sigma=sigma, ndraws=ndraws)
  while(p_value < p){
    y <- y + 1
    p_value <- ppLnorm(y, mean, sigma=sigma, ndraws=ndraws)
  }
  return(y)
})


#' @rdname Poisson-Lognormal
#' @export
rpLnorm <- function(n, mean=1, sigma=1, ndraws=1500) {
  if(mean<=0  || sigma<=0){
    print('The values of `mean` and `sigma` have to have values greater than 0.')
    stop()
  }
  
  u <- stats::runif(n)
  y <- sapply(u, function(p) qpLnorm(p, mean, sigma=sigma, ndraws=ndraws))
  return(y)
}