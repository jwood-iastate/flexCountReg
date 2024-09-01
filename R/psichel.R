#' Sichel Distribution
#'
#' These functions provide the density function, distribution function, quantile function, and random number generation for the Sichel (SI) Distribution
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mu numeric value or vector of mean values for the distribution (the values have to be greater than 0).
#' @param sigma single value or vector of values for the scale parameter of the distribution (the values have to be greater than 0).
#' @param gamma single value or vector of values for the shape parameter of the distribution.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dsichel} computes the density (PDF) of the Sichel Distribution.
#'
#' \code{psichel} computes the CDF of the Sichel Distribution.
#'
#' \code{qsichel} computes the quantile function of the Sichel Distribution.
#'
#' \code{rsichel} generates random numbers from the Sichel Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Sichel distribution uses the formulation from Zhou et al. (2011) and Rigby et al. (2008):
#' \deqn{f(y|\mu, \sigma, \gamma)=\frac{\left(\frac{\mu}{c}\right)^y K_{y+\gamma}(\alpha)}{K_\gamma(1/\sigma)y!(\alpha\sigma)^{y+\gamma}}}
#' 
#' Where \eqn{\sigma} and \eqn{\gamma} are distribution parameters with \eqn{-\infty < \gamma < \infty} and \eqn{\sigma>0}, 
#' \eqn{c=\frac{K_{\gamma+1}(1/\sigma)}{K_\gamma(1/\sigma)}}, \eqn{\alpha^2=\sigma^{-2}+2\mu(c\sigma)^{-1}}, 
#' a mean value of \eqn{\mu}, \eqn{y} is a non-negative integer, and \eqn{K_j(x)} is a modified Bessel function of the third kind with order \eqn{j} and argument \eqn{x}.
#'
#' The variance of the distribution is:
#' \deqn{\sigma^2=\mu+\left(\frac{2\sigma(\gamma+1)}{c}+\frac{1}{c^2}-1\right)\mu^2}
#' 
#' @references 
#' Zou, Y., Lord, D., & Zhang, Y. (2012). Analyzing highly dispersed crash data using the Sichel generalized additive models for location, scale and shape. Paper submitted for publication.
#' 
#' Rigby, R. A., Stasinopoulos, D. M., & Akantziliotou, C. (2008). A framework for modelling overdispersed count data, including the Poisson-shifted generalized inverse Gaussian distribution. Computational Statistics & Data Analysis, 53(2), 381-393.
#' 
#' @examples
#' dsichel(1, mu=0.75, sigma=1, gamma=-3)
#' psichel(c(0,1,2,3,5,7,9,10), mu=0.75, sigma=2, gamma=3)
#' qsichel(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, sigma=1, gamma=15)
#' rsichel(30, mu=0.75, sigma=0.5, gamma=1)
#'
#' @import stats
#' @export
#' @name SichelDistribution

#' @rdname SichelDistribution
#' @export
dsichel <- Vectorize(function(x, mu=1, sigma = 1, gamma=1, log=FALSE){
  #test to make sure the value of x is an integer
  tst <- ifelse(is.na(nchar(strsplit(as.character(x), "\\.")[[1]][2])>0),FALSE, TRUE)
  if(tst || x < 0){
    print("The value of `x` must be a non-negative whole number")
    stop()
  }
  if(sigma<=0){
    print("The value of `sigma` must be greater than 0.")
    stop()
  }

  K_g1_invsig <- besselK(x=(1/sigma), nu=(gamma+1)) 
  K_g_invsig <-  besselK(x=(1/sigma), nu=(gamma)) 
  c <- K_g1_invsig/K_g_invsig
  a <- sqrt(sigma^(-2)+2*mu*(c*sigma)^(-1))
  K_yg_a <- besselK(x=a, nu=(x+gamma))
  p <- (mu/c)^x * K_yg_a/(K_g_invsig * factorial(x) * (a*sigma)^(x+gamma))
  
  if (log) return(log(p))
  else return(p)
})

#' @rdname SichelDistribution
#' @export
psichel <- Vectorize(function(q, mu=1, sigma = 1, gamma=1, lower.tail=TRUE, log.p=FALSE){
  
  y <- seq(0,q,1)
  probs <- dsichel(y, mu, sigma, gamma)
  p <- sum(probs)

  if(!lower.tail) p <- 1-p

  if (log.p) return(log(p))
  else return(p)
})

#' @rdname SichelDistribution
#' @export
qsichel <- Vectorize(function(p, mu=1, sigma = 1, gamma=1) {
  y <- 0
  p_value <- psichel(y, mu, sigma, gamma)
  if(is.na(p_value) || is.null(p_value)){
    print('Error computing probability')
    stop()
  }
  while(p_value < p){
    y <- y + 1
    p_value <- psichel(y, mu, sigma, gamma)
  }
  return(y)
})

#' @rdname SichelDistribution
#' @export
rsichel <- function(n, mu=1, sigma = 1, gamma=1) {
  u <- runif(n)
  y <- sapply(u, qsichel, mu=mu, sigma=sigma, gamma=gamma)
  return(y)
}
