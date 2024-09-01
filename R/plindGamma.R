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
#' @param ndraws the number of Halton draws to use for the integration.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
#' @param hdraws and optional vector of Halton draws to use for the integration.
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
#' \deqn{f(x|\mu,\theta,\alpha)\frac{\alpha  (\theta+2) ^2 \Gamma (x+\alpha ) }{(\mu)^2(\theta +1)^3 \Gamma (\alpha )}\left(\frac{\mu\theta(\theta+1)}{\theta+2} U\left(x+1,2-\alpha ,\frac{\alpha  (\theta+2) }{\mu(\theta+1)}\right)+\alpha(x+1) U\left(x+2,3-\alpha ,\frac{\alpha  (\theta+2) }{\mu(\theta+1)}\right)\right)}
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
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @name NegativeBinomialLindley
#' 
#' @rdname NegativeBinomialLindley
#' @export
dplindGamma <- function(x, mean=1, theta = 1, alpha=1, log=FALSE, ndraws=1000, hdraws=NULL){
  if (is.null(hdraws)){
    h <- randtoolbox::halton(ndraws)
  }
  
  p <- dplindgamma_cpp(x, mean, theta, alpha, h)

  if (log) return(log(p))
  else return(p)
}

#' @rdname NegativeBinomialLindley
#' @export
pplindGamma <- Vectorize(function(q, mean=1, theta = 1, alpha=1, lower.tail=TRUE, log.p=FALSE){
  if(mean<=0 || theta<=0  || alpha<=0){
    print('The values of `mean`, `theta`, and `alpha` all have to have values greater than 0.')
    stop()
  }
  
  y <- seq(0,q,1)
  probs <- dplindGamma(y, mean, theta, alpha)
  p <- sum(probs)
  
  if(!lower.tail) p <- 1-p
  
  if (log.p) return(log(p))
  else return(p)
})

#' @rdname NegativeBinomialLindley
#' @export
qplindGamma <- Vectorize(function(p, mean=1, theta=1, alpha=1) {
  if(p < 0){
    print("The value of `p` must be a value greater than 0 and less than 1.")
    stop()
  }
  if(is.na(p)){
    print("The value of `p` cannot be an `NA` value")
    stop()
  }

  if(mean<=0 || theta<=0 || alpha<=0){
    print('The values of `mean`, `theta`, and `alpha` all have to have values greater than 0.')
    stop()
  }
 
  
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

  if(mean<=0 || theta<=0  || alpha<=0){
    print('The values of `mean`, `theta`, and `alpha` all have to have values greater than 0.')
    stop()
  }
  
  u <- runif(n)
  y <- lapply(u, function(p) qplindGamma(p, mean, theta, alpha=alpha))
  return(unlist(y))
}
