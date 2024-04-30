#' Poisson-Generalized-Exponential Distribution
#'
#' These functions provide density, distribution function, quantile function, and random number generation for the Poisson-Generalized-Exponential (PGE) Distribution
#'
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mean numeric value or vector of mean values for the distribution (the values have to be greater than 0). This is NOT the value of \eqn{\lambda}.
#' @param shape numeric value or vector of shape values for the shape parameter of the generalized exponential distribution (the values have to be greater than 0).
#' @param scale single value or vector of values for the scale parameter of the generalized exponential distribution (the values have to be greater than 0).
#' @param ndraws the number of Halton draws to use for the integration.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
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
#' The Generalized Exponential distribution can be written as a function with a shape parameter \eqn{\alpha>0} and scale parameter \eqn{\gamma>0}. The distribution has striclty positive continuous values. The PDF of the distribution is:
#' \deqn{f(x|\alph,\gamma)=\frac{\alpha}{\gamma}\left(1-e^{-\frac{x}{\gamma}}\right)^{\alpha-1}e^{-\frac{x}{\gamma}}} 
#' 
#' Thus, the compound Probability Mass Function(PMF) for the PGE distribution is:
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
#' @import stats randtoolbox
#' @export
#' @name Poisson-Generalized-Exponential

#' @rdname Poisson-Generalized-Exponential
#' @export
dpge <- Vectorize(function(x, mean=1, shape=1, scale=1, ndraws=1500, log=FALSE){
  
  if(mean<=0 || scale<=0 || shape <=0){
    print('The values of `mean`, `shape`, and `scale` have to have values greater than 0.')
    stop()
  }
  
  qge <- function(p, shape, scale){
    q <- log((p*shape+1)^(1/shape)+1)/scale
    return(q)
  }
  
  lambda <- mean * scale /(digamma(shape+1)-digamma(1))

  # Generate Halton draws to use as quantile values
  h <- randtoolbox::halton(ndraws)

  # Evaluate the density of the normal distribution at those quantiles and use the exponent to transform to lognormal values
  gedist <- qge(h, shape, scale)
  mu_i <- lambda * gedist

  p_pge.i <- sapply(mu_i, stats::dpois, x=x)

  p <- mean(p_pge.i)

  if (log) return(log(p))
  else return(p)
})

#' @rdname Poisson-Generalized-Exponential
#' @export
ppge <- Vectorize(function(q, mean=1, shape=1, scale=1, ndraws=1500, lower.tail=TRUE, log.p=FALSE){
  if(mean<=0 || scale<=0 || shape <=0){
    print('The values of `mean`, `shape`, and `scale` have to have values greater than 0.')
    stop()
  }

  y <- seq(0,q,1)
  probs <- dpge(y, mean, shape, scale, ndraws=ndraws)
  p <- sum(probs)

  if(!lower.tail) p <- 1-p

  if (log.p) return(log(p))
  else return(p)
})

#' @rdname Poisson-Generalized-Exponential
#' @export
qpge <- Vectorize(function(p, mean=1, shape=1, scale=1, ndraws=1500) {
  if(mean<=0 || scale<=0 || shape <=0){
    print('The values of `mean`, `shape`, and `scale` have to have values greater than 0.')
    stop()
  }

  y <- 0
  p_value <- ppge(y, mean, shape, scale, ndraws=ndraws)
  while(p_value < p){
    y <- y + 1
    p_value <- ppge(y, mean, shape, scale, ndraws=ndraws)
  }
  return(y)
})


#' @rdname Poisson-Generalized-Exponential
#' @export
rpge <- function(n, mean=1, shape=1, scale=1, ndraws=1500) {
  if(mean<=0 || scale<=0 || shape <=0){
    print('The values of `mean`, `shape`, and `scale` have to have values greater than 0.')
    stop()
  }

  u <- runif(n)
  y <- sapply(u, function(p) qpge(p, mean, shape, scale, ndraws))
  return(y)
}
