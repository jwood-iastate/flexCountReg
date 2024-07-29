#' Poisson-Lindley-Lognormal Distribution
#'
#' These functions provide density, distribution function, quantile function, and random number generation for the Poisson-Lindley-Lognormal (PLL) Distribution
#'
#' The Poisson-Lindley-Lognormal is a 3-parameter count distribution that captures high densities for small integer values and provides flexibility for heavy tails.
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mean numeric value or vector of mean values for the distribution (the values have to be greater than 0).
#' @param theta single value or vector of values for the theta parameter of the distribution (the values have to be greater than 0).
#' @param sigma single value or vector of values for the sigma parameter of the lognormal distribution (the values have to be greater than 0).
#' @param lambda alternative parameterization (use instead of the mean); numeric value or vector of values for lambda parameter of the distribution (the values have to be greater than 0).
#' @param ndraws the number of Halton draws to use for the integration.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dplind} computes the density (PDF) of the Poisson-Lindley Distribution.
#'
#' \code{pplind} computes the CDF of the Poisson-Lindley Distribution.
#'
#' \code{qplind} computes the quantile function of the Poisson-Lindley Distribution.
#'
#' \code{rplind} generates random numbers from the Poisson-Lindley Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Poisson-Lindley-Lognormal (PLL) distribution is:
#' \deqn{f(y|\mu,\theta,\alpha)=\int_0^\infty \frac{\theta^2\mu^y x^y(\theta+\mu\cdot x+y+1)}{(\theta+1)(\theta+\mu\cdot x)^{y+2}} \frac{\exp\left(-\frac{ln^2(x)}{2\sigma^2} \right)}{x\sigma\sqrt{2\pi}}dx}
#'
#' Where \eqn{\theta} and \eqn{\lambda} are distribution parameters from the Poisson-Lindley distribution with the restrictions that \eqn{\theta>0} and \eqn{\lambda>0}, \eqn{\sigma} is a parameter for the lognormal distribution with the restriction \eqn{\sigma>0}, and \eqn{y} is a non-negative integer.
#'
#' The expected value of the distribution is:
#' \deqn{E[y]=\mu=e^{X\beta}=\frac{\lambda(\theta+2)e^{\frac{\sigma^2}{2}}}{\theta(\theta+1)}}
#'
#' The default is to use the input mean value for the distribution. However, the lambda parameter can be used as an alternative to the mean value.
#'
#' Halton draws are used to perform simulation over the lognormal distribution to solve the integral.
#'
#' @examples
#' dplindLnorm(0, mean=0.75, theta=7, sigma=2, ndraws=10)
#' pplindLnorm(c(0,1,2,3,5,7,9,10), mean=0.75, theta=7, sigma=2, ndraws=10)
#' qplindLnorm(c(0.1,0.3,0.5,0.9,0.95), lambda=4.67, theta=7, sigma=2, ndraws=10)
#' rplindLnorm(3, mean=0.75, theta=7, sigma=2, ndraws=10)
#'
#' @import stats randtoolbox
#' @include plind.R
#' @export
#' @name Poisson-Lindley-Lognormal

#' @rdname Poisson-Lindley-Lognormal
#' @export
dplindLnorm <- Vectorize(function(x, mean=1, theta = 1, sigma=1, lambda=NULL, ndraws=1500, log=FALSE){
  #test to make sure the value of x is an integer
  tst <- ifelse(is.na(nchar(strsplit(as.character(x), "\\.")[[1]][2])>0),FALSE, TRUE)
  if(tst || x < 0){
    print("The value of `x` must be a non-negative whole number")
    stop()
  }
  if(is.null(lambda)){
    if(mean<=0 || theta<=0 || sigma<=0){
      print('The values of `mean`, `theta`, and `sigma` all have to have values greater than 0.')
      stop()
    }
    else{
      lambda <- mean*theta*(theta+1)/((theta+2) * exp(sigma^2/2))
    }
  }
  else{
    if(lambda<=0 || theta<=0  || sigma<=0){
      print('The values of `lambda`, `theta`, and `sigma` all have to have values greater than 0.')
      stop()
    }
  }
  
  # Generate Halton draws to use as quantile values
  h <- randtoolbox::halton(ndraws)
  
  # Evaluate the density of the normal distribution at those quantiles and use the exponent to transform to lognormal values
  lnormdist <- exp(stats::qnorm(h, 0, sigma))
  
  mu <- lambda*(theta+2)/(theta*(theta+1))
  mu_i <- outer(mu, lnormdist)
  
  p_plind.i <- sapply(mu_i, function(y) dplind(x=x, mean=y, theta=theta))
  
  p <- mean(p_plind.i)
  
  if (log) return(log(p))
  else return(p)
})

#' @rdname Poisson-Lindley-Lognormal
#' @export
pplindLnorm <- Vectorize(function(q, mean=1, theta = 1, lambda=NULL, sigma=1, ndraws=1500, lower.tail=TRUE, log.p=FALSE){
  if(is.null(lambda)){
    if(mean<=0 || theta<=0  || sigma<=0){
      print('The values of `mean`, `theta`, and `sigma` all have to have values greater than 0.')
      stop()
    }
  }
  else{
    if(lambda<=0 || theta<=0  || sigma<=0){
      print('The values of `lambda`, `theta`, and `sigma` all have to have values greater than 0.')
      stop()
    }
    else{
      mean <- lambda*(theta+2)/(theta*(theta+1))*exp(sigma^2/2)
    }
  }
  
  y <- seq(0,q,1)
  probs <- dplindLnorm(y, mean, theta, sigma=sigma, ndraws=ndraws)
  p <- sum(probs)
  
  if(!lower.tail) p <- 1-p
  
  if (log.p) return(log(p))
  else return(p)
})

#' @rdname Poisson-Lindley-Lognormal
#' @export
qplindLnorm <- Vectorize(function(p, mean=1, theta=1, sigma=1, ndraws=1500, lambda=NULL) {
  if(is.null(lambda)){
    if(mean<=0 || theta<=0  || sigma<=0){
      print('The values of `mean`, `theta`, and `sigma` all have to have values greater than 0.')
      stop()
    }
  }
  else{
    if(lambda<=0 || theta<=0  || sigma<=0){
      print('The values of `lambda`, `theta`, and `sigma` all have to have values greater than 0.')
      stop()
    }
    else{
      mean <- lambda*(theta+2)/(theta*(theta+1))*exp(sigma^2/2)
    }
  }
  
  y <- 0
  p_value <- pplindLnorm(y, mean, theta, sigma=sigma, ndraws=ndraws)
  while(p_value < p){
    y <- y + 1
    p_value <- pplindLnorm(y, mean, theta, sigma=sigma, ndraws=ndraws)
  }
  return(y)
})


#' @rdname Poisson-Lindley-Lognormal
#' @export
rplindLnorm <- function(n, mean=1, theta=1, sigma=1, ndraws=1500, lambda=NULL) {
  if(is.null(lambda)){
    if(mean<=0 || theta<=0  || sigma<=0){
      print('The values of `mean`, `theta`, and `sigma` allboth have to have values greater than 0.')
      stop()
    }
  }
  else{
    if(lambda<=0 || theta<=0  || sigma<=0){
      print('The values of `lambda`, `theta`, and `sigma` all have to have values greater than 0.')
      stop()
    }
    else{
      mean <- lambda*(theta+2)/(theta*(theta+1))*exp(sigma^2/2)
    }
  }
  
  u <- runif(n)
  y <- sapply(u, function(p) qplindLnorm(p, mean, theta, sigma=sigma, ndraws=ndraws))
  return(y)
}