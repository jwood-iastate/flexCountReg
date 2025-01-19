#' Poisson-Inverse-Gamma Distribution
#'
#' These functions provide the density function, distribution function, quantile 
#' function, and random number generation for the Poisson-Inverse-Gamma 
#' (PInvGamma) Distribution

#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mu numeric value or vector of mean values for the distribution (the 
#' values have to be greater than 0).
#' @param eta single value or vector of values for the scale parameter of the 
#' distribution (the values have to be greater than 0).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dpinvgamma} computes the density (PDF) of the Poisson-Inverse-Gamma Distribution.
#'
#' \code{ppinvgamma} computes the CDF of the Poisson-Inverse-Gama Distribution.
#'
#' \code{qpinvgamma} computes the quantile function of the Poisson-Inverse-Gamma Distribution.
#'
#' \code{rpinvgamma} generates random numbers from the Poisson-Inverse-Gamma Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Poisson-Inverse-Gamma 
#' distribution is:
#' \deqn{f(x|\eta,\mu)=\frac{2\left(\mu\left(\frac{1}{\eta}+1\right)\right)^{\frac{x+\frac{1}{eta}+2}{2}}}{x!\Gamma\left(\frac{1}{\eta}+2\right)}K_{x-\frac{1}{\eta}-2}\left(2\sqrt{\mu\left(\frac{1}{\eta}+1\right)}\right)}
#' 
#' Where \eqn{\eta} is a shape parameter with the restriction that \eqn{\eta>0}, 
#' \eqn{\mu>0} is the mean value,  \eqn{y} is a non-negative integer, and 
#' \eqn{K_i(z)} is the modified Bessel function of the second kind. This 
#' formulation uses the mean directly.
#'
#' The variance of the distribution is:
#' \deqn{\sigma^2=\mu+\eta\mu^2}
#' 
#' 
#' @examples
#' dpinvgamma(1, mu=0.75, eta=1)
#' ppinvgamma(c(0,1,2,3,5,7,9,10), mu=0.75, eta=3)
#' qpinvgamma(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, eta=0.5)
#' rpinvgamma(30, mu=0.75, eta=1.5)
#'
#' @importFrom stats runif
#' @export
#' @name PoissonInverseGamma
#' @rdname PoissonInverseGamma
#' @export
dpinvgamma <- Vectorize(function(x, mu=1, eta = 1,  log=FALSE){
  #test to make sure the value of x is an integer
  tst <- ifelse(is.na(nchar(strsplit(as.character(x), "\\.")[[1]][2])>0),FALSE, TRUE)
  if(tst || x < 0){
    print("The value of `x` must be a non-negative whole number")
    stop()
  }
  if(eta<=0){
    print("The value of `eta` must be greater than 0.")
    stop()
  }
  
  p <- 2*(mu*(1/eta+1))^((x+1/eta+2)/2)*besselK(2*sqrt(mu*(1/eta+1)),
                                                x-1/eta-2)/(factorial(x)*gamma(1/eta+2))
  
  if (log) return(log(p))
  else return(p)
})

#' @rdname PoissonInverseGamma
#' @export
ppinvgamma <- Vectorize(function(q, mu=1, eta = 1, lower.tail=TRUE, log.p=FALSE){
  y <- seq(0,q,1)
  probs <- dpinvgamma(y, mu, eta)
  p <- sum(probs)

  if(!lower.tail) p <- 1-p

  if (log.p) return(log(p))
  else return(p)
})

#' @rdname PoissonInverseGamma
#' @export
qpinvgamma <- Vectorize(function(p, mu=1, eta = 1) {
  y <- 0
  p_value <- ppinvgamma(y, mu, eta)
  while(p_value < p){
    y <- y + 1
    p_value <- ppinvgamma(y, mu, eta)
  }
  return(y)
})


#' @rdname PoissonInverseGamma
#' @export
rpinvgamma <- function(n, mu=1, eta = 1) {
  u <- runif(n)
  y <- sapply(u, function(p) qpinvgamma(p, mu, eta))
  return(y)
}
