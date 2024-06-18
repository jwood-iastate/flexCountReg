#' Poisson-Inverse-Gaussian Distribution
#'
#' These functions provide the density function, distribution function, quantile function, and random number generation for the Poisson-Inverse-Gaussian (PInvGaus) Distribution
#'
#' The Poisson-Inverse-Gaussian distribution is a special case of the Sichel distribution, as noted by Cameron & Trivedi (2013). It is also known as a uivariate Sichel distribution (Hilbe, 2011).
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mu numeric value or vector of mean values for the distribution (the values have to be greater than 0).
#' @param eta single value or vector of values for the scale parameter of the distribution (the values have to be greater than 0).
#' @param form optional parameter indicating which formulation to use. Options include "Type 1" which is the standard form or "Type 2" which follows the formulation by Dean et. al. (1987).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dpinvgaus} computes the density (PDF) of the Poisson-Inverse-Gaussian Distribution.
#'
#' \code{ppinvgaus} computes the CDF of the Poisson-Inverse-Gaussian Distribution.
#'
#' \code{qpinvgaus} computes the quantile function of the Poisson-Inverse-Gaussian Distribution.
#'
#' \code{rpinvgaus} generates random numbers from the Poisson-Inverse-Gamma Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Poisson-Inverse-Gaussian distribution (Type 1) is (Cameron & Trivedi, 2013):
#' \deqn{f(y|\eta,\mu)=\begin{cases}
#'                              f(y=0)=exp\left(\frac{\mu}{\eta}\left(1-\sqrt{1+2\eta}\right)\right) \\
#'                              f(y|y>0)=f(y=0)\frac{\mu^y}{y!}(1+2\eta)^{-y/2}\cdot\sum_{j=0}^{y-1}\frac{\Gamma(y+j)}{\Gamma(y-j)\Gamma(j+1)}\left(\frac{\eta}{2\mu}\right)^2(1+2\eta)^{-j/2}
#'                              \end{cases}}
#' 
#' Where \eqn{\eta} is a scale parameter with the restriction that \eqn{eta>0}, \eqn{\mu} is the mean value, and \eqn{y} is a non-negative integer.
#'
#' The variance of the distribution is:
#' \deqn{\sigma^2=\mu+\eta\mu}
#' 
#' The alternative parameterization by Dean et. al. (1987) replaces \eqn{\eta} with \eqn{\eta\mu}. This version(Type 2) has the PMF:
#' \deqn{f(y|\eta,\mu)=\begin{cases}
#'                              f(y=0)=exp\left(\frac{1}{\eta}\left(1-\sqrt{1+2\eta\mu}\right)\right) \\
#'                              f(y|y>0)=f(y=0)\frac{\mu^y}{y!}(1+2\eta\mu)^{-y/2}\cdot\sum_{j=0}^{y-1}\frac{\Gamma(y+j)}{\Gamma(y-j)\Gamma(j+1)}\left(\frac{\eta}{2}\right)^2(1+2\eta\mu)^{-j/2}
#'                              \end{cases}}
#' 
#'  This results in the variance of:
#' \deqn{\sigma^2=\mu+\eta\mu^2}
#'
#' @references 
#' Cameron, A. C., & Trivedi, P. K. (2013). Regression analysis of count data, 2nd Edition. Cambridge university press.
#' 
#' Dean, C., Lawless, J. F., & Willmot, G. E. (1989). A mixed Poisson–Inverse‐Gaussian regression model. Canadian Journal of Statistics, 17(2), 171-181.
#' 
#' Hilbe, J. M. (2011). Negative binomial regression. Cambridge University Press.
#' 
#' @examples
#' dpinvgaus(1, mu=0.75, eta=1)
#' ppinvgaus(c(0,1,2,3,5,7,9,10), mu=0.75, eta=3, form="Type 2")
#' qpinvgaus(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, eta=0.5, form="Type 2")
#' rpinvgaus(30, mu=0.75, eta=1.5)
#'
#' @import stats
#' @export
#' @name Poisson-Inverse-Gaussian

#' @rdname Poisson-Inverse-Gaussian
#' @export
dpinvgaus <- Vectorize(function(x, mu=1, eta = 1, form="Type 1", log=FALSE){
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
  if(form=='Type 2'){
    eta <- eta*mu # make adjusted value for Type 2 formulation
  }
  p0 <- exp(mu/eta*(1-sqrt(1+2*eta)))
  
  if(x>0){
    if(x>1){
      j <- seq(0,(x-1),1)
    }
    else{
      j <- 0
    }
    e2 <- sum(gamma(x+j)/(gamma(x-j)*gamma(j+1))*(eta/(2*mu))^j *(1+2*eta)^(-j/2))
    p <- p0 * (mu^x)/gamma(x+1)*(1+2*eta)^(-x/2) * e2
  }
  else{
    p <- p0
  }
  if (log) return(log(p))
  else return(p)
})

#' @rdname Poisson-Inverse-Gaussian
#' @export
ppinvgaus <- Vectorize(function(q, mu=1, eta = 1, form="Type 1", lower.tail=TRUE, log.p=FALSE){
  y <- seq(0,q,1)
  probs <- dpinvgaus(y, mu, eta, form)
  p <- sum(probs)

  if(!lower.tail) p <- 1-p

  if (log.p) return(log(p))
  else return(p)
})

#' @rdname Poisson-Inverse-Gaussian
#' @export
qpinvgaus <- Vectorize(function(p, mu=1, eta = 1, form="Type 1") {
  y <- 0
  p_value <- ppinvgaus(y, mu, eta, form)
  while(p_value < p){
    y <- y + 1
    p_value <- ppinvgaus(y, mu, eta, form)
  }
  return(y)
})


#' @rdname Poisson-Inverse-Gaussian
#' @export
rpinvgaus <- function(n, mu=1, eta = 1, form="Type 1") {
  u <- runif(n)
  y <- sapply(u, function(p) qpinvgaus(p, mu, eta, form))
  return(y)
}
