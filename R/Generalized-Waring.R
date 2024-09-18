#' Generalized Waring Distribution
#'
#' These functions provide density, distribution function, quantile function, and random number generation for the Generalized Waring Distribution.
#'
#' The Generalized Waring distribution is a 3-parameter count distribution that is used to model overdispersed count data.
#'
#' @param y non-negative integer vector of count outcomes.
#' @param q non-negative integer vector of quantiles.
#' @param p numeric vector of probabilities.
#' @param n integer number of random numbers to generate.
#' @param mu numeric vector of means of the distribution.
#' @param k non-negative numeric parameter of the distribution.
#' @param rho non-negative numeric parameter of the distribution.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dgwar} computes the density (PMF) of the Generalized Waring Distribution.
#'
#' \code{pgwar} computes the CDF of the Generalized Waring Distribution.
#'
#' \code{qwaring} computes the quantile function of the Generalized Waring Distribution.
#'
#' \code{rwaring} generates random numbers from the Generalized Waring Distribution.
#'
#' The Probability Mass Function (PMF) for the Generalized Waring (GW) distribution is:
#' \deqn{f(y|a_x,k,\rho)=\frac{\Gamma(a_x+\rho)\Gamma(k+\rho)\left(a_x\right)_y(k)_y}{y!\Gamma(\rho)\Gamma(a_x+k+\rho)(a_x+k+\rho)_y}}
#' Where \eqn{(\alpha)_r=\frac{\Gamma(\alpha+r)}{\Gamma(\alpha)}}, and \eqn{a_x, \ k, \ \rho)>0}.
#'
#' The mean value is:
#' \deqn{E[Y]=\frac{a_x K}{\rho-1}}
#'
#' Thus, we can use:
#' \deqn{a_x=\frac{\mu(\rho-1)}{k}}
#' 
#' This results in a regression model where:
#' \deqn{\mu=e^{X\beta}}
#' \deqn{\sigma^2=\mu\left(1-\frac{1}{\alpha+\rho+1}\right)+\mu^2\frac{(\alpha+\rho)^2}{\alpha\rho(\alpha+\rho+1)}}
#'
#' @examples
#' dgwar(0, mu=1, k=2, rho=3)
#' pgwar(c(0,1,2,3), mu=1, k=2, rho=3)
#' qgwar(0.8, mu=1, k=2, rho=3)
#' rgwar(10, mu=1, k=2, rho=3)
#'
#' @import stats
#' @export
#' @name Generalized-Waring
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg

#' @rdname Generalized-Waring
#' @export
dgwar <- Vectorize(function(y, mu, k, rho, log = FALSE) {

  pmf <- genWaring_cpp(y, mu, k, rho)
  
  if (log) pmf <- log(pmf)
  return(pmf)
})


#' @rdname Generalized-Waring
#' @export
pgwar <- Vectorize(function(q, mu, k, rho, lower.tail = TRUE, log.p = FALSE) {
  if (any(q < 0) || !all(q == floor(q))) {
    stop("All q values must be non-negative integers.")
  }
  y <- seq(0, q, 1)
  probs <- dgwar(y, mu, k = k, rho = rho)
  p <- sum(probs)
  
  if (!lower.tail) p <- 1 - p
  
  if (log.p) return(log(p))
  else return(p)
})


#' @rdname Generalized-Waring
#' @export
qgwar <- Vectorize(function(p, mu, k, rho) {
  if (any(p < 0 | p > 1)) {
    stop("All p values must be in the interval [0, 1].")
  }
  
  y <- 0
  p_value <- pgwar(y, mu, k, rho)
  while (p_value < p) {
    y <- y + 1
    p_value <- pgwar(y, mu, k, rho)
  }
  return(y)
})


#' @rdname Generalized-Waring
#' @export
rgwar <- function(n, mu, k, rho) {
  p <- runif(n)
  random_counts <- qgwar(p, mu, k, rho)
  return(random_counts)
}