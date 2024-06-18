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
#' @param alpha non-negative numeric parameter of the distribution.
#' @param rho non-negative numeric parameter of the distribution.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dwaring} computes the density (PMF) of the Generalized Waring Distribution.
#'
#' \code{pwaring} computes the CDF of the Generalized Waring Distribution.
#'
#' \code{qwaring} computes the quantile function of the Generalized Waring Distribution.
#'
#' \code{rwaring} generates random numbers from the Generalized Waring Distribution.
#'
#' The Probability Mass Function (PMF) for the Generalized Waring (GW) distribution is:
#' \deqn{PMF=\frac{\Gamma(\alpha+\rho)}{\Gamma(\alpha)\Gamma(\rho)}
#'        \left(\frac{\alpha}{\mu(\rho-1)}\right)^{\alpha-1}
#'        \left(1+\frac{\alpha}{\mu(\rho-1)}\right)^{-(\alpha+\rho)}}
#'
#' Where \eqn{\alpha} and \eqn{\rho} are distribution parameters with the constraints that \eqn{\alpha \geq 0} and \eqn{\rho > 1}, \eqn{\mu} is the mean of the distribution.
#'
#'
#' @examples
#' dgwar(0, mu=1, alpha=2, rho=3)
#' pgwar(c(0,1,2,3), mu=1, alpha=2, rho=3)
#' qgwar(c(0.1, 0.5, 0.9), mu=1, alpha=2, rho=3)
#' rgwar(10, mu=1, alpha=2, rho=3)
#'
#' @import stats
#' @export
#' @name Generalized-Waring

#' @rdname Generalized-Waring
#' @export
dgwar <- function(y, mu, alpha, rho, log = FALSE) {
  if (any(y < 0) || !all(y == floor(y))) {
    stop("All y values must be non-negative integers.")
  }
  if (alpha < 0 || rho <= 1) {
    stop("Alpha must be non-negative and rho must be greater than 1.")
  }
  
  # Compute PMF for each y value
  v <- alpha / (mu * (rho - 1))
  p <- gamma(alpha + rho) / (gamma(alpha) * gamma(rho)) * v^(alpha - 1) * (1 + v)^(-1 * (alpha + rho))
  pmf <- p * sapply(y, function(yi) v^yi / factorial(yi))
  
  if (log) pmf <- log(pmf)
  return(pmf)
}


#' @rdname Generalized-Waring
#' @export
pgwar <- function(q, mu, alpha, rho, lower.tail = TRUE, log.p = FALSE) {
  if (any(q < 0) || !all(q == floor(q))) {
    stop("All q values must be non-negative integers.")
  }
  # Compute the cumulative probability for each quantile
  p <- sapply(q, function(qi) sum(dgwar(0:qi, mu, alpha, rho)))
  
  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  return(p)
}


#' @rdname Generalized-Waring
#' @export
qgwar <- function(p, mu, alpha, rho) {
  if (any(p < 0 | p > 1)) {
    stop("All p values must be in the interval [0, 1].")
  }
  
  quantiles <- sapply(p, function(pi) {
    low <- 0
    high <- 1000  # Set a reasonable upper bound; adjust based on your data's range
    while (low < high) {
      mid <- (low + high) %/% 2
      P <- pgwar(mid, mu, alpha, rho)
      if (P < pi) {
        low <- mid + 1
      } else {
        high <- mid
      }
    }
    return(low)
  })
  return(quantiles)
}


#' @rdname Generalized-Waring
#' @export
rgwar <- function(n, mu, alpha, rho) {
  p <- runif(n)
  random_counts <- qgwar(p, mu, alpha, rho)
  return(random_counts)
}
