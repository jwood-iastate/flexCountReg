#' Poisson-Lindley Distribution
#'
#' These functions provide density, distribution function, quantile function,
#' and random number generation for the Poisson-Lindley (PL) Distribution
#'
#' The Poisson-Lindley is a 2-parameter count distribution that captures high
#' densities for small integer values. This makes it ideal for data that are
#' zero-inflated.
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mean numeric value or vector of mean values for the distribution (the
#'   values have to be greater than 0).
#' @param theta single value or vector of values for the theta parameter of the
#'   distribution (the values have to be greater than 0).
#' @param lambda alternative parameterization (use instead of the mean); numeric
#'   value or vector of values for lambda parameter of the distribution (the
#'   values have to be greater than 0).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dplind} computes the density (PDF) of the Poisson-Lindley Distribution.
#'
#' \code{pplind} computes the CDF of the Poisson-Lindley Distribution.
#'
#' \code{qplind} computes the quantile function of the Poisson-Lindley
#' Distribution.
#'
#' \code{rplind} generates random numbers from the Poisson-Lindley Distribution.
#'
#' The compound Probability Mass Function(PMF) for the Poisson-Lindley (PL)
#' distribution is:
#' \deqn{f(y| \theta, \lambda) = 
#'   \frac{\theta^2 \lambda^y (\theta + \lambda + y + 1)}
#'     {(\theta + 1)(\theta + \lambda)^{y + 2}}}
#'
#' Where \eqn{\theta} and \eqn{\lambda} are distribution parameters with the
#' restrictions that \eqn{\theta > 0} and \eqn{\lambda > 0}, and \eqn{y} is a
#' non-negative integer.
#'
#' The expected value of the distribution is:
#' \deqn{\mu=\frac{\lambda(\theta + 2)}{\theta(\theta + 1)}}
#'
#' The default is to use the input mean value for the distribution. However, the
#' lambda parameter can be used as an alternative to the mean value.
#'
#' @examples
#' dplind(0, mean = 0.75, theta = 7)
#' pplind(c(0, 1, 2, 3, 5, 7, 9, 10), mean = 0.75, theta = 7)
#' qplind(c(0.1, 0.3, 0.5, 0.9, 0.95), lambda = 4.67, theta = 7)
#' rplind(30, mean = 0.75, theta = 7)
#'
#' @import stats
#' @export
#' @name Poisson-Lindley

#' @rdname Poisson-Lindley
#' @export
dplind <- Vectorize(function(
    x, mean = 1, theta = 1, lambda = NULL, log = FALSE){
  
  # Test to make sure the value of x is an integer
  tst <- !(is.na(nchar(strsplit(as.character(x), "\\.")[[1]][2]) > 0))
  if(tst | x < 0){
    msg <- ("The value of `x` must be a non-negative whole number")
    stop(msg)
  }
  if(is.null(lambda)){
    if(mean <= 0 | theta <= 0){
      msg <- paste0('The values of `mean` and `theta` both ',
                    'have to have values greater than 0.')
      stop(msg)
    }
    else{
      lambda <- mean * theta * (theta + 1) / (theta + 2)
      p <- (theta^2 * lambda^x * (theta + lambda + x + 1))/
        ((theta + 1) * (theta + lambda)^(x + 2))
    }
  }
  else{
    if(lambda <= 0 | theta <= 0){
      msg <- paste0('The values of `lambda` and `theta` both',
                    ' have to have values greater than 0.')
      stop(msg)
    }
    else{
      p <- (theta^2 * lambda^x * (theta + lambda + x + 1)) /
        ((theta + 1) * (theta + lambda)^(x + 2))
    }
  }
  if (log) return(log(p))
  else return(p)
})

#' @rdname Poisson-Lindley
#' @export
pplind <- Vectorize(function(
    q, mean = 1, theta = 1, lambda = NULL, lower.tail = TRUE, log.p = FALSE){
  
  if (is.null(lambda)){
    if(mean <= 0 | theta <= 0){
      msg <- paste0('The values of `mean` and `theta` both ',
                    'have to have values greater than 0.')
      stop(msg)
    }
  }
  else{
    if(lambda <= 0 | theta <= 0){
      msg <- paste0('The values of `lambda` and `theta` both have ',
                    'to have values greater than 0.')
      stop(msg)
    }
    else{
      mean <- lambda * (theta + 2) / (theta * (theta + 1))
    }
  }

  y <- seq(0, q, 1)
  probs <- dplind(y, mean, theta)
  p <- sum(probs)

  if(!lower.tail) p <- 1 - p

  if (log.p) return(log(p))
  else return(p)
})

#' @rdname Poisson-Lindley
#' @export
qplind <- Vectorize(function(p, mean=1, theta=1, lambda=NULL) {
  if (is.null(lambda)){
    if(mean<=0 || theta<=0){
      msg <- paste0('The values of `mean` and `theta` both have to ',
                    'have values greater than 0.')
      stop(msg)
    }
  }
  else{
    if (lambda<=0 | theta <= 0) {
      msg <- paste0('The values of `lambda` and `theta` both have to ',
                    'have values greater than 0.')
      stop(msg)
    }
    else{
      mean <- lambda * (theta + 2) / (theta * (theta + 1))
    }
  }

  y <- 0
  p_value <- pplind(y, mean, theta)
  while(p_value < p){
    y <- y + 1
    p_value <- pplind(y, mean, theta)
  }
  return(y)
})


#' @rdname Poisson-Lindley
#' @export
rplind <- function(n, mean=1, theta=1, lambda=NULL) {
  
  if(is.null(lambda)){
    if(mean <= 0 | theta <= 0){
      msg <- paste0('The values of `mean` and `theta` both have to ',
                    'have values greater than 0.')
      stop(msg)
    }
  }
  else{
    if(lambda <= 0 | theta <= 0){
      msg <- paste0('The values of `lambda` and `theta` both have to ',
                    'have values greater than 0.')
      stop(msg)
    }
    else{
      mean <- lambda * (theta + 2) / (theta * (theta + 1))
    }
  }

  u <- runif(n)
  y <- sapply(u, function(p) qplind(p, mean, theta))
  return(y)
}
