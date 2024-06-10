#' Triangle Distribution
#'
#' These functions provide density, distribution function, quantile function,
#' and random number generation for the Triangle Distribution, specified by its
#' mean, standard deviation, and optional lower and upper bounds.
#'
#' The Triangle Distribution is defined by three points: a (minimum), b
#' (maximum), and c (mode), where the density is zero outside the interval [a,
#' b], increases linearly from a to c, and decreases linearly from c to b.
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mode numeric value or vector of mode values for the distribution.
#' @param sigma single value or vector indicating both the positive and negative
#'   max differences from the mean (if the difference is the same).
#' @param upper single value or vector for the upper limit of the distribution
#'   (must be used with `lower`).
#' @param lower single value or vector for the lower limit of the distribution
#'   (must be used with `upper`).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dtri} computes the density (PDF) of the Triangle Distribution.
#'
#' \code{ptri} computes the CDF of the Triangle Distribution.
#'
#' \code{qtri} computes the quantile function of the Triangle Distribution.
#'
#' \code{rtri} generates random numbers from the Triangle Distribution.
#'
#' The mode and standard deviation parameters define the distribution's location
#' and scale, respectively, while the lower and upper bounds explicitly set the
#' minimum and maximum values of the distribution.
#'
#' @examples
#' dtri(4, mode=8, upper=13, lower=1)
#' ptri(c(0, 1, 2, 3, 5, 7, 9, 10), mode = 3, upper=9, lower = 1)
#' qtri(c(0.1, 0.3, 0.5, 0.9, 0.95), mode = 3, upper = 9, lower = 1)
#' rtri(30, mode = 5, sigma = 3)
#'
#' @import stats
#' @export
#' @name Triangular

#' @rdname Triangular
#' @export
dtri <- Vectorize(function(
    x, mode = 0, sigma = 1, upper = NULL, lower = NULL, log = FALSE){

  c <- mode

  if (!is.null(upper) && !is.null(lower)) {
    if (lower >= upper){
      msg <-'The value of `lower` must be smaller than the value of `upper`' 
      stop(msg)
    }
    else if (lower > c){
      msg <- paste0('The value of `mode` must be greater than or ', 
                    'equal to the value of `lower`')
        
      stop(msg)
    }
    else if (c > upper){
      msg <- paste0('The value of `mode` must be smaller than or ',
                    'equal to the value of `upper`') 
      stop(msg)
    }
    a <- lower
    b <- upper
  }
  else{
    a <- c - sigma
    b <- c + sigma
  }
  
  # Compute probabilities
  if (x <= a | x >= b) {
    p <- 0
  } else if (x < c) {
    p <- 2 * (x - a) / ((b - a) * (c - a))
  }
  else if (x == c){
    p <- 2 / (b - a)
  }
  else {
    p <- 2 * (b - x) / ((b - a) * (b - c))
  } 
  
  if (log) return(log(p))
  else return(p)
})


#' @rdname Triangular
#' @export
ptri <- Vectorize(function(
    q, mode=0, sigma = 1, upper = NULL, lower=NULL, 
    lower.tail=TRUE, log.p=FALSE){
  c <- mode
  if (!is.null(upper) && !is.null(lower)){
    if (lower>=upper){
      print('The value of `lower` must be smaller than the value of `upper`')
      stop()
    }
    else if (lower>c){
      print('The value of `mode` must be greater than or equal to the value of `lower`')
      stop()
    }
    else if (c>upper){
      print('The value of `mode` must be smaller than or equal to the value of `upper`')
      stop()
    }
    a <- lower
    b <- upper
  }
  else{
    a <- c-sigma
    b <- c+sigma
  }
  # compute cumulative probabilities
  if(q<a){
    p <- 0
  }
  else if(q>b){
    p <- 1
  }
  else if (q<=c){
    p <- (q-a)^2/((b-a)*(c-a))
  }
  else{
    p <- 1-(b-q)^2/((b-a)*(b-c))
  }
  if(!lower.tail) p <- 1-p
  
  if (log.p) return(log(p))
  else return(p)
})

#' @rdname Triangular
#' @export
qtri <- Vectorize(function(p, mode = 0, sigma = 1, upper = NULL, lower = NULL) {
  c <- mode
  if (!is.null(upper) && !is.null(lower)){
    if (lower>=upper){
      print('The value of `lower` must be smaller than the value of `upper`')
      stop()
    }
    else if (lower>c){
      print('The value of `mode` must be greater than or equal to the value of `lower`')
      stop()
    }
    else if (c>upper){
      print('The value of `mode` must be smaller than or equal to the value of `upper`')
      stop()
    }
    a <- lower
    b <- upper
    p_mode <- ptri(q=mode, mode=mode, upper=upper, lower=lower)
  } else {
    a <- c - sigma
    b <- c + sigma
    p_mode <- ptri(mode, mode=mode, sigma=sigma)
  }
  if (p<=p_mode) {
    q <- a + sqrt(p * (b - a) * (c - a))
  } else {
    q <- b - sqrt((1 - p) * (b - a) * (b - c))
  }
  return(q)
})

#' @rdname Triangular
#' @export
rtri <- function(n, mode = 0, sigma = 1, upper = NULL, lower = NULL) {
  if(is.null(upper) || is.null(lower)){
    return(qtri(stats::runif(n), mode, sigma))
  }
  else(
    return(qtri(stats::runif(n), mode, upper, lower))
  )
}
