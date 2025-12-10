#' Conway-Maxwell-Poisson (COM) Distribution
#'
#' These functions provide the density function, distribution function,
#' quantile function, and random number generation for the
#' Conway-Maxwell-Poisson (COM) Distribution
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mu optional. Numeric value or vector of mean values for the
#'   distribution (the values have to be greater than 0).
#' @param lambda optional. Numeric value or vector of values for the rate
#'   parameter of the distribution (the values have to be greater than 0).
#'   If `mu` is provided, `lambda` is ignored.
#' @param nu optional. Numeric value or vector of values for the decay
#'   parameter of the distribution ((the values have to be greater than 0).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dcom} computes the density (PDF) of the COM Distribution.
#'
#' \code{pcom} computes the CDF of the COM Distribution.
#'
#' \code{qcom} computes the quantile function of the COM Distribution.
#'
#' \code{rcom} generates random numbers from the COM Distribution.
#'
#' The Probability Mass Function (PMF) for the Conway-Maxwell-Poisson
#' distribution is:
#' \deqn{f(x|\lambda, \nu) = \frac{\lambda^x}{(x!)^\nu Z(\lambda,\nu)}}
#'
#' Where \eqn{\lambda} and \eqn{\nu} are distribution parameters with
#' \eqn{\lambda>0} and \eqn{\nu>0}, and \eqn{Z(\lambda,\nu)} is the
#' normalizing constant.
#'
#' The normalizing constant is given by:
#' \deqn{Z(\lambda,\nu)=\sum_{n=0}^{\infty}\frac{\lambda^n}{(n!)^\nu}}
#'
#' The mean and variance of the distribution are given by:
#' \deqn{E[x]=\mu=\lambda \frac{\delta}{\delta \lambda} \log(Z(\lambda,\nu))}
#' \deqn{Var(x)=\lambda \frac{\delta}{\delta \lambda} \mu}
#'
#' When the mean value is given, the rate parameter (\eqn{\lambda}) is
#' computed using the mean and the decay parameter (\eqn{\nu}). This is
#' useful to allow the calculation of the rate parameter when the mean is
#' known (e.g., in regression))
#'
#' @examples
#' dcom(1, mu=0.75, nu=3)
#' pcom(c(0,1,2,3,5,7,9,10), lambda=0.75, nu=0.75)
#' qcom(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, nu=0.75)
#' rcom(30, mu=0.75, nu=0.5)
#'
#' @importFrom stats uniroot runif
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @export
#' @name COMDistribution

#' @rdname COMDistribution
#' @export
dcom <- function(x, mu = NULL, lambda = 1, nu = 1, log = FALSE){
  # test to make sure the value of x is an integer
  if (any(!is.numeric(x) | x < 0 | floor(x) != x)) {
    stop("The value of `x` must be a non-negative whole number.")
  }
  x <- as.integer(x)
  N_obs <- length(x)
  
  # ensure that the vectors of parameter values is the same length as
  # the length of x
  if (any(lambda <= 0) | any(nu <= 0) | any(mu <=0))
    stop("The values of `mu`, lambda`, and `nu` must all be greater than 0.")
  
  if (N_obs > 1){
    if(!is.null(mu)){
      if (length(mu)==1){
        if (length(nu)==1){
          lambda <- find_lambda_cpp(mu, nu)
          lambda <- rep(lambda, N_obs)
          nu <- rep(nu, N_obs)
        } else if(length(nu)==N_obs){
          mu <- rep(mu, N_obs)
          lambda <- find_lambda_vec_cpp(mu, nu)
        } else {
          msg <- paste("`nu` must be a single value or a vector",
                       "with the same length as `x`")
          warning(msg)
        }
      } else if (length(mu) == N_obs){
        if (length(nu)==1){
          nu <- rep(nu, N_obs)
        } else if (length(nu)!=N_obs){
          msg <- paste("`nu` must be a single value or",
                       "a vector with the same length as `x`")
          warning(msg)
        }
        lambda <- find_lambda_vec_cpp(mu, nu)
      } else {
        msg <- 
          "`mu` must be a single value or a vector with the same length as `x`."
        warning(msg)
      }
    } else {
      if (length(lambda)==1){
        lambda <- rep(lambda, N_obs)
      } else {
        if (length(lambda)!=N_obs) {
          msg <- paste("`lambda` must be a single value or a vector",
                       "with the same length as `x`")
          warning(msg)
        }
      }
    }
    if(length(nu)==1 & N_obs>1){
      nu <- rep(nu, N_obs)
    } else {
      if (length(nu) != N_obs) {
        msg <- paste("`nu` must be a single value or a",
                     "vector with the same length as `x`")
        warning(msg)
      }
    }
  } else{ # if x is a single value and mu is provided
    if (!is.null(mu) & length(mu) == 1 & length(nu) == 1)
      lambda <- find_lambda_cpp(mu, nu)
  }
  
  probabilities <- dcom_vec_cpp(x, lambda, nu, log)
  
  return(probabilities)
}

#' @rdname COMDistribution
#' @export
pcom <- function(q, mu = NULL, lambda = 1, nu = 1, lower.tail = TRUE,
                 log.p = FALSE){
  
  if (any(!is.numeric(q)) | any(q < 0) | any(floor(q) != q))
    warning("The value of `q` must be a non-negative whole number.")
  
  if (any(lambda <= 0) | any(nu <= 0) | any(mu <=0)) {
    msg <- "The values of `mu`, lambda`, and `nu` must all be greater than 0."
    warning(msg)
  }
  
  q <- as.integer(q)
  N_obs <- length(q)
  
  # ensure that the vectors of parameter values is the same length as
  # the length of q
  if (N_obs > 1){
    if(!is.null(mu)){
      if (length(mu)==1){
        if (length(nu)==1){
          lambda <- find_lambda_cpp(mu, nu)
          lambda <- rep(lambda, N_obs)
          nu <- rep(nu, N_obs)
        } else if(length(nu)==N_obs){
          mu <- rep(mu, N_obs)
          lambda <- find_lambda_vec_cpp(mu, nu)
        } else{
          msg <- paste("`nu` must be a single value or a vector",
                       "with the same length as `q`")
          warning(msg)
        }
      } else if (length(mu)==N_obs){
        if (length(nu)==1){
          nu <- rep(nu, N_obs)
        } else if (length(nu)!=N_obs){
          msg <- paste("`nu` must be a single value or a",
                       "vector with the same length as `q`")
          warning(msg)
        }
        lambda <- find_lambda_vec_cpp(mu, nu)
      } else {
        msg <- paste("`mu` must be a single value or a vector",
                     "with the same length as `q`.")
        warning(msg)
      }
    } else {
      if (length(lambda) == 1){
        lambda <- rep(lambda, N_obs)
      } else {
        if (length(lambda) != N_obs) {
          msg <- paste("`lambda` must be a single value or a",
                       "vector with the same length as `q`")
          warning(msg)
        }
      }
    }
    if(length(nu)==1 & N_obs>1){
      nu <- rep(nu, N_obs)
    } else {
      if (length(nu) != N_obs) {
        msg <- paste("`nu` must be a single value or a",
                     "vector with the same length as `q`")
        warning(msg)
      }
    }
  } else{ # if q is a single value and mu is provided
    if (!is.null(mu) & length(mu) == 1 & length(nu) == 1)
      lambda <- find_lambda_cpp(mu, nu)
  }
  
  p <- pcom_vec_cpp(q, lambda, nu, lower.tail, log.p)
  
  return(p)
}

#' @rdname COMDistribution
#' @export
qcom <- function(p, mu = NULL, lambda = 1, nu = 1) {
  if (any(p <= 0) | any(p >= 1)) {
    msg <- "The values for `p` must be greater than 0 and less than 1."
    warning(msg)
  }
  
  if (any(lambda <= 0) | any(nu <= 0) | any(mu <=0)) {
    warning("The values of `mu`, lambda`, and `nu` must all be greater than 0.")
  }
  
  N_obs <- length(p)
  
  # ensure that the vectors of parameter values is the same length as
  # the length of p
  if (N_obs > 1){
    if(!is.null(mu)){
      if (length(mu)==1){
        if (length(nu)==1){
          lambda <- find_lambda_cpp(mu, nu)
          lambda <- rep(lambda, N_obs)
          nu <- rep(nu, N_obs)
        } else if(length(nu)==N_obs){
          mu <- rep(mu, N_obs)
          lambda <- find_lambda_vec_cpp(mu, nu)
        } else{
          msg <- paste("`nu` must be a single value or a vector",
                       "with the same length as `p`")
          warning(msg)
        }
      } else if (length(mu)==N_obs){
        if (length(nu)==1){
          nu <- rep(nu, N_obs)
        } else if (length(nu)!=N_obs){
          msg <- paste("`nu` must be a single value or a",
                       "vector with the same length as `p`")
          warning(msg)
        }
        lambda <- find_lambda_vec_cpp(mu, nu)
      } else {
        msg <- paste("`mu` must be a single value or a vector",
                     "with the same length as `p`.")
        warning(msg)
      }
    } else {
      if (length(lambda) == 1) {
        lambda <- rep(lambda, N_obs)
      } else {
        if (length(lambda) != N_obs) {
          msg <- paste("`lambda` must be a single value or a vector",
                       "with the same length as `p`")
          warning(msg)
        }
      }
    }
    if(length(nu) == 1 & N_obs > 1){
      nu <- rep(nu, N_obs)
    } else {
      if (length(nu)!=N_obs) {
        msg <- "`nu` must be a single value or a vector with the same length"
        warning(msg)
      }
    }
  } else{
    # if p is a single value and mu is provided
    if (!is.null(mu) & length(mu) == 1 & length(nu) == 1)
      lambda <- find_lambda_cpp(mu, nu)
  }
  
  quantiles <- qcom_vec_cpp(p, lambda, nu)
  
  return(quantiles)
  
}

#' @rdname COMDistribution
#' @export
rcom <- function(n, mu = NULL, lambda = 1, nu = 1) {
  if (length(n) != 1 || !is.numeric(n) || n <= 0 || floor(n) != n)
    warning("`n` must be a single positive integer.")
  
  if ((lambda <= 0) | (nu <= 0) | (mu <=0)) {
    msg <- "The values of `mu`, lambda`, and `nu` must all be greater than 0."
    warning(msg)
  }
  
  if (!is.null(mu)) {
    lambda <- find_lambda_cpp(mu, nu)
  }
  
  random_numbers <- rcom_cpp(n, lambda, nu)
  return(random_numbers)
}
