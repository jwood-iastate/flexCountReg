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
#' @returns dcom gives the density, pcom gives the distribution function, qcom 
#'  gives the quantile function, and rcom generates random deviates.
#' 
#'  The length of the result is determined by n for rcom, and is the maximum of 
#'  the lengths of the numerical arguments for the other functions.
#' 
#'  The numerical arguments other than n are recycled to the length of the 
#'  result. Only the first elements of the logical arguments are used.
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
  
  pars <- .com_params(N_obs, mu, lambda, nu, "x")
  lambda <- pars$lambda
  nu <- pars$nu
  
  probabilities <- dcom_vec_cpp(x, lambda, nu, log)
  
  return(probabilities)
}

#' @rdname COMDistribution
#' @export
pcom <- function(q, mu = NULL, lambda = 1, nu = 1, lower.tail = TRUE,
                 log.p = FALSE){
  # Same recycling rule as dcom() and qcom(): length 1 or length(q)
  for (nm in c("mu", "lambda", "nu")) {
    v <- get(nm)
    if (!is.null(v) && !length(v) %in% c(1L, length(q))) {
      stop("`", nm, "` must be a single value or a vector with the same ",
           "length as `q`.", call. = FALSE)
    }
  }
  .fcr_p_contract(.pcom_core, q, list(mu = mu, lambda = lambda, nu = nu),
                  lower.tail, log.p)
}

.pcom_core <- function(q, mu = NULL, lambda = 1, nu = 1, lower.tail = TRUE,
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
  pars <- .com_params(N_obs, mu, lambda, nu, "q")
  lambda <- pars$lambda
  nu <- pars$nu
  
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
  pars <- .com_params(N_obs, mu, lambda, nu, "p")
  lambda <- pars$lambda
  nu <- pars$nu
  
  quantiles <- qcom_vec_cpp(p, lambda, nu)
  
  return(quantiles)
  
}

#' @rdname COMDistribution
#' @export
rcom <- function(n, mu = NULL, lambda = 1, nu = 1) {
  n <- .fcr_rn(n)
  nu <- rep_len(nu, n)
  if (!is.null(mu)) {
    mu <- rep_len(mu, n)
    bad <- !(is.finite(mu) & mu > 0 & is.finite(nu) & nu > 0)
    .fcr_warn_na(bad, "`mu` and `nu` must be greater than 0.")
    return(.fcr_r_grouped(n, list(mu = mu, nu = nu), !bad,
                          function(k, mu, nu) {
                            rcom_cpp(k, find_lambda_cpp(mu, nu), nu)
                          }))
  }
  lambda <- rep_len(lambda, n)
  bad <- !(is.finite(lambda) & lambda > 0 & is.finite(nu) & nu > 0)
  .fcr_warn_na(bad, "`lambda` and `nu` must be greater than 0.")
  .fcr_r_grouped(n, list(lambda = lambda, nu = nu), !bad,
                 function(k, lambda, nu) rcom_cpp(k, lambda, nu))
}

# Recycle the COM parameters to length n (base-R style: length 1 or n), and
# find the rate lambda from the mean mu when mu is given. lambda is solved
# once per distinct (mu, nu) pair.
.com_params <- function(n, mu, lambda, nu, arg) {
  check <- function(v, nm) {
    if (!length(v) %in% c(1L, n)) {
      stop("`", nm, "` must be a single value or a vector with the same ",
           "length as `", arg, "`.", call. = FALSE)
    }
  }
  check(nu, "nu")
  nu <- rep_len(nu, n)
  if (!is.null(mu)) {
    check(mu, "mu")
    mu <- rep_len(mu, n)
    key <- paste(mu, nu)
    first <- !duplicated(key)
    lambda <- find_lambda_vec_cpp(mu[first], nu[first])[match(key, key[first])]
  } else {
    check(lambda, "lambda")
    lambda <- rep_len(lambda, n)
  }
  list(lambda = lambda, nu = nu)
}
