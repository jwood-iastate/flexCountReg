#' Conway-Maxwell-Poisson (COM) Distribution
#'
#' These functions provide the density function, distribution function, quantile function, and random number generation for the Sichel (SI) Distribution
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mu optional. Numeric value or vector of mean values for the distribution (the values have to be greater than 0).
#' @param lambda optional. Numeric value or vector of values for the rate parameter of the distribution (the values have to be greater than 0). If `mu` is provided, `lambda` is ignored.
#' @param nu optional. Numeric value or vector of values for the decay parameter of the distribution ((the values have to be greater than 0).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
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
#' The compound Probability Mass Function (PMF) for the Conway-Maxwell-Poisson distribution uses:
#' \deqn{f(x|\lambda, \nu)=\frac{\lambda^x}{(x!)^\nu Z(\lambda,\nu)}}
#' #' 
#' Where \eqn{\lambda} and \eqn{\nu} are distribution parameters with \eqn{\lambda>0} and \eqn{\nu>0}, and \eqn{Z(\lambda,\nu)} is the normalizing constant.
#' 
#' The normalizing constant is given by:
#' \deqn{Z(\lambda,\nu)=\sum_{n=0}^{\infty}\frac{\lambda^n}{(n!)^\nu}}
#'
#' The mean and variance of the distribution are given by:
#' \deqn{E[x]=\mu=\lambda \frac{\delta}{\delta \lambda} \log(Z(\lambda,\nu))}
#' \deqn{Var(x)=\lambda \frac{\delta}{\delta \lambda} \mu}
#' 
#' When the mean value is given, the rate parameter (\eqn{\lambda}) is computed using the mean and the decay parameter (\eqn{\nu}). This is useful to allow the calculation of the rate parameter when the mean is known (e.g., in regression))
#' 
#' @examples
#' dcom(1, mu=0.75, nu=3)
#' pcom(c(0,1,2,3,5,7,9,10), lambda=0.75, nu)
#' qcom(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, nu=0.75)
#' rcom(30, mu=0.75, nu=0.5)
#'
#' @importFrom stats factorial uniroot runif
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @export
#' @name COMDistribution

#' @rdname COMDistribution
#' @export
dcom <- Vectorize(function(x, mu=NULL, lambda = 1, nu=1, log=FALSE){
  #test to make sure the value of x is an integer
  tst <- ifelse(is.na(nchar(strsplit(as.character(x), "\\.")[[1]][2])>0),FALSE, TRUE)
  if(tst || x < 0 || lambda<0 || nu<0){
    stop("The value of `x` must be a non-negative whole number, `lambda` must be greater than 0, and `nu` must be greater than 0.")
  }
  
  if(!is.null(mu)){
    lambda <- find_lambda(mean, nu, maxval=1000, tol=1e-8, max_iter=1000)
  }
  
  p <- lambda^x/(factorial(x)^nu * cmp_normalizer_cpp(lambda, nu, maxval=1000))
  
  if (log) return(log(p))
  else return(p)
})

#' @rdname COMDistribution
#' @export
pcom <- Vectorize(function(q, mu=NULL, lambda = 1, nu=1, lower.tail=TRUE, log.p=FALSE){
  
  y <- seq(0,q,1)
  probs <- dcom(y, mu, lambda, nu)
  p <- sum(probs)

  if(!lower.tail) p <- 1-p

  if (log.p) return(log(p))
  else return(p)
})

#' @rdname COMDistribution
#' @export
qcom<- Vectorize(function(p, mu=NULL, lambda = 1, nu=1) {
  y <- 0
  p_value <- pcom(y, mu, lambda, nu)
  while(p_value < p){
    y <- y + 1
    p_value <- pcom(y, mu, lambda, nu)
  }
  return(y)
})

#' @rdname COMDistribution
#' @export
rcom <- function(n, mu=NULL, lambda = 1, nu=1) {
  u <- runif(n)
  y <- sapply(u, qcom, mu=mu, lambda=lambda, nu=nu)
  return(y)
}

# Helper function to get lambda
find_lambda <- function(mu, nu, maxval=100, tol=1e-8, max_iter=100) {
  # Function to compute the difference between computed mean and target mean
  func_to_solve <- function(lambda) {
    expected <- com_expect_cpp(lambda, nu, maxval)
    return(expected - mu)^2
  }
  
  # Start with initial guesses for lambda
  lambda_lower <- 1e-8
  lambda_upper <- mu + 10 * mu  # heuristic upper bound
  
  # Adjust upper bound if necessary
  while (func_to_solve(lambda_upper) < 0 && lambda_upper < 1e8) {
    lambda_upper <- lambda_upper * 2
  }
  
  # Use uniroot to find the root efficiently
  result <- uniroot(
    f = func_to_solve,
    lower = lambda_lower,
    upper = lambda_upper,
    tol = tol,
    maxiter = max_iter
  )
  
  return(result$root)
}
