#' @title One-Parameter Lindley Distribution
#' 
#' @description Distribution function for the one-parameter Lindley distribution with parameter theta.
#'
#' @param x a single value or vector of positive values.
#' @param p a single value or vector of probabilities.
#' @param q a single value or vector of quantiles.
#' @param n number of random values to generate.
#' @param theta distribution parameter value. Default is 1.
#' @param log,log.p logical; If TRUE, probabilities p are given as log(p). If FALSE, probabilities p are given directly. Default is FALSE.
#' @param lower.tail logical; If TRUE, (default), \eqn{P(X \leq x)} are returned, otherwise \eqn{P(X > x)} is returned. Default is TRUE.
#'
#' @details
#' Probability density function (PDF)
#' \deqn{f(x\mid \theta )=\frac{\theta ^{2}}{(1+\theta )}(1+x)e^{-\theta x}}
#'
#' Cumulative distribution function (CDF)
#' \deqn{F(x\mid \theta ) =1 - \left(1+ \frac{\theta x}{1+\theta }\right)e^{-\theta x}}
#'
#' Quantile function (Inverse CDF)
#' \deqn{Q(p\mid \theta )=-1-\frac{1}{\theta }-\frac{1}{\theta }W_{-1}\left((1+\theta)( p-1)e^{-(1+\theta) }\right)}
#'
#' where \eqn{W_{-1}()} is the negative branch of the Lambert W function.
#' 
#' The moment generating function (MGF) is:
#' \deqn{M_X(t)=\frac{\theta^2(\theta-t+1)}{(\theta+1)(\theta-t)^2}}
#' 
#' The distribution mean and variance are:
#' \deqn{\mu=\frac{\theta+2}{\theta(1+\theta)}}
#' \deqn{\sigma^2=\frac{\mu}{\theta+2}\left(\frac{6}{\theta}-4\right)-\mu^2}
#' 
#' @importFrom lamW lambertWm1
#'
#' @examples
#' x <- seq(0, 5, by = 0.1)
#' p <- seq(0.1, 0.9, by = 0.1)
#' q <- c(0.2, 3, 0.2)
#' dlindley(x, theta = 1.5)
#' dlindley(x, theta=0.5, log=TRUE)
#' plindley(q, theta = 1.5)
#' plindley(q, theta = 0.5, lower.tail = FALSE)
#' qlindley(p, theta = 1.5)
#' qlindley(p, theta = 0.5)
#' 
#' set.seed(123154)
#' rlindley(5, theta = 1.5)
#' rlindley(5, theta = 0.5)
#'
#' @rdname Lindley
#' @export
dlindley <- Vectorize(function(x, theta=1, log = FALSE){
  if (theta<=0) stop("Requirements: theta > 0 & x > 0")
  p <- theta^2 / (1 + theta) * (1+x) * exp(-theta * x)
  if(log) return(log(p)) else return(p)
})

#' @rdname Lindley
#' @export
plindley <- Vectorize(function(q, theta=1, lower.tail = TRUE, log.p = FALSE){
  if (theta<=0 || q<=0) stop("Requirements: theta > 0 & q > 0")
  CDF <- 1-(1+theta+theta*q)/(1+theta)*exp(-theta*q)
  UpperTail <- 1-CDF
  return_val <- ifelse(lower.tail, CDF, UpperTail)
  if(log.p) return(log(return_val)) else return(return_val)
})

#' @rdname Lindley
#' @export
qlindley <- Vectorize(function(p, theta=1, log.p = FALSE){
  if (log.p) p <- exp(p) else p <- p
  if (theta<=0 || p<=0 || p>=1) stop("Requirements: theta > 0 & p > 0 & p < 1")
  
  qval <- -1-1/theta-lambertWm1((1+theta)*(p-1)*exp(-(1+theta)))/theta
  return(qval)
})

#' @rdname Lindley
#' @export
rlindley <- function(n, theta=1){
  if (theta<=0 || n<=0) stop("Requirements: theta > 0 & n > 0")
  qlindley(runif(n), theta)
}
