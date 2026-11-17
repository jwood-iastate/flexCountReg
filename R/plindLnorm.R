#' Poisson-Lindley-Lognormal Distribution
#'
#' These functions provide density, distribution, quantile, and random
#' generation for the Poisson-Lindley-Lognormal (PLL) Distribution.
#'
#' The PLL is a 3-parameter count distribution that captures high mass at
#' small y and allows flexible heavy tails.
#'
#' @param x numeric value or vector of values.
#' @param q quantile or vector of quantiles.
#' @param p probability or vector of probabilities.
#' @param n number of random draws.
#' @param mean the Poisson-Lindley mean before the lognormal factor (>0); the
#'   mean of the distribution is \code{mean * exp(sigma^2 / 2)}.
#' @param theta Poisson-Lindley theta parameter (>0).
#' @param sigma lognormal sigma parameter (>= 0; 0 gives the
#'   Poisson-Lindley distribution).
#' @param lambda optional Poisson scale parameter (>0), an alternative to
#'   \code{mean}: \code{mean = lambda * (theta + 2) / (theta * (theta + 1))}.
#' @param ndraws the number of Halton draws (only with \code{engine = "halton"}).
#' @param log return log-density.
#' @param log.p return log-CDF.
#' @param lower.tail TRUE returns \eqn{P[X \leq x]}.
#' @param hdraws optional Halton draws (only with \code{engine = "halton"}).
#' @param engine the numerical method: \code{"quadrature"} (default) computes
#'   the log-probability by adaptive quadrature over the log of the mixing
#'   variable, accurate for every count including far tails;
#'   \code{"halton"} averages the Poisson probabilities over Halton draws
#'   (the method of flexCountReg 0.1.2 and earlier, kept for
#'   reproducibility; badly biased far in the tail; the default when
#'   \code{hdraws} are given).
#'
#' @details
#' \code{dplindLnorm} computes the PLL density.
#' \code{pplindLnorm} computes the PLL CDF.
#' \code{qplindLnorm} computes quantiles.
#' \code{rplindLnorm} generates random draws.
#'
#' The PMF is:
#' \deqn{
#' f(y|\mu,\theta,\sigma)=\int_0^\infty
#'   \frac{\theta^2(\lambda x)^y(\theta+\lambda x+y+1)}
#'        {(\theta+1)(\theta+\lambda x)^{y+2}}
#'   \frac{\exp\left(-\frac{\ln^2(x)}{2\sigma^2}\right)}
#'        {x\sigma\sqrt{2\pi}}dx
#' }
#'
#' Where \eqn{\mu} is the \code{mean} argument (the Poisson-Lindley mean) and
#' \eqn{\lambda=\mu\theta(\theta+1)/(\theta+2)} is the corresponding
#' Poisson-Lindley scale. The mean of the distribution is:
#' \deqn{
#' E[y]=\mu e^{\sigma^2/2}=\frac{\lambda(\theta+2)e^{\sigma^2/2}}
#'             {\theta(\theta+1)}
#' }
#'
#' Halton draws are used to evaluate the integral.
#' 
#' @returns dplindLnorm gives the density, pplindLnorm gives the distribution 
#'  function, qplindLnorm gives the quantile function, and rplindLnorm generates
#'  random  deviates.
#' 
#'  The length of the result is determined by n for rplindLnorm, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' dplindLnorm(0, mean=0.75, theta=7, sigma=2, ndraws=10)
#' pplindLnorm(0:10, mean=0.75, theta=7, sigma=2, ndraws=10)
#' qplindLnorm(c(0.1,0.5,0.9), lambda=4.67, theta=7, sigma=2)
#' rplindLnorm(5, mean=0.75, theta=7, sigma=2)
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @name PoissonLindleyLognormal
#'
#' @rdname PoissonLindleyLognormal
#' @export
dplindLnorm <- function(x, mean=1, theta=1, sigma=1,
                        ndraws=1500, log=FALSE, hdraws=NULL,
                        engine = c("quadrature", "halton")) {
  # Draws given without an engine: the draw-based method
  if (missing(engine) && !is.null(hdraws)) engine <- "halton"
  engine <- match.arg(engine)
  if (engine == "quadrature") {
    lp <- .fcr_mix_logpmf(x, mean, theta, sigma, "PLL")
    return(if (log) lp else exp(lp))
  }
  
  if(!is.null(hdraws)) {
    h <- qnorm(hdraws)
  } else {
    h <- randtoolbox::halton(ndraws, normal=TRUE)
  }
  
  p <- dplindlogn_cpp(x, mean, theta, sigma, h)
  
  if(log) return(log(p))
  return(p)
}

#' @rdname PoissonLindleyLognormal
#' @export
pplindLnorm <- function(q, mean=1, theta=1, lambda=NULL,
                        sigma=1, ndraws=1500,
                        lower.tail=TRUE, log.p=FALSE,
                        engine = c("quadrature", "halton")) {
  engine <- match.arg(engine)
  if (!is.null(lambda)) {
    mean <- lambda * (theta + 2) / (theta * (theta + 1))
  }
  .fcr_p_contract(.pplindLnorm_core, q,
                  list(mean = mean, theta = theta, sigma = sigma),
                  lower.tail, log.p, ndraws = ndraws, engine = engine,
                  logdens = if (engine == "quadrature") function(x, ...)
                    dplindLnorm(x, ..., log = TRUE))
}

.pplindLnorm_core <- function(q, mean=1, theta=1, lambda=NULL,
                        sigma=1, ndraws=1500,
                        lower.tail=TRUE, log.p=FALSE,
                        engine = "quadrature") {
  
  if(!is.null(lambda)) {
    mean <- lambda * (theta + 2) / (theta * (theta + 1))
  }
  if (engine == "quadrature") {
    return(.fcr_mix_cdf(q, mean, theta, sigma, "PLL", lower.tail, log.p))
  }
  
  # FIX: Set normal=FALSE.
  # dplindLnorm applies qnorm() internally.
  h <- randtoolbox::halton(ndraws, normal=FALSE)
  
  n <- max(length(q), length(mean), length(theta), length(sigma))
  q <- rep_len(q, n)
  mean <- rep_len(mean, n)
  theta <- rep_len(theta, n)
  sigma <- rep_len(sigma, n)
  
  cdf <- numeric(n)
  
  for(i in seq_len(n)) {
    if(q[i] < 0) {
      cdf[i] <- 0
      next
    }
    y_seq <- 0:floor(q[i])
    probs <- dplindLnorm(
      y_seq, mean[i], theta[i], sigma[i],
      ndraws = ndraws, hdraws = h
    )
    cdf[i] <- sum(probs)
  }
  
  if(!lower.tail) cdf <- 1 - cdf
  if(log.p) return(log(cdf))
  return(cdf)
}

#' @rdname PoissonLindleyLognormal
#' @export
qplindLnorm <- function(p, mean = 1, theta = 1, sigma = 1, ndraws = 1500,
                        lambda = NULL, engine = c("quadrature", "halton")) {
  engine <- match.arg(engine)
  if (!is.null(lambda)) {
    mean <- lambda * (theta + 2) / (theta * (theta + 1))
  }
  if (any(mean <= 0 | theta <= 0 | sigma < 0, na.rm = TRUE)) {
    warning(paste("The values of `mean` (or `lambda`) and `theta` must be",
                  "greater than 0, and `sigma` must be non-negative."))
  }
  h <- randtoolbox::halton(ndraws)
  .fcr_qcount(p, function(y, mean, theta, sigma) {
    if (is.na(mean) || is.na(theta) || is.na(sigma)) {
      return(rep(NA_real_, length(y)))
    }
    # sigma = 0 is the Poisson-Lindley limit
    if (!(mean > 0 && theta > 0 && sigma >= 0)) return(rep(NaN, length(y)))
    cumsum(dplindLnorm(y, mean, theta, sigma, hdraws = h, engine = engine))
  }, mean, theta, sigma)
}

#' @rdname PoissonLindleyLognormal
#' @export
rplindLnorm <- function(n, mean = 1, theta = 1, sigma = 1,
                        ndraws = 1500, lambda = NULL) {
  n <- .fcr_rn(n)
  theta <- rep_len(theta, n); sigma <- rep_len(sigma, n)
  if (is.null(lambda)) {
    mean <- rep_len(mean, n)
    lambda <- mean * theta * (theta + 1) / (theta + 2)
  } else {
    lambda <- rep_len(lambda, n)
  }
  bad <- !(is.finite(lambda) & lambda > 0 & is.finite(theta) & theta > 0 &
             is.finite(sigma) & sigma >= 0)
  .fcr_warn_na(bad, paste("`mean` (or `lambda`) and `theta` must be > 0 and",
                          "`sigma` >= 0."))
  out <- rep(NA_real_, n); ok <- !bad; k <- sum(ok)
  out[ok] <- stats::rpois(k, lambda[ok] * .fcr_rlindley_mix(k, theta[ok]) *
                            exp(stats::rnorm(k, 0, sigma[ok])))
  out
}
