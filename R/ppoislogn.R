#' Poisson-Lognormal Distribution
#'
#' These functions provide density, distribution function, quantile function,
#' and random number generation for the Poisson-Lognormal (PLogN) Distribution
#'
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mean numeric value or vector of mean values (not the expected value). 
#'  This is the mean for the Poisson portion of the distribution. Note that
#'  the mean values have to be greater than 0.
#' @param sigma single value or vector of values for the sigma parameter of the
#'   lognormal distribution (the values have to be greater than 0).
#' @param ndraws the number of Halton draws to use for the integration.
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#' @param hdraws and optional vector of Halton draws to use for the integration.
#' @param engine the engine to use for the integration. Options include "halton"
#'  which uses Halton draws or "poilog" which uses the poilog package.
#'
#' @details
#' \code{dpLnorm} computes the density (PDF) of the Poisson-Lognormal
#' Distribution.
#'
#' \code{ppLnorm} computes the CDF of the Poisson-Lognormal Distribution.
#'
#' \code{qpLnorm} computes the quantile function of the Poisson-Lognormal
#' Distribution.
#'
#' \code{rpLnorm} generates random numbers from the Poisson-Lognormal
#' Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Poisson-Lognormal
#' distribution is:
#' \deqn{f(y|\mu,\sigma)=
#'   \int_0^\infty \frac{\mu^y x^y e^{-\mu x}}{y!} 
#'   \frac{exp\left(-\frac{ln^2(x)}{2\sigma^2} \right)}{x\sigma\sqrt{2\pi}}dx}
#'
#' Where \eqn{\sigma} is a parameter for the lognormal distribution with the
#' restriction \eqn{\sigma>0}, and \eqn{y} is a non-negative integer.
#'
#' The expected value of the distribution is:
#' \deqn{E[y]=e^{X\beta+\sigma^2/2} = \mu e^{\sigma^2/2}}
#' 
#' The variance for the distribution is:
#' \deqn{V[Y]=E[Y]+\left(e^{\sigma^2}-1\right)E[Y]^2}
#' 
#' Halton draws are used to perform simulation over the lognormal distribution
#' to solve the integral if the engine is set to "halton". The poilog package is 
#' used to solve the integral if the engine is set to "poilog".
#' 
#' @returns dpLnorm gives the density, ppLnorm gives the distribution 
#'  function, qpLnorm gives the quantile function, and rpLnorm generates
#'  random  deviates.
#' 
#'  The length of the result is determined by n for rpLnorm, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' dpLnorm(0, mean=0.75, sigma=2, ndraws=10)
#' dpLnorm(0, mean=0.75, sigma=2, engine="poilog")
#' ppLnorm(c(0,1,2,3,5,7,9,10), mean=0.75, sigma=2, ndraws=10)
#' qpLnorm(c(0.1,0.3,0.5,0.9,0.95), mean=0.75, sigma=2, ndraws=10)
#' rpLnorm(30, mean=0.75,  sigma=2)
#'
#' @importFrom stats runif qnorm
#' @importFrom randtoolbox halton
#' @export
#' @name PoissonLognormal
#' 
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @rdname PoissonLognormal
#' @export
dpLnorm <- function(x, mean = 1, sigma = 1, ndraws = 1500, log = FALSE,
                    hdraws = NULL, engine = c("halton", "poilog")) {
  engine <- match.arg(engine)
  
  if (any(mean <= 0, na.rm = TRUE)) {
    stop("The values of `mean` must be greater than 0.")
  }
  if (any(sigma <= 0, na.rm = TRUE)) {
    stop("The values of `sigma` must be greater than 0.")
  }
  
  r <- .fcr_recycle(x, mean, sigma)
  if (length(r[[1]]) == 0L) return(numeric(0))
  
  if (engine == "halton") {
    h <- if (!is.null(hdraws)) {
      stats::qnorm(hdraws)
    } else {
      randtoolbox::halton(n = ndraws, dim = 1, normal = TRUE)
    }
    # The C++ kernel is vectorized over x, mean and sigma
    p <- dpLnorm_cpp(x = as.numeric(r[[1]]), mean = r[[2]], sigma = r[[3]],
                     h = as.numeric(h))
  } else {
    rlang::check_installed("poilog", reason = "for `engine = \"poilog\"`.")
    p <- mapply(function(xi, mi, si) poilog::dpoilog(xi, base::log(mi), si),
                r[[1]], r[[2]], r[[3]], USE.NAMES = FALSE)
  }
  
  if (log) base::log(p) else p
}
#' @rdname PoissonLognormal
#' @export
ppLnorm <- function(q, mean = 1, sigma = 1, ndraws = 1500,
                    lower.tail = TRUE, log.p = FALSE,
                    engine = c("halton", "poilog")) {
  engine <- match.arg(engine)
  if (any(mean <= 0, na.rm = TRUE) || any(sigma <= 0, na.rm = TRUE)) {
    warning("The values of `mean` and `sigma` must be greater than 0.")
  }
  r <- .fcr_recycle(q, mean, sigma)
  q <- r[[1]]; mean <- r[[2]]; sigma <- r[[3]]
  # One set of Halton draws for every element
  h <- if (engine == "halton") randtoolbox::halton(ndraws) else NULL
  cdf <- rep(NA_real_, length(q))
  cdf[!is.na(q) & q < 0] <- 0
  for (i in which(!is.na(q) & q >= 0 & mean > 0 & sigma > 0)) {
    cdf[i] <- sum(dpLnorm(0:floor(q[i]), mean[i], sigma[i], hdraws = h,
                          engine = engine))
  }
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) base::log(cdf) else cdf
}

#' @rdname PoissonLognormal
#' @export
qpLnorm <- function(p, mean = 1, sigma = 1, ndraws = 1500,
                    engine = c("halton", "poilog")) {
  engine <- match.arg(engine)
  if (any(mean <= 0, na.rm = TRUE) || any(sigma <= 0, na.rm = TRUE)) {
    warning("The values of `mean` and `sigma` must be greater than 0.")
  }
  h <- if (engine == "halton") randtoolbox::halton(ndraws) else NULL
  .fcr_qcount(p, function(y, mean, sigma) {
    cumsum(dpLnorm(y, mean, sigma, hdraws = h, engine = engine))
  }, mean, sigma)
}


#' @rdname PoissonLognormal
#' @export
rpLnorm <- function(n, mean = 1, sigma = 1) {
  n <- .fcr_rn(n)
  mean <- rep_len(mean, n); sigma <- rep_len(sigma, n)
  bad <- !(is.finite(mean) & mean > 0 & is.finite(sigma) & sigma > 0)
  .fcr_warn_na(bad, "`mean` and `sigma` must be greater than 0.")
  out <- rep(NA_real_, n); ok <- !bad
  out[ok] <- stats::rpois(sum(ok),
                          mean[ok] * exp(stats::rnorm(sum(ok), 0, sigma[ok])))
  out
}
