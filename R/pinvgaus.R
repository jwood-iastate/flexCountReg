#' Generate inverse-Gaussian random numbers
#'
#' Generates random variables from an inverse-Gaussian distribution using
#' the Michael, Schucany, and Haas (1976) method.
#'
#' Parameterization:
#' IG(mean = mu, shape = lambda)
#'
#' @param n Number of observations.
#' @param mean Mean parameter (> 0).
#' @param shape Shape parameter (> 0).
#' @returns A numeric vector of length n.
#' @keywords internal
rinvgaus <- function(n, mean = 1, shape = 1) {
  
  if(length(n) != 1L || is.na(n) || n < 0) {
    stop("`n` must be a non-negative integer.")
  }
  
  n <- as.integer(n)
  
  if(n == 0L) {
    return(numeric(0))
  }
  
  if(any(mean <= 0, na.rm = TRUE)) {
    stop("`mean` must be positive.")
  }
  
  if(any(shape <= 0, na.rm = TRUE)) {
    stop("`shape` must be positive.")
  }
  
  mean <- rep_len(mean, n)
  shape <- rep_len(shape, n)
  
  z <- rnorm(n)
  v <- z^2
  m <- mean
  lambda <- shape
  
  x <- m +
    (m^2 * v) / (2 * lambda) -
    (m / (2 * lambda)) *
    sqrt(4 * m * lambda * v + m^2 * v^2)
  
  u <- runif(n)
  
  ifelse(
    u <= m / (m + x),
    x,
    m^2 / x
  )
}

#' Poisson-Inverse-Gaussian Distribution
#'
#' These functions provide the density function, distribution function,
#' quantile function, and random number generation for the
#' Poisson-Inverse-Gaussian (PInvGaus) Distribution.
#'
#' The Poisson-Inverse-Gaussian distribution is a special case of the
#' Sichel distribution (Cameron & Trivedi, 2013). It is also known as a
#' univariate Sichel distribution (Hilbe, 2011).
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mu numeric value or vector of mean values for the distribution
#'   (values must be > 0).
#' @param eta single value or vector of scale parameter values (must be > 0).
#' @param form optional parameter specifying the formulation. Options:
#'   "Type 1" (standard) or "Type 2" (Dean et al., 1987).
#' @param log logical; if TRUE, probabilities p are returned as log(p).
#' @param log.p logical; if TRUE, probabilities p are returned as log(p).
#' @param lower.tail logical; if TRUE, returns \eqn{P[X\leq x]}, otherwise 
#' \eqn{P[X > x]}.
#'
#' @details
#' \code{dpinvgaus} computes the PDF of the Poisson-Inverse-Gaussian dist.
#'
#' \code{ppinvgaus} computes the CDF of the Poisson-Inverse-Gaussian dist.
#'
#' \code{qpinvgaus} computes quantiles of the Poisson-Inverse-Gaussian dist.
#'
#' \code{rpinvgaus} generates random numbers from the distribution.
#'
#' The PMF (Type 1) is:
#' \deqn{
#' f(y|\eta,\mu)=\begin{cases}
#' f(0)=\exp\left(\frac{\mu}{\eta}(1-\sqrt{1+2\eta})\right)\\
#' f(y>0)=f(0)\frac{\mu^y}{y!}(1+2\eta)^{-y/2}
#'     \sum_{j=0}^{y-1}\frac{\Gamma(y+j)}
#'     {\Gamma(y-j)\Gamma(j+1)}
#'     \left(\frac{\eta}{2\mu}\right)^j(1+2\eta)^{-j/2}
#' \end{cases}}
#'
#' The variance is:
#' \deqn{\sigma^2=\mu+\eta\mu}
#'
#' Type 2 modifies \eqn{\eta} -> \eqn{\eta\mu}:
#' \deqn{
#' f(0)=\exp\left(\frac{1}{\eta}(1-\sqrt{1+2\eta\mu})\right)
#' }
#' \deqn{
#' f(y>0)=f(0)\frac{\mu^y}{y!}(1+2\eta\mu)^{-y/2}
#'   \sum_{j=0}^{y-1}\frac{\Gamma(y+j)}
#'   {\Gamma(y-j)\Gamma(j+1)}
#'   \left(\frac{\eta}{2}\right)^j(1+2\eta\mu)^{-j/2}
#' }
#'
#' Resulting variance:
#' \deqn{\sigma^2=\mu+\eta\mu^2}
#'
#' @references
#' Cameron & Trivedi (2013). Regression Analysis of Count Data.
#' Dean, Lawless & Willmot (1989). Mixed Poisson-Inverse Gaussian Models.
#' Hilbe (2011). Negative Binomial Regression.
#' 
#' @returns description dpinvgaus gives the density, ppinvgaus gives the distribution 
#'  function, qpinvgaus gives the quantile function, and rpinvgaus generates
#'  random  deviates.
#' 
#'  The length of the result is determined by n for rpinvgaus, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'
#' @examples
#' dpinvgaus(1, mu=0.75, eta=1)
#' ppinvgaus(c(0,1,2,3,5,7,9,10), mu=0.75, eta=3, form="Type 2")
#' qpinvgaus(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, eta=0.5, form="Type 2")
#' rpinvgaus(30, mu=0.75, eta=1.5)
#'
#' @importFrom stats runif
#' @export
#' @name PoissonInverseGaussian

#' @rdname PoissonInverseGaussian
#' @export
dpinvgaus <- function(x, mu = 1, eta = 1, form = "Type 1", log = FALSE) {
  form <- match.arg(form, c("Type 1", "Type 2"))
  r <- .fcr_recycle(x, mu, eta)
  x <- r[[1]]; mu <- r[[2]]; eta <- r[[3]]
  
  is_count <- .fcr_is_count(x)
  if (any(!is_count & !is.na(x))) {
    warning("The value of `x` must be a non-negative whole number.")
  }
  bad_par <- !is.na(mu) & !is.na(eta) & (mu <= 0 | eta <= 0)
  if (any(bad_par)) {
    warning("The values of `mu` and `eta` must be greater than 0.")
  }
  
  lp <- suppressWarnings(
    .ldpinvgaus(ifelse(is_count, round(x), 0), mu, eta, form))
  lp[!is_count] <- -Inf
  lp[bad_par] <- NaN
  lp[is.na(x) | is.na(mu) | is.na(eta)] <- NA
  if (log) lp else exp(lp)
}

# Log PMF of the Poisson-inverse-Gaussian, computed with lgamma() and
# log-sum-exp so it stays finite for large counts.
.ldpinvgaus <- function(x, mu, eta, form = "Type 1") {
  if (form == "Type 2") eta <- eta * mu
  l1 <- log1p(2 * eta)
  out <- mu / eta * (1 - sqrt(1 + 2 * eta))          # log P(0)
  for (i in which(is.finite(x) & x > 0 & mu > 0 & eta > 0)) {
    j <- 0:(x[i] - 1)
    terms <- lgamma(x[i] + j) - lgamma(x[i] - j) - lgamma(j + 1) +
      j * log(eta[i] / (2 * mu[i])) - j / 2 * l1[i]
    m <- max(terms)
    out[i] <- out[i] + x[i] * log(mu[i]) - lgamma(x[i] + 1) -
      x[i] / 2 * l1[i] + m + log(sum(exp(terms - m)))
  }
  out
}

#' @rdname PoissonInverseGaussian
#' @export
ppinvgaus <- function(q, mu=1, eta=1, form="Type 1",
                      lower.tail=TRUE, log.p=FALSE) {
  
  if(any(eta <= 0, na.rm=TRUE)) warning("'eta' must be positive")
  if(any(mu <= 0, na.rm=TRUE)) warning("'mu' must be positive")
  
  form <- match.arg(form, c("Type 1", "Type 2"))
  
  n <- max(length(q), length(mu), length(eta))
  q <- rep_len(as.integer(floor(q)), n)
  mu <- rep_len(mu, n)
  eta <- rep_len(eta, n)
  
  cdf <- rep(NA_real_, n)
  invalid_q <- is.na(q) | q < 0
  cdf[invalid_q] <- 0
  bad_par <- !(is.finite(mu) & mu > 0 & is.finite(eta) & eta > 0)
  cdf[bad_par & !invalid_q] <- NaN
  valid <- !invalid_q & !bad_par
  
  if(any(valid)) {
    for(i in which(valid)) {
      mu_i <- mu[i]
      eta_i <- eta[i]
      q_i <- q[i]
      
      eta_adj <- if(form == "Type 2") eta_i * mu_i else eta_i
      
      log_pmf <- numeric(q_i + 1)
      sqrt_term <- sqrt(1 + 2 * eta_adj)
      log_p0 <- mu_i / eta_adj * (1 - sqrt_term)
      log_pmf[1] <- log_p0
      
      if(q_i > 0) {
        log_1p2e <- log(1 + 2 * eta_adj)
        
        for(y in 1:q_i) {
          j_seq <- 0:(y - 1)
          
          log_sum_terms <- lgamma(y + j_seq) -
            lgamma(y - j_seq) - lgamma(j_seq + 1) +
            j_seq * log(eta_adj / (2 * mu_i)) -
            j_seq / 2 * log_1p2e
          
          max_log <- max(log_sum_terms[is.finite(log_sum_terms)])
          log_sum <- if(is.finite(max_log)) {
            max_log + log(sum(exp(log_sum_terms - max_log)))
          } else -Inf
          
          log_pmf[y + 1] <- log_p0 +
            y * log(mu_i) -
            lgamma(y + 1) -
            y / 2 * log_1p2e +
            log_sum
        }
      }
      
      finite_mask <- is.finite(log_pmf)
      if(!any(finite_mask)) {
        cdf[i] <- NA_real_
      } else {
        max_log <- max(log_pmf[finite_mask])
        log_cdf <- max_log +
          log(sum(exp(log_pmf[finite_mask] - max_log)))
        cdf[i] <- exp(log_cdf)
      }
    }
  }
  
  if(!lower.tail) cdf <- 1 - cdf
  if(log.p) cdf <- log(cdf)
  
  return(cdf)
}

#' @rdname PoissonInverseGaussian
#' @export
qpinvgaus <- function(p, mu = 1, eta = 1, form = "Type 1") {
  form <- match.arg(form, c("Type 1", "Type 2"))
  .fcr_qcount(p, function(y, mu, eta) {
    cumsum(dpinvgaus(y, mu = mu, eta = eta, form = form))
  }, mu, eta)
}

#' @rdname PoissonInverseGaussian
#' @export
rpinvgaus <- function(n, mu = 1, eta = 1, form = "Type 1") {
  
  if(length(n) != 1L || is.na(n) || n < 0) {
    stop("`n` must be a non-negative integer.")
  }
  
  n <- as.integer(n)
  
  if(any(mu <= 0, na.rm = TRUE)) {
    stop("`mu` must be positive.")
  }
  
  if(any(eta <= 0, na.rm = TRUE)) {
    stop("`eta` must be positive.")
  }
  
  form <- match.arg(form, c("Type 1", "Type 2"))
  
  mu <- rep_len(mu, n)
  eta <- rep_len(eta, n)
  
  # Mean-one inverse-Gaussian multiplier
  #
  # Type 1:
  #   W ~ IG(1, mu / eta)
  #
  # Type 2:
  #   W ~ IG(1, 1 / eta)
  shape <- if(form == "Type 1") {
    mu / eta
  } else {
    1 / eta
  }
  
  w <- rinvgaus(
    n = n,
    mean = 1,
    shape = shape
  )
  
  # Random Poisson mean
  lambda <- mu * w
  
  # Poisson mixture
  rpois(n, lambda = lambda)
}

