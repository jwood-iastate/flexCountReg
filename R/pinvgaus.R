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
#' @param lower.tail logical; if TRUE, returns \eqn{P[X\leq x]}, otherwise P[X > x].
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
#' Type 2 modifies \eqn{\eta} → \eqn{\eta\mu}:
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
#' Dean, Lawless & Willmot (1989). Mixed Poisson–Inverse Gaussian Models.
#' Hilbe (2011). Negative Binomial Regression.
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
dpinvgaus <- Vectorize(function(x, mu=1, eta=1, form="Type 1", log=FALSE) {
  
  # integer check
  tst <- ifelse(
    is.na(nchar(strsplit(as.character(x), "\\.")[[1]][2]) > 0),
    FALSE, TRUE
  )
  
  if(tst || x < 0) warning(paste(
    "The value of `x` must be a non-negative",
    "whole number."
  ))
  
  if(eta <= 0) warning(paste(
    "The value of `eta` must be greater",
    "than 0."
  ))
  
  if(form == "Type 2") {
    eta <- eta * mu
  }
  
  p0 <- exp(mu / eta * (1 - sqrt(1 + 2 * eta)))
  
  if(x > 0) {
    j <- if(x > 1) seq(0, x - 1, 1) else 0
    e2 <- sum(
      gamma(x + j) /
        (gamma(x - j) * gamma(j + 1)) *
        (eta / (2 * mu))^j *
        (1 + 2 * eta)^(-j / 2)
    )
    p <- p0 * mu^x / gamma(x + 1) * (1 + 2 * eta)^(-x / 2) * e2
  } else {
    p <- p0
  }
  
  if(log) return(log(p))
  return(p)
})

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
  valid <- !invalid_q
  
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
qpinvgaus <- Vectorize(function(p, mu=1, eta=1, form="Type 1") {
  y <- 0
  p_value <- ppinvgaus(y, mu, eta, form)
  while(p_value < p) {
    y <- y + 1
    p_value <- ppinvgaus(y, mu, eta, form)
  }
  return(y)
})

#' @rdname PoissonInverseGaussian
#' @export
rpinvgaus <- function(n, mu=1, eta=1, form="Type 1") {
  u <- runif(n)
  y <- vapply(
    X = u, 
    FUN = \(p) qpinvgaus(p, mu, eta, form), 
    FUN.VALUE = numeric(1))
  return(y)
}

