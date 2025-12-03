#' Poisson-Inverse-Gaussian Distribution
#'
#' These functions provide the density function, distribution function, quantile function, and random number generation for the Poisson-Inverse-Gaussian (PInvGaus) Distribution
#'
#' The Poisson-Inverse-Gaussian distribution is a special case of the Sichel
#' distribution, as noted by Cameron & Trivedi (2013). It is also known as a
#' univariate Sichel distribution (Hilbe, 2011).
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mu numeric value or vector of mean values for the distribution (the values have to be greater than 0).
#' @param eta single value or vector of values for the scale parameter of the distribution (the values have to be greater than 0).
#' @param form optional parameter indicating which formulation to use. Options include "Type 1" which is the standard form or "Type 2" which follows the formulation by Dean et. al. (1987).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]} otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dpinvgaus} computes the density (PDF) of the Poisson-Inverse-Gaussian Distribution.
#'
#' \code{ppinvgaus} computes the CDF of the Poisson-Inverse-Gaussian Distribution.
#'
#' \code{qpinvgaus} computes the quantile function of the Poisson-Inverse-Gaussian Distribution.
#'
#' \code{rpinvgaus} generates random numbers from the Poisson-Inverse-Gamma Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Poisson-Inverse-Gaussian distribution (Type 1) is (Cameron & Trivedi, 2013):
#' \deqn{f(y|\eta,\mu)=\begin{cases}
#'                              f(y=0)=\exp\left(\frac{\mu}{\eta}\left(1-\sqrt{1+2\eta}\right)\right) \\
#'                              f(y|y>0)=f(y=0)\frac{\mu^y}{y!}(1+2\eta)^{-y/2}\cdot\sum_{j=0}^{y-1}\frac{\Gamma(y+j)}{\Gamma(y-j)\Gamma(j+1)}\left(\frac{\eta}{2\mu}\right)^2(1+2\eta)^{-j/2}
#'                              \end{cases}}
#' 
#' Where \eqn{\eta} is a scale parameter with the restriction that \eqn{\eta>0}, \eqn{\mu} is the mean value, and \eqn{y} is a non-negative integer.
#'
#' The variance of the distribution is:
#' \deqn{\sigma^2=\mu+\eta\mu}
#' 
#' The alternative parameterization by Dean et. al. (1987) replaces \eqn{\eta} with \eqn{\eta\mu}. This version (Type 2) has the PMF:
#' \deqn{f(y|\eta,\mu)=\begin{cases}
#'                              f(y=0)=\exp\left(\frac{1}{\eta}\left(1-\sqrt{1+2\eta\mu}\right)\right) \\
#'                              f(y|y>0)=f(y=0)\frac{\mu^y}{y!}(1+2\eta\mu)^{-y/2}\cdot\sum_{j=0}^{y-1}\frac{\Gamma(y+j)}{\Gamma(y-j)\Gamma(j+1)}\left(\frac{\eta}{2}\right)^2(1+2\eta\mu)^{-j/2}
#'                              \end{cases}}
#' 
#'  This results in the variance of:
#' \deqn{\sigma^2=\mu+\eta\mu^2}
#'
#' @references 
#' Cameron, A. C., & Trivedi, P. K. (2013). Regression analysis of count data, 2nd Edition. Cambridge university press.
#' 
#' Dean, C., Lawless, J. F., & Willmot, G. E. (1989). A mixed Poisson–Inverse‐Gaussian regression model. Canadian Journal of Statistics, 17(2), 171-181.
#' 
#' Hilbe, J. M. (2011). Negative binomial regression. Cambridge University Press.
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
dpinvgaus <- Vectorize(function(x, mu=1, eta = 1, form="Type 1", log=FALSE){
  #test to make sure the value of x is an integer
  tst <- ifelse(is.na(nchar(strsplit(as.character(x), "\\.")[[1]][2])>0),FALSE, TRUE)
  if(tst || x < 0) warning("The value of `x` must be a non-negative whole number")
  if(eta<=0) warning("The value of `eta` must be greater than 0.")
  if(form=='Type 2'){
    eta <- eta*mu # make adjusted value for Type 2 formulation
  }
  p0 <- exp(mu/eta*(1-sqrt(1+2*eta)))
  
  if(x>0){
    if(x>1){
      j <- seq(0,(x-1),1)
    }
    else{
      j <- 0
    }
    e2 <- sum(gamma(x+j)/(gamma(x-j)*gamma(j+1))*(eta/(2*mu))^j *(1+2*eta)^(-j/2))
    p <- p0 * (mu^x)/gamma(x+1)*(1+2*eta)^(-x/2) * e2
  }
  else{
    p <- p0
  }
  if (log) return(log(p))
  else return(p)
})

#' @rdname PoissonInverseGaussian
#' @export
ppinvgaus <- function(q, mu = 1, eta = 1, form = "Type 1", 
                      lower.tail = TRUE, log.p = FALSE) {
  
  
  # --- Input Validation ---
  
  if (any(eta <= 0, na.rm = TRUE)) warning("'eta' must be positive")
  if (any(mu <= 0, na.rm = TRUE)) warning("'mu' must be positive")
  form <- match.arg(form, c("Type 1", "Type 2"))
  
  
  # --- Vectorization Setup ---
  
  n <- max(length(q), length(mu), length(eta))
  q <- rep_len(as.integer(floor(q)), n)
  mu <- rep_len(mu, n)
  eta <- rep_len(eta, n)
  
  
  # --- Initialize Result ---
  
  cdf <- rep(NA_real_, n)
  
  
  # --- Handle Special Cases ---
  
  invalid_q <- is.na(q) | q < 0
  cdf[invalid_q] <- 0
  
  
  # --- Compute CDF for Valid Cases ---
  
  valid <- !invalid_q
  
  if (any(valid)) {
    for (i in which(valid)) {
      mu_i <- mu[i]
      eta_i <- eta[i]
      q_i <- q[i]
      
      # Adjust eta for Type 2 formulation
      eta_adj <- if (form == "Type 2") eta_i * mu_i else eta_i
      
      # Compute log-PMF for 0:q_i
      log_pmf <- numeric(q_i + 1)
      
      # P(Y=0) = exp(mu/eta * (1 - sqrt(1 + 2*eta)))
      sqrt_term <- sqrt(1 + 2 * eta_adj)
      log_p0 <- mu_i / eta_adj * (1 - sqrt_term)
      
      log_pmf[1] <- log_p0  # y = 0
      
      if (q_i > 0) {
        # For y > 0:
        # P(y) = P(0) * (mu^y / y!) * (1 + 2*eta)^(-y/2) * sum_j(...)
        
        log_1_plus_2eta <- log(1 + 2 * eta_adj)
        
        for (y in 1:q_i) {
          # Compute the summation term
          j_seq <- 0:(y - 1)
          
          # log of each term in sum:
          # log(gamma(y+j)) - log(gamma(y-j)) - log(gamma(j+1)) + j*log(eta/(2*mu)) - j/2*log(1+2*eta)
          log_sum_terms <- lgamma(y + j_seq) - lgamma(y - j_seq) - lgamma(j_seq + 1) +
            j_seq * log(eta_adj / (2 * mu_i)) -
            j_seq / 2 * log_1_plus_2eta
          
          # Log-sum-exp for the summation
          max_log_term <- max(log_sum_terms[is.finite(log_sum_terms)])
          if (!is.finite(max_log_term)) {
            log_sum <- -Inf
          } else {
            log_sum <- max_log_term + log(sum(exp(log_sum_terms - max_log_term)))
          }
          
          # Full log-PMF for y
          # log(P(y)) = log(P(0)) + y*log(mu) - lgamma(y+1) - y/2*log(1+2*eta) + log(sum)
          log_pmf[y + 1] <- log_p0 + y * log(mu_i) - lgamma(y + 1) -
            y / 2 * log_1_plus_2eta + log_sum
        }
      }
      
      # Log-sum-exp for CDF
      finite_mask <- is.finite(log_pmf)
      if (!any(finite_mask)) {
        cdf[i] <- NA_real_
      } else {
        max_log <- max(log_pmf[finite_mask])
        log_cdf <- max_log + log(sum(exp(log_pmf[finite_mask] - max_log)))
        cdf[i] <- exp(log_cdf)
      }
    }
  }
  
  
  # --- Apply Tail and Log Options ---
  
  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  
  return(cdf)
}

#' @rdname PoissonInverseGaussian
#' @export
qpinvgaus <- Vectorize(function(p, mu=1, eta = 1, form="Type 1") {
  y <- 0
  p_value <- ppinvgaus(y, mu, eta, form)
  while(p_value < p){
    y <- y + 1
    p_value <- ppinvgaus(y, mu, eta, form)
  }
  return(y)
})


#' @rdname PoissonInverseGaussian
#' @export
rpinvgaus <- function(n, mu=1, eta = 1, form="Type 1") {
  u <- runif(n)
  y <- sapply(u, function(p) qpinvgaus(p, mu, eta, form))
  return(y)
}
