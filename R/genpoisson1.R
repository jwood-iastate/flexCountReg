# ------------------------------------------------------------
# Generalized Poisson helper functions
# ------------------------------------------------------------

.gp1_validate <-  function(mu, phi) {
  if (any(!is.finite(mu)) || any(mu <= 0, na.rm = TRUE)) {
    warning("`mu` must be greater than 0.")
  }
  if (any(!is.finite(phi)) || any(phi <= -1, na.rm = TRUE)) {
    warning("`phi` must be greater than -1.")
  }
}

.gp1_ymax <- function(mu, phi) {
  if (phi < 0) {
    floor(-mu / phi)
  } else {
    Inf
  }
}

.gp1_logpmf_raw <- function(y, mu, phi) {
  if (is.na(y) || is.na(mu) || is.na(phi)) return(NA_real_)
  
  if (y < 0 || abs(y - round(y)) > .Machine$double.eps^0.5) {
    return(-Inf)
  }
  
  if (mu <= 0 || phi <= -1) return(-Inf)
  
  if (phi == 0) {
    return(dpois(y, lambda = mu, log = TRUE))
  }
  
  if (phi < 0) {
    ymax <- .gp1_ymax(mu, phi)
    if (y > ymax) return(-Inf)
  }
  
  a <- mu + phi * y
  den <- 1 + phi
  
  if (a <= 0 || den <= 0) return(-Inf)
  
  log(mu) +
    (y - 1) * log(a) -
    a / den -
    y * log1p(phi) -
    lgamma(y + 1)
}

.gp1_log_norm_const <- function(mu, phi) {
  if (phi >= 0) return(0)
  
  ymax <- .gp1_ymax(mu, phi)
  y_seq <- 0:ymax
  logp <- vapply(y_seq, function(yy) .gp1_logpmf_raw(yy, mu, phi), numeric(1))
  
  m <- max(logp)
  m + log(sum(exp(logp - m)))
}

.gp1_logpmf <- function(y, mu, phi) {
  lp <- .gp1_logpmf_raw(y, mu, phi)
  if (is.infinite(lp) && lp < 0) return(lp)
  
  if (phi < 0) {
    lp - .gp1_log_norm_const(mu, phi)
  } else {
    lp
  }
}

.gp1_cdf_scalar <- function(q, mu, phi) {
  if (is.na(q) || is.na(mu) || is.na(phi)) return(NA_real_)
  if (q < 0) return(0)
  
  q <- floor(q)
  
  if (phi == 0) {
    return(ppois(q, lambda = mu))
  }
  
  if (phi < 0) {
    ymax <- .gp1_ymax(mu, phi)
    if (q >= ymax) return(1)
    
    y_seq <- 0:q
    logp <- vapply(y_seq, function(yy) .gp1_logpmf_raw(yy, mu, phi), numeric(1))
    logp <- logp - .gp1_log_norm_const(mu, phi)
  } else {
    y_seq <- 0:q
    logp <- vapply(y_seq, function(yy) .gp1_logpmf_raw(yy, mu, phi), numeric(1))
  }
  
  cum_log <- numeric(length(logp))
  cum_log[1] <- logp[1]
  
  if (length(logp) > 1) {
    for (k in 2:length(logp)) {
      a <- cum_log[k - 1]
      b <- logp[k]
      if (a >= b) {
        cum_log[k] <- a + log1p(exp(b - a))
      } else {
        cum_log[k] <- b + log1p(exp(a - b))
      }
    }
  }
  
  exp(cum_log[length(cum_log)])
}


#' Generalized Poisson Version 1 Distribution
#'
#' These functions provide the density function, distribution function,
#' quantile function, and random number generation for the
#' Generalized Poisson Version 1 (GP-1) Distribution
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mu numeric value or vector of mean values for the distribution (the
#' values have to be greater than 0).
#' @param phi single value or vector of values for the scale parameter of the
#' distribution (the values have to be greater than -1).
#' @param log logical; if TRUE, probabilities p are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#' otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dgp1} computes the density (PDF) of the Generalized Poisson
#' Distribution.
#'
#' \code{pgp1} computes the CDF of the Generalized Poisson Distribution.
#'
#' \code{qgp1} computes the quantile function of the
#' Generalized Poisson Distribution.
#'
#' \code{rgp1} generates random numbers from the Generalized Poisson
#' Distribution.
#'
#' The compound Probability Mass Function (PMF) for the Generalized Poisson
#' distribution, version 1, (GP-1) is:
#' \deqn{
#' f(y|\phi,\mu)=\frac{\mu(\mu+\phi y)^{y-1} exp\left(-\frac{\mu+\phi y}
#' {1+\phi}\right)}{(1+\phi)^y y!}
#' }
#'
#' Where \eqn{\phi} is a scale parameter with the restriction that
#' \eqn{\eta>0}, \eqn{\mu>0} is the mean value, and \eqn{y} is a non-negative
#' integer. This formulation uses the mean directly. 
#'
#' The variance of the GP-1 distribution is:
#' \deqn{\sigma^2=(1+\phi)^2 \mu}
#' 
#' If \eqn{\phi>0}, the distribution is overdispersed. If \eqn{\phi=0}, the 
#' distribution is equidispersed. If \eqn{\phi<0}, the distribution is 
#' underdispersed.
#' 
#' Furthermore, \eqn{phi>-1} is required for this distribution. When 
#' \eqn{\phi<0}, there is also a maximum value of support for the integar value 
#' \eqn{y}. This is \eqn{y_max=\left\lfloor -\frac{\mu}{\phi} \right\rfloor}.
#' 
#' @returns dgp1 gives the density, pgp1 gives the distribution 
#'  function, qgp1 gives the quantile function, and rgp1 generates random 
#'  deviates.
#' 
#' The length of the result is determined by n for rgp1, and is the 
#'  maximum of the lengths of the numerical arguments for the other functions.
#'  
#' @references  Consul, PoC, and Felix Famoye. "Generalized Poisson regression model." 
#'  Communications in Statistics-Theory and Methods 21.1 (1992): 89-109.
#'  
#'  Zamani, Hossein, and Noriszura Ismail. "Functional form for the generalized 
#'  Poisson regression model." Communications in Statistics-Theory and Methods 
#'  41.20 (2012): 3666-3675.
#'  
#' @examples
#' dgp1(1, mu=0.75, phi=-0.1)
#' pgp1(c(0,1,2,3,5,7,9,10), mu=0.75, phi=3)
#' qgp1(c(0.1,0.3,0.5,0.9,0.95), mu=0.75, phi=0.5)
#' rgp1(30, mu=0.75, phi=1.5)
#'
#' @importFrom stats runif
#' @name GeneralizedPoisson
#' @rdname GeneralizedPoisson
#' @export
dgp1 <- Vectorize(function(x, mu = 1, phi = 1, log = FALSE) {
  .gp1_validate(mu, phi)
  
  if (is.na(x) || is.na(mu) || is.na(phi)) {
    return(NA_real_)
  }
  
  if (x < 0 || abs(x - round(x)) > .Machine$double.eps^0.5) {
    return(if (log) -Inf else 0)
  }
  
  lp <- .gp1_logpmf(x, mu, phi)
  if (log) lp else exp(lp)
})

#' @rdname GeneralizedPoisson
#' @export
pgp1 <- function(q, mu = 1, phi = 1, lower.tail = TRUE, log.p = FALSE) {
  .gp1_validate(mu, phi)
  
  n <- max(length(q), length(mu), length(phi))
  q <- rep_len(q, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  
  out <- vapply(seq_len(n), function(i) {
    .gp1_cdf_scalar(q[i], mu[i], phi[i])
  }, numeric(1))
  
  if (!lower.tail) out <- 1 - out
  if (log.p) out <- log(out)
  
  out
}

#' @rdname GeneralizedPoisson
#' @export
qgp1 <- Vectorize(function(p, mu = 1, phi = 1, lower.tail = TRUE, log.p = FALSE) {
  .gp1_validate(mu, phi)
  
  if (is.na(p) || is.na(mu) || is.na(phi)) return(NA_real_)
  
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  
  if (p < 0 || p > 1) return(NaN)
  if (p == 0) return(0)
  
  if (phi == 0) {
    return(qpois(p, lambda = mu))
  }
  
  if (phi < 0) {
    ymax <- .gp1_ymax(mu, phi)
    y_seq <- 0:ymax
    
    logp <- vapply(y_seq, function(yy) .gp1_logpmf_raw(yy, mu, phi), numeric(1))
    logp <- logp - .gp1_log_norm_const(mu, phi)
    
    cdf <- numeric(length(logp))
    cdf[1] <- exp(logp[1])
    
    if (length(logp) > 1) {
      for (k in 2:length(logp)) {
        cdf[k] <- cdf[k - 1] + exp(logp[k])
      }
    }
    
    idx <- which(cdf >= p)[1]
    return(y_seq[idx])
  }
  
  if (p == 1) return(Inf)
  
  y <- 0L
  cdf <- .gp1_cdf_scalar(y, mu, phi)
  
  max_iter <- 1e6L
  iter <- 0L
  
  while (cdf < p) {
    y <- y + 1L
    cdf <- .gp1_cdf_scalar(y, mu, phi)
    iter <- iter + 1L
    if (iter > max_iter) {
      stop("Quantile search exceeded the iteration limit.")
    }
  }
  
  y
})

#' @rdname GeneralizedPoisson
#' @export
rgp1 <- function(n, mu = 1, phi = 1) {
  .gp1_validate(mu, phi)
  
  if (phi == 0) {
    return(rpois(n, lambda = mu))
  }
  
  if (phi < 0) {
    ymax <- .gp1_ymax(mu, phi)
    y_seq <- 0:ymax
    
    logp <- vapply(y_seq, function(yy) .gp1_logpmf_raw(yy, mu, phi), numeric(1))
    logp <- logp - .gp1_log_norm_const(mu, phi)
    prob <- exp(logp)
    
    return(sample(y_seq, size = n, replace = TRUE, prob = prob))
  }
  
  u <- runif(n)
  vapply(u, function(p) qgp1(p, mu = mu, phi = phi), numeric(1))
}
