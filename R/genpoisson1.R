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

# Upper limit of the support used for underdispersed (phi < 0) GP models.
# The exact support ends at floor(-mu / phi), which grows without bound as
# phi -> 0-. The distribution is then underdispersed (variance below mu), so
# the probability beyond mu + 40 sqrt(mu) + 40 is negligible (far below
# double precision) and the support is truncated there.
.gp_cap <- function(mu) ceiling(mu + 40 * sqrt(mu) + 40)

.gp1_ymax <- function(mu, phi) {
  if (phi < 0) {
    min(floor(-mu / phi), .gp_cap(mu))
  } else {
    Inf
  }
}

# Vectorized unnormalized GP1 log-PMF (GP2 is GP1 with phi = alpha * mu).
.gp1_lraw_vec <- function(y, mu, phi) {
  a <- mu + phi * y
  lp <- suppressWarnings(log(mu) + (y - 1) * log(a) - a / (1 + phi) -
                           y * log1p(phi) - lgamma(y + 1))
  pois <- !is.na(phi) & phi == 0
  lp[pois] <- stats::dpois(y[pois], mu[pois], log = TRUE)
  bad <- !pois & (is.na(a) | a <= 0)
  lp[bad] <- -Inf
  lp[!is.na(mu) & !is.na(phi) & (mu <= 0 | phi <= -1)] <- -Inf
  lp[is.na(y) | is.na(mu) | is.na(phi)] <- NA
  lp
}

# Vectorized log normalizing constant (0 unless phi < 0).
.gp1_lnorm_vec <- function(mu, phi) {
  out <- numeric(length(mu))
  for (i in which(!is.na(phi) & phi < 0 & phi > -1 & !is.na(mu) & mu > 0)) {
    K <- min(floor(-mu[i] / phi[i]), .gp_cap(mu[i]))
    lp <- .gp1_lraw_vec(0:K, rep(mu[i], K + 1), rep(phi[i], K + 1))
    m <- max(lp)
    out[i] <- m + log(sum(exp(lp - m)))
  }
  out
}

# Vectorized GP1 log-PMF for counts y (support beyond the truncation is -Inf).
.gp1_logpmf_vec <- function(y, mu, phi) {
  lp <- .gp1_lraw_vec(y, mu, phi)
  neg <- !is.na(phi) & phi < 0
  if (any(neg)) {
    cap <- ifelse(neg, pmin(floor(-mu / phi), .gp_cap(mu)), Inf)
    lp[neg & y > cap] <- -Inf
    lp[neg] <- lp[neg] - .gp1_lnorm_vec(mu[neg], phi[neg])
  }
  lp
}

.gp1_log_norm_const <- function(mu, phi) {
  if (phi >= 0) return(0)
  .gp1_lnorm_vec(mu, phi)
}

# Unnormalized log-PMF over the counts y for scalar mu and phi.
.gp1_lraw <- function(y, mu, phi) {
  .gp1_lraw_vec(y, rep_len(mu, length(y)), rep_len(phi, length(y)))
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
    logp <- .gp1_lraw(y_seq, mu, phi)
    logp <- logp - .gp1_log_norm_const(mu, phi)
  } else {
    y_seq <- 0:q
    logp <- .gp1_lraw(y_seq, mu, phi)
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
#' Where \eqn{\phi} is a scale parameter with the restrictions that
#' \eqn{\mu>0} is the mean value, \deqn{\phi>-1}, and \eqn{y} is a non-negative
#' integer. This formulation uses the mean directly. 
#'
#' The variance of the GP-1 distribution is:
#' \deqn{\sigma^2=(1+\phi)^2 \mu}
#' 
#' If \eqn{\phi>0}, the distribution is overdispersed. If \eqn{\phi=0}, the 
#' distribution is equidispersed. If \eqn{\phi<0}, the distribution is 
#' underdispersed.
#' 
#' Furthermore, \eqn{\phi>-1} is required for this distribution. When 
#' \eqn{\phi<0}, there is also a maximum value of support for the integer value 
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
dgp1 <- function(x, mu = 1, phi = 1, log = FALSE) {
  .gp1_validate(mu, phi)
  r <- .fcr_recycle(x, mu, phi)
  x <- r[[1]]; mu <- r[[2]]; phi <- r[[3]]
  is_count <- .fcr_is_count(x)
  lp <- rep(-Inf, length(x))
  lp[is_count] <- .gp1_logpmf_vec(round(x[is_count]), mu[is_count],
                                  phi[is_count])
  lp[!is.na(mu) & !is.na(phi) & (mu <= 0 | phi <= -1 | !is.finite(mu) |
                                    !is.finite(phi))] <- NaN
  lp[is.na(x) | is.na(mu) | is.na(phi)] <- NA
  if (log) lp else exp(lp)
}

#' @rdname GeneralizedPoisson
#' @export
pgp1 <- function(q, mu = 1, phi = 1, lower.tail = TRUE, log.p = FALSE) {
  .gp1_validate(mu, phi)
  
  n <- max(length(q), length(mu), length(phi))
  q <- rep_len(q, n)
  mu <- rep_len(mu, n)
  phi <- rep_len(phi, n)
  
  out <- vapply(seq_len(n), function(i) {
    if (!is.na(mu[i]) && !is.na(phi[i]) &&
        (!is.finite(mu[i]) || mu[i] <= 0 || !is.finite(phi[i]) ||
         phi[i] <= -1)) {
      return(NaN)
    }
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
  if (!is.finite(mu) || mu <= 0 || !is.finite(phi) || phi <= -1) return(NaN)
  
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
    
    logp <- .gp1_lraw(y_seq, mu, phi)
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
  n <- .fcr_rn(n)
  mu <- rep_len(mu, n); phi <- rep_len(phi, n)
  ok <- is.finite(mu) & mu > 0 & is.finite(phi) & phi > -1
  .fcr_warn_na(!ok, "`mu` must be > 0 and `phi` > -1.")
  .fcr_r_grouped(n, list(mu = mu, phi = phi), ok, .rgp1_one)
}

.rgp1_one <- function(n, mu = 1, phi = 1) {
  
  if (phi == 0) {
    return(rpois(n, lambda = mu))
  }
  
  if (phi < 0) {
    ymax <- .gp1_ymax(mu, phi)
    y_seq <- 0:ymax
    
    logp <- .gp1_lraw(y_seq, mu, phi)
    logp <- logp - .gp1_log_norm_const(mu, phi)
    prob <- exp(logp)
    
    return(sample(y_seq, size = n, replace = TRUE, prob = prob))
  }
  
  u <- runif(n)
  vapply(u, function(p) qgp1(p, mu = mu, phi = phi), numeric(1))
}
