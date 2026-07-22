# ------------------------------------------------------------
# Generalized Poisson helper functions
# ------------------------------------------------------------

.gp2_validate <- function(mu, alpha) {
  ok <- is.finite(mu) & is.finite(alpha) & (mu > 0) & (1 + alpha * mu > 0)
  ok[is.na(ok)] <- FALSE
  
  if (any(!ok)) {
    warning(
      "Invalid Generalized Poisson-2 parameter values. Valid values require ",
      "mu > 0 and 1 + alpha * mu > 0.",
      call. = FALSE
    )
  }
  
  ok
}

.gp2_ymax <- function(alpha) {
  if (alpha < 0) {
    ceiling(-1 / alpha) - 1L
  } else {
    Inf
  }
}

.gp2_logpmf_raw <- function(y, mu, alpha) {
  if (is.na(y) || is.na(mu) || is.na(alpha)) return(NA_real_)
  
  if (y < 0 || abs(y - round(y)) > .Machine$double.eps^0.5) {
    return(-Inf)
  }
  
  if (mu <= 0 || !is.finite(alpha)) return(-Inf)
  
  if (alpha == 0) {
    return(dpois(y, lambda = mu, log = TRUE))
  }
  
  a <- mu + alpha * mu * y
  den <- 1 + alpha * mu
  
  if (a <= 0 || den <= 0) return(-Inf)
  
  log(mu) +
    (y - 1) * log(a) -
    a / den -
    y * log1p(alpha * mu) -
    lgamma(y + 1)
}

.gp2_log_norm_const <- function(mu, alpha) {
  if (alpha >= 0) return(0)
  
  ymax <- .gp2_ymax(alpha)
  y_seq <- 0:ymax
  logp <- vapply(y_seq, function(yy) .gp2_logpmf_raw(yy, mu, alpha), numeric(1))
  
  finite <- is.finite(logp)
  if (!any(finite)) return(0)
  
  m <- max(logp[finite])
  m + log(sum(exp(logp[finite] - m)))
}

.gp2_logpmf <- function(y, mu, alpha) {
  lp <- .gp2_logpmf_raw(y, mu, alpha)
  if (is.na(lp) || (is.infinite(lp) && lp < 0)) return(lp)
  
  if (alpha < 0) {
    lp - .gp2_log_norm_const(mu, alpha)
  } else {
    lp
  }
}

.gp2_cdf_scalar <- function(q, mu, alpha) {
  if (is.na(q) || is.na(mu) || is.na(alpha)) return(NA_real_)
  if (q < 0) return(0)
  
  q <- floor(q)
  
  if (alpha == 0) {
    return(ppois(q, lambda = mu))
  }
  
  if (alpha < 0) {
    ymax <- .gp2_ymax(alpha)
    if (q >= ymax) return(1)
  }
  
  y_seq <- 0:q
  logp <- vapply(y_seq, function(yy) .gp2_logpmf_raw(yy, mu, alpha), numeric(1))
  
  if (alpha < 0) {
    logp <- logp - .gp2_log_norm_const(mu, alpha)
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

#' Generalized Poisson Version 2 Distribution
#'
#' These functions provide the density function, distribution function,
#' quantile function, and random number generation for the
#' Generalized Poisson Version 2 (GP-2) Distribution.
#'
#' @param x numeric value or a vector of values.
#' @param q quantile or a vector of quantiles.
#' @param p probability or a vector of probabilities.
#' @param n the number of random numbers to generate.
#' @param mu numeric value or vector of mean values for the distribution (the
#'   values have to be greater than 0).
#' @param alpha numeric value or vector of values for the dispersion parameter
#'   of the distribution. Values may be negative, zero, or positive, provided
#'   that \eqn{1 + \alpha \mu > 0}.
#' @param log logical; if TRUE, probabilities are given as log(p).
#' @param log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE, probabilities p are \eqn{P[X\leq x]}
#'   otherwise, \eqn{P[X>x]}.
#'
#' @details
#' \code{dgp2} computes the density (PMF) of the Generalized Poisson
#' Version 2 Distribution.
#'
#' \code{pgp2} computes the CDF of the Generalized Poisson
#' Version 2 Distribution.
#'
#' \code{qgp2} computes the quantile function of the
#' Generalized Poisson Version 2 Distribution.
#'
#' \code{rgp2} generates random numbers from the Generalized Poisson
#' Version 2 Distribution.
#'
#' The probability mass function (PMF) for the Generalized Poisson
#' Version 2 distribution (GP-2) is:
#' \deqn{
#' f(y|\mu,\alpha)=
#' \frac{\mu(\mu+\alpha\mu y)^{y-1}
#' \exp\left(-\frac{\mu+\alpha\mu y}{1+\alpha\mu}\right)}
#' {(1+\alpha\mu)^y y!}
#' }
#'
#' where \eqn{\mu>0} is the mean and \eqn{y} is a non-negative integer.
#' This formulation uses the mean directly.
#'
#' The variance of the GP-2 distribution is:
#' \deqn{\sigma^2=\mu(1+\alpha\mu)^2}
#'
#' If \eqn{\alpha > 0}, the distribution is overdispersed. If
#' \eqn{\alpha = 0}, the distribution is equidispersed and reduces to the
#' ordinary Poisson distribution. If \eqn{\alpha < 0} and
#' \eqn{1 + \alpha\mu > 0}, the distribution is underdispersed and has
#' finite support.
#'
#' Under the underdispersed case, the largest possible value of \eqn{y}
#' is bounded above by \eqn{\left\lceil -1/\alpha \right\rceil - 1},
#' i.e., the largest integer satisfying \eqn{1 + \alpha y > 0}.
#'
#' @returns dgp2 gives the density, pgp2 gives the distribution
#'  function, qgp2 gives the quantile function, and rgp2 generates random
#'  deviates.
#'
#' The length of the result is determined by n for rgp2, and is the
#' maximum of the lengths of the numerical arguments for the other functions.
#' 
#' @references Consul, PoC, and Felix Famoye. "Generalized Poisson regression model."
#' Communications in Statistics-Theory and Methods 21.1 (1992): 89-109.
#'
#' Wang, Weizhen, and Felix Famoye. "Modeling household fertility decisions
#' with generalized Poisson regression." Journal of Population Economics 10.3
#' (1997): 273-283.
#'
#' Yang, Zhao, James W. Hardin, and Cheryl L. Addy. "A score test for
#' overdispersion in Poisson regression based on the generalized
#' Poisson-2 model." Journal of Statistical Planning and Inference 139.4
#' (2009): 1514-1521.
#'
#' Harris, Tammy, Zhao Yang, and James W. Hardin. "Modeling underdispersed
#' count data with generalized Poisson regression." Stata Journal 12.4
#' (2012): 736-747.
#'
#' @examples
#' dgp2(1, mu = 0.75, alpha = 0.3)
#' pgp2(c(0, 1, 2, 3, 5, 7, 9, 10), mu = 0.75, alpha = 0.25)
#' qgp2(c(0.1, 0.3, 0.5, 0.9, 0.95), mu = 0.75, alpha = -0.2)
#' rgp2(30, mu = 0.75, alpha = 0.2)
#'
#' @importFrom stats runif dpois ppois qpois rpois
#' @name GeneralizedPoisson2
#' @rdname GeneralizedPoisson2
#' @export
dgp2 <- function(x, mu = 1, alpha = 0, log = FALSE) {
  n <- max(length(x), length(mu), length(alpha), length(log))
  x <- rep_len(x, n)
  mu <- rep_len(mu, n)
  alpha <- rep_len(alpha, n)
  log <- rep_len(log, n)
  
  ok <- .gp2_validate(mu, alpha)
  out <- rep(NA_real_, n)
  
  valid_idx <- which(ok)
  if (length(valid_idx) > 0) {
    for (i in valid_idx) {
      xi <- x[i]
      
      if (is.na(xi) || is.na(log[i])) {
        out[i] <- NA_real_
        next
      }
      
      if (xi < 0 || abs(xi - round(xi)) > .Machine$double.eps^0.5) {
        out[i] <- if (isTRUE(log[i])) -Inf else 0
        next
      }
      
      lp <- .gp2_logpmf(xi, mu[i], alpha[i])
      out[i] <- if (isTRUE(log[i])) lp else exp(lp)
    }
  }
  
  out
}

#' @rdname GeneralizedPoisson2
#' @export
pgp2 <- function(q, mu = 1, alpha = 0, lower.tail = TRUE, log.p = FALSE) {
  n <- max(length(q), length(mu), length(alpha), length(lower.tail), length(log.p))
  q <- rep_len(q, n)
  mu <- rep_len(mu, n)
  alpha <- rep_len(alpha, n)
  lower.tail <- rep_len(lower.tail, n)
  log.p <- rep_len(log.p, n)
  
  ok <- .gp2_validate(mu, alpha)
  out <- rep(NA_real_, n)
  
  valid_idx <- which(ok)
  if (length(valid_idx) > 0) {
    for (i in valid_idx) {
      if (is.na(q[i]) || is.na(lower.tail[i]) || is.na(log.p[i])) {
        out[i] <- NA_real_
        next
      }
      
      out[i] <- .gp2_cdf_scalar(q[i], mu[i], alpha[i])
      
      if (!isTRUE(lower.tail[i])) out[i] <- 1 - out[i]
      if (isTRUE(log.p[i])) out[i] <- log(out[i])
    }
  }
  
  out
}

#' @rdname GeneralizedPoisson2
#' @export
qgp2 <- function(p, mu = 1, alpha = 0, lower.tail = TRUE, log.p = FALSE) {
  n <- max(length(p), length(mu), length(alpha), length(lower.tail), length(log.p))
  p <- rep_len(p, n)
  mu <- rep_len(mu, n)
  alpha <- rep_len(alpha, n)
  lower.tail <- rep_len(lower.tail, n)
  log.p <- rep_len(log.p, n)
  
  ok <- .gp2_validate(mu, alpha)
  out <- rep(NA_real_, n)
  
  valid_idx <- which(ok)
  if (length(valid_idx) > 0) {
    for (i in valid_idx) {
      pi <- p[i]
      
      if (is.na(pi) || is.na(lower.tail[i]) || is.na(log.p[i])) {
        out[i] <- NA_real_
        next
      }
      
      if (isTRUE(log.p[i])) pi <- exp(pi)
      if (!isTRUE(lower.tail[i])) pi <- 1 - pi
      
      if (!is.finite(pi) || pi < 0 || pi > 1) {
        out[i] <- NA_real_
        next
      }
      
      if (pi == 0) {
        out[i] <- 0
        next
      }
      
      if (pi == 1) {
        if (alpha[i] < 0) {
          out[i] <- .gp2_ymax(alpha[i])
        } else {
          out[i] <- Inf
        }
        next
      }
      
      if (alpha[i] == 0) {
        out[i] <- qpois(pi, lambda = mu[i])
        next
      }
      
      if (alpha[i] < 0) {
        ymax <- .gp2_ymax(alpha[i])
        y_seq <- 0:ymax
        
        logp_vec <- vapply(
          y_seq,
          function(yy) .gp2_logpmf_raw(yy, mu[i], alpha[i]),
          numeric(1)
        )
        logp_vec <- logp_vec - .gp2_log_norm_const(mu[i], alpha[i])
        
        p_y <- exp(logp_vec)
        cdf <- cumsum(p_y)
        
        idx <- which(cdf >= pi)[1]
        out[i] <- y_seq[idx]
        next
      }
      
      y <- 0L
      cdf <- .gp2_cdf_scalar(y, mu[i], alpha[i])
      
      max_iter <- 1e6L
      iter <- 0L
      
      while (cdf < pi) {
        y <- y + 1L
        cdf <- .gp2_cdf_scalar(y, mu[i], alpha[i])
        iter <- iter + 1L
        if (iter > max_iter) {
          warning("Quantile search exceeded the iteration limit.", call. = FALSE)
          out[i] <- NA_real_
          break
        }
      }
      
      if (is.na(out[i]) && iter <= max_iter) {
        out[i] <- y
      }
    }
  }
  
  out
}

#' @rdname GeneralizedPoisson2
#' @export
rgp2 <- function(n, mu = 1, alpha = 0) {
  if (length(n) != 1L || is.na(n) || !is.finite(n) || n < 0) {
    warning("`n` must be a single non-negative integer.", call. = FALSE)
    return(rep(NA_integer_, 1L))
  }
  
  n <- as.integer(n)
  
  if (length(mu) != 1L || length(alpha) != 1L) {
    warning("`mu` and `alpha` must be scalar values for `rgp2()`.", call. = FALSE)
    return(rep(NA_integer_, n))
  }
  
  ok <- .gp2_validate(mu, alpha)
  if (!ok) {
    return(rep(NA_integer_, n))
  }
  
  if (alpha == 0) {
    return(rpois(n, lambda = mu))
  }
  
  if (alpha < 0) {
    ymax <- .gp2_ymax(alpha)
    y_seq <- 0:ymax
    
    logp <- vapply(y_seq, function(yy) .gp2_logpmf_raw(yy, mu, alpha), numeric(1))
    logp <- logp - .gp2_log_norm_const(mu, alpha)
    prob <- exp(logp)
    
    return(sample(y_seq, size = n, replace = TRUE, prob = prob))
  }
  
  u <- runif(n)
  vapply(u, function(p) qgp2(p, mu = mu, alpha = alpha), numeric(1))
}
