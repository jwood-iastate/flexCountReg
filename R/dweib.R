# Discrete Weibull (type I, Nakagawa and Osaki, 1975)

.dweib_valid <- function(p) {
  is.finite(p$theta) & p$theta > 0 & p$theta < 1 &
    is.finite(p$shape) & p$shape > 0
}

# c = -log(theta) > 0; P(Y >= y) = exp(-c y^shape)
.dweib_c <- function(theta) -log(theta)

.dweib_logpmf <- function(y, p) .dweib_logpmf_c(y, .dweib_c(p$theta), p$shape)

.dweib_logcdf <- function(y, p, lower) {
  .dweib_logcdf_c(y, .dweib_c(p$theta), p$shape, lower)
}

# The same in terms of c = -log(theta) (used by the mean-parameterized
# regression family, where theta can be too close to 1 to store)
.dweib_logpmf_c <- function(y, cc, b) {
  # (y + 1)^b - y^b without cancellation: y^b expm1(b log1p(1 / y))
  diff <- ifelse(y == 0, 1, y^b * expm1(b * log1p(1 / pmax(y, 1))))
  -cc * y^b + .fcr_log1mexp(-cc * diff)
}

.dweib_logcdf_c <- function(y, cc, b, lower) {
  a <- cc * (y + 1)^b
  if (lower) .fcr_log1mexp(-a) else -a
}

# Exact inversion of the survival function, checked by the bisection
.dweib_qstart <- function(p, llow) {
  t <- -.fcr_log1mexp(llow)          # -log P(Y > y) target
  y <- ceiling((t / .dweib_c(p$theta))^(1 / p$shape) - 1)
  y[!is.finite(y)] <- 0
  pmax(y, 0)
}

# Y = floor(X), X continuous Weibull with survival exp(-c x^shape)
.dweib_rand <- function(n, p) {
  floor((stats::rexp(n) / .dweib_c(p$theta))^(1 / p$shape))
}

# Mean and variance. E[Y] = sum_{k >= 1} f(k) and
# E[Y^2] = sum_{k >= 1} (2k - 1) f(k) with f(k) = exp(-c k^b). The sums are
# taken directly up to K, where the terms are below exp(-50) or K = 1e5; the
# remainder is added with the Euler-Maclaurin formula
#   sum_{k > K} g(k) = int_K^Inf g - g(K) / 2 - g'(K) / 12 + ...,
# using closed-form incomplete gamma integrals. Heavy tails (small shape)
# therefore do not require summing an enormous number of terms.
.dweib_moments <- function(theta, shape) {
  r <- .fcr_recycle(theta, shape)
  theta <- r[[1]]; shape <- r[[2]]
  m <- v <- rep(NA_real_, length(theta))
  ok <- .dweib_valid(list(theta = theta, shape = shape))
  m[!ok & !is.na(theta) & !is.na(shape)] <- NaN
  v[!ok & !is.na(theta) & !is.na(shape)] <- NaN
  for (i in which(ok)) {
    cc <- -log(theta[i]); b <- shape[i]
    K <- min(ceiling((50 / cc)^(1 / b)), 1e5)
    k <- seq_len(K)
    f <- exp(-cc * k^b)
    s1 <- sum(f)
    s2 <- sum((2 * k - 1) * f)
    if (K < ceiling((50 / cc)^(1 / b))) {
      u <- cc * K^b
      fK <- exp(-u)
      dfK <- -cc * b * K^(b - 1) * fK
      # int_K^Inf x^(j - 1) exp(-c x^b) dx = c^(-j/b) Gamma(j/b) Q(j/b, cK^b) / b
      I <- function(j) exp(lgamma(j / b) - (j / b) * log(cc) - log(b) +
                             stats::pgamma(u, j / b, lower.tail = FALSE,
                                           log.p = TRUE))
      s1 <- s1 + I(1) - fK / 2 - dfK / 12
      gK <- (2 * K - 1) * fK
      dgK <- 2 * fK + (2 * K - 1) * dfK
      s2 <- s2 + 2 * I(2) - I(1) - gK / 2 - dgK / 12
    }
    m[i] <- s1
    v[i] <- s2 - s1^2
  }
  list(mean = m, var = v)
}

#' Discrete Weibull Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the (type I) discrete Weibull distribution of Nakagawa and Osaki
#' (1975).
#'
#' @param x vector of (non-negative integer) quantiles.
#' @param q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is
#'   taken to be the number required.
#' @param theta vector of values of the parameter \eqn{\theta}, with
#'   \eqn{0 < \theta < 1}. It equals \eqn{P(Y \geq 1)}.
#' @param shape vector of values of the shape parameter \eqn{\beta > 0}.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#'   \eqn{P(Y \leq y)}, otherwise \eqn{P(Y > y)}.
#'
#' @details
#' The discrete Weibull distribution has survival function
#' \deqn{P(Y \geq y) = \theta^{y^\beta}, \quad y = 0, 1, 2, \ldots}
#' so that
#' \deqn{P(Y = y) = \theta^{y^\beta} - \theta^{(y+1)^\beta}, \qquad
#'   P(Y \leq y) = 1 - \theta^{(y+1)^\beta}.}
#' It is the distribution of \eqn{\lfloor X \rfloor} where \eqn{X} has a
#' continuous Weibull distribution with shape \eqn{\beta} and survival
#' function \eqn{\exp(-c x^\beta)}, \eqn{c = -\log\theta}; \code{rdweib()}
#' uses that representation. With \eqn{\beta = 1} it is the geometric
#' distribution with success probability \eqn{1 - \theta}
#' (\code{dgeom(x, 1 - theta)}).
#'
#' The moments have no closed form:
#' \deqn{E[Y] = \sum_{k \geq 1} \theta^{k^\beta}, \qquad
#'   E[Y^2] = \sum_{k \geq 1} (2k - 1)\theta^{k^\beta}.}
#' The distribution can be under- or overdispersed relative to the Poisson,
#' which makes it a flexible model for counts (Klakattawi et al., 2018).
#'
#' This is a discrete response distribution; it is unrelated to the
#' Poisson-Weibull mixture of \code{\link{dpoisweibull}}.
#'
#' All computations use the logarithms of the tail probabilities, so the
#' upper tail \code{lower.tail = FALSE} keeps its relative accuracy far into
#' the tail. \code{qdweib()} inverts the survival function exactly and
#' returns the smallest \eqn{y} with \eqn{P(Y \leq y) \geq p}.
#'
#' Invalid parameter values give \code{NaN} with a warning, and non-integer
#' \code{x} has probability 0 (with a warning), as in base R.
#'
#' @returns \code{ddweib} gives the probability mass, \code{pdweib} the
#'   distribution function, \code{qdweib} the quantile function, and
#'   \code{rdweib} random deviates. The length of the result is determined
#'   by \code{n} for \code{rdweib}, and is the maximum of the lengths of the
#'   numerical arguments for the other functions (zero if any is zero).
#'
#' @references
#' Nakagawa, T. and Osaki, S. (1975). The discrete Weibull distribution.
#' \emph{IEEE Transactions on Reliability}, R-24(5), 300-301.
#'
#' Klakattawi, H., Vinciotti, V. and Yu, K. (2018). A simple and adaptive
#' dispersion regression model for count data. \emph{Entropy}, 20(2), 142.
#'
#' @examples
#' ddweib(0:5, theta = 0.7, shape = 1.5)
#' # shape = 1 gives the geometric distribution
#' all.equal(ddweib(0:10, theta = 0.6, shape = 1), dgeom(0:10, prob = 0.4))
#' pdweib(3, theta = 0.7, shape = 1.5)
#' pdweib(60, theta = 0.7, shape = 1.5, lower.tail = FALSE, log.p = TRUE)
#' qdweib(c(0.1, 0.5, 0.9), theta = 0.7, shape = 1.5)
#' set.seed(1)
#' y <- rdweib(1000, theta = 0.9, shape = 0.8)
#' c(mean = mean(y), var = var(y))   # overdispersed
#'
#' @name DiscreteWeibull
#' @rdname DiscreteWeibull
#' @export
ddweib <- function(x, theta, shape, log = FALSE) {
  .fcr_d(x, list(theta = theta, shape = shape), .dweib_valid,
         .dweib_logpmf, log)
}

#' @rdname DiscreteWeibull
#' @export
pdweib <- function(q, theta, shape, lower.tail = TRUE, log.p = FALSE) {
  .fcr_p(q, list(theta = theta, shape = shape), .dweib_valid, .dweib_logcdf,
         lower.tail, log.p)
}

#' @rdname DiscreteWeibull
#' @export
qdweib <- function(p, theta, shape, lower.tail = TRUE, log.p = FALSE) {
  .fcr_q(p, list(theta = theta, shape = shape), .dweib_valid, .dweib_logcdf,
         .dweib_qstart, lower.tail, log.p)
}

#' @rdname DiscreteWeibull
#' @export
rdweib <- function(n, theta, shape) {
  .fcr_r(n, list(theta = theta, shape = shape), .dweib_valid, .dweib_rand)
}
