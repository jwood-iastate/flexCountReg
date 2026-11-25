# Gamma renewal count distribution (Winkelmann, 1995)

.gcount_valid <- function(p) {
  is.finite(p$shape) & p$shape > 0 & is.finite(p$rate) & p$rate > 0 &
    is.finite(p$time) & p$time > 0
}

# log G(s; x) = log pgamma(s, x) (or its upper tail), recycled. `logs` is
# log(s), given when s itself may underflow (the mean-parameterized
# GAMMACOUNT family with small shapes); below s = exp(-700) the series
# s^x e^(-s) / Gamma(x + 1) sum_n s^n / ((x + 1)...(x + n)) is
# s^x / Gamma(x + 1) to double precision.
.gcount_lpg <- function(s, x, lower, logs = NULL) {
  out <- stats::pgamma(s, x, lower.tail = lower, log.p = TRUE)
  if (!is.null(logs)) {
    r <- .fcr_recycle(logs, x)
    tiny <- !is.na(r[[1]]) & r[[1]] <= -700
    if (any(tiny)) {
      out <- rep_len(out, length(r[[1]]))
      lg <- r[[2]][tiny] * r[[1]][tiny] - lgamma(r[[2]][tiny] + 1)
      out[tiny] <- if (lower) lg else .fcr_log1mexp(lg)
    }
  }
  out
}

# With s = rate * time and G(s; a) = pgamma(s, a) the gamma(a, 1) CDF:
#   P(N <= y) = 1 - G(s; (y + 1) shape),  P(N > y) = G(s; (y + 1) shape)
.gcount_logcdf <- function(y, p, lower) {
  .gcount_lpg(p$rate * p$time, (y + 1) * p$shape, !lower, p$logs)
}

# P(N = y) = G(s; y a) - G(s; (y + 1) a) = Q(s; (y + 1) a) - Q(s; y a), with
# Q = 1 - G. The difference is taken between the pair of tail probabilities
# that are not close to 1. `p$logs` (optional) is log(rate * time).
.gcount_logpmf <- function(y, p) {
  s <- p$rate * p$time
  a <- p$shape
  logs <- p$logs
  if (!is.null(logs)) {
    r <- .fcr_recycle(y, s, a, logs)
    y <- r[[1]]; s <- r[[2]]; a <- r[[3]]; logs <- r[[4]]
  }
  out <- .gcount_lpg(s, a, FALSE, logs)   # y = 0
  pos <- y > 0
  if (any(pos)) {
    s <- s[pos]; a <- a[pos]; yy <- y[pos]
    lgs <- if (is.null(logs)) NULL else logs[pos]
    la <- .gcount_lpg(s, yy * a, TRUE, lgs)          # log G(s; y a)
    lb <- .gcount_lpg(s, (yy + 1) * a, TRUE, lgs)    # log G(s; (y+1) a)
    ua <- .gcount_lpg(s, yy * a, FALSE, lgs)
    ub <- .gcount_lpg(s, (yy + 1) * a, FALSE, lgs)
    use_low <- la <= log(0.5)
    res <- ifelse(use_low, la + .fcr_log1mexp(lb - la),
                  ub + .fcr_log1mexp(ua - ub))
    res[(use_low & la == -Inf) | (!use_low & ub == -Inf)] <- -Inf
    out[pos] <- res
  }
  out
}

# Mean and variance:
#   E[N] = sum_{k >= 1} G(s; k a),  E[N^2] = sum_{k >= 1} (2k - 1) G(s; k a).
# The terms decrease in k; once k a > s they decrease faster than
# geometrically, and the sum stops when they are below 1e-17 of the total
# with a ratio below 1/2.
.gcount_moments <- function(shape, rate, time = 1) {
  r <- .fcr_recycle(shape, rate, time)
  p <- list(shape = r[[1]], rate = r[[2]], time = r[[3]])
  n <- length(r[[1]])
  m <- v <- rep(NA_real_, n)
  ok <- .gcount_valid(p)
  ok[is.na(ok)] <- FALSE
  bad <- !ok & !.fcr_any_na(p, n)
  m[bad] <- NaN
  v[bad] <- NaN
  for (i in which(ok)) {
    s <- p$rate[i] * p$time[i]
    a <- p$shape[i]
    s1 <- s2 <- 0
    k0 <- 1
    block <- max(64, ceiling(2 * s / a))
    repeat {
      k <- k0:(k0 + block - 1)
      g <- stats::pgamma(s, k * a)
      s1 <- s1 + sum(g)
      s2 <- s2 + sum((2 * k - 1) * g)
      last <- g[block]
      if (last < 1e-17 * max(s1, 1e-300) &&
          (last == 0 || last / g[block - 1] < 0.5)) {
        break
      }
      k0 <- k0 + block
    }
    m[i] <- s1
    v[i] <- s2 - s1^2
  }
  list(mean = m, var = v)
}

.gcount_qstart <- function(p, llow) {
  # Renewal-theory approximation: mean ~ s / a, variance ~ s / a^2
  s <- p$rate * p$time
  .fcr_qstart_normal(s / p$shape, s / p$shape^2 + 1, llow)
}

.gcount_rand <- function(n, p) {
  .fcr_q(stats::runif(n), p, .gcount_valid, .gcount_logcdf, .gcount_qstart,
         TRUE, FALSE)
}

#' Gamma Count Distribution
#'
#' Density, distribution function, quantile function, and random generation
#' for the gamma renewal count distribution (Winkelmann, 1995): the number of
#' events in a time window when the times between events are independent
#' gamma variables.
#'
#' @inheritParams DiscreteWeibull
#' @param shape vector of (positive) shape parameters \eqn{a} of the gamma
#'   inter-arrival times. \eqn{a = 1} is the Poisson distribution, \eqn{a > 1}
#'   gives underdispersion and \eqn{a < 1} overdispersion.
#' @param rate vector of (positive) rate parameters \eqn{b} of the gamma
#'   inter-arrival times.
#' @param time vector of (positive) lengths \eqn{t} of the observation
#'   window. Default 1.
#'
#' @details
#' Let the times between events be independent Gamma(shape = \eqn{a},
#' rate = \eqn{b}) variables and let \eqn{N(t)} be the number of events
#' completed by time \eqn{t} (the process starts with an event at time 0).
#' With \eqn{G(s; c)} the distribution function of a Gamma(\eqn{c}, 1)
#' variable at \eqn{s = bt} (\code{pgamma(s, c)}),
#' \deqn{P(N(t) \leq y) = 1 - G(bt; (y + 1)a), \qquad
#'   P(N(t) = y) = G(bt; ya) - G(bt; (y + 1)a),}
#' with \eqn{G(bt; 0) = 1}. Only the product \eqn{bt} enters the
#' distribution. With \eqn{a = 1} the inter-arrival times are exponential and
#' \eqn{N(t)} is Poisson with mean \eqn{bt}.
#'
#' The mean is \eqn{E[N(t)] = \sum_{k \geq 1} G(bt; ka)}, which is
#' approximately (but not exactly) \eqn{bt / a} for large \eqn{bt}; the
#' expected count is generally not proportional to the window length
#' \eqn{t}, so lengthening the window is not the same as a log-exposure
#' offset.
#'
#' Probabilities are differences of gamma distribution functions, taken
#' between the pair of tail probabilities that is not close to one. When the
#' shape \eqn{a} is very small the two terms are nearly equal and some
#' relative accuracy (about \eqn{\log_{10}(1/a)} digits) is lost.
#' Both tails of the distribution function are computed directly with
#' \code{pgamma()}. \code{rgcount()} uses inversion of the distribution
#' function.
#'
#' @returns \code{dgcount} gives the probability mass, \code{pgcount} the
#'   distribution function, \code{qgcount} the quantile function, and
#'   \code{rgcount} random deviates. The length of the result is determined
#'   by \code{n} for \code{rgcount}, and is the maximum of the lengths of
#'   the numerical arguments for the other functions.
#'
#' @references
#' Winkelmann, R. (1995). Duration dependence and dispersion in count-data
#' models. \emph{Journal of Business & Economic Statistics}, 13(4), 467-474.
#'
#' Kharrat, T., Boshnakov, G. N., McHale, I. and Baker, R. (2019). Flexible
#' regression models for count data based on renewal processes: the Countr
#' package. \emph{Journal of Statistical Software}, 90(13), 1-35.
#'
#' @examples
#' dgcount(0:6, shape = 2, rate = 6)
#' # shape = 1 gives the Poisson distribution with mean rate * time
#' all.equal(dgcount(0:10, shape = 1, rate = 1.5, time = 2), dpois(0:10, 3))
#' pgcount(3, shape = 0.5, rate = 1.5)
#' qgcount(c(0.1, 0.5, 0.9), shape = 2, rate = 6)
#' set.seed(6)
#' y <- rgcount(2000, shape = 3, rate = 12)
#' c(mean = mean(y), var = var(y))   # underdispersed
#'
#' @name GammaCount
#' @rdname GammaCount
#' @export
dgcount <- function(x, shape, rate, time = 1, log = FALSE) {
  .fcr_d(x, list(shape = shape, rate = rate, time = time), .gcount_valid,
         .gcount_logpmf, log)
}

#' @rdname GammaCount
#' @export
pgcount <- function(q, shape, rate, time = 1, lower.tail = TRUE,
                    log.p = FALSE) {
  .fcr_p(q, list(shape = shape, rate = rate, time = time), .gcount_valid,
         .gcount_logcdf, lower.tail, log.p)
}

#' @rdname GammaCount
#' @export
qgcount <- function(p, shape, rate, time = 1, lower.tail = TRUE,
                    log.p = FALSE) {
  .fcr_q(p, list(shape = shape, rate = rate, time = time), .gcount_valid,
         .gcount_logcdf, .gcount_qstart, lower.tail, log.p)
}

#' @rdname GammaCount
#' @export
rgcount <- function(n, shape, rate, time = 1) {
  .fcr_r(n, list(shape = shape, rate = rate, time = time), .gcount_valid,
         .gcount_rand)
}
