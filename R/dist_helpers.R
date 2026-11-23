# Internal helpers shared by the distribution functions.

# TRUE where x is a non-negative whole number.
.fcr_is_count <- function(x) {
  ok <- is.finite(x)
  out <- rep(FALSE, length(x))
  out[ok] <- x[ok] >= 0 & abs(x[ok] - round(x[ok])) <= sqrt(.Machine$double.eps)
  out
}

# Recycle a set of arguments to a common length (base-R style). Returns a list.
.fcr_recycle <- function(...) {
  args <- list(...)
  n <- max(vapply(args, length, integer(1)))
  if (any(vapply(args, length, integer(1)) == 0L)) n <- 0L
  lapply(args, rep_len, length.out = n)
}

# log K_nu(z), the log of the modified Bessel function of the second kind,
# vectorized over z and nu, and finite for large orders where besselK()
# overflows. K_nu = K_{-nu}, so the order is taken as |nu| = a + m with
# 0 <= a < 1. K_a and K_{a+1} come from besselK(..., expon.scaled = TRUE);
# higher orders use the upward recurrence
#   K_{v+1}(z) = K_{v-1}(z) + (2 v / z) K_v(z),
# which is numerically stable for K, carried out on the log scale.
.fcr_log_besselK <- function(z, nu) {
  r <- .fcr_recycle(z, nu)
  z <- r[[1]]
  v <- abs(r[[2]])
  out <- rep(NA_real_, length(z))
  ok <- is.finite(z) & z > 0 & is.finite(v)
  if (!any(ok)) return(out)
  z <- z[ok]
  v <- v[ok]
  m <- floor(v)
  a <- v - m
  lk0 <- log(besselK(z, a, expon.scaled = TRUE)) - z
  lk1 <- log(besselK(z, a + 1, expon.scaled = TRUE)) - z
  res <- ifelse(m == 0, lk0, lk1)
  if (any(m >= 2)) {
    for (k in 1:(max(m) - 1)) {
      act <- m > k
      # lk2 = log K_{a+k+1} from K_{a+k-1} (lk0) and K_{a+k} (lk1)
      lk2 <- lk1 + log(exp(lk0 - lk1) + 2 * (a + k) / z)
      lk0 <- lk1
      lk1 <- lk2
      done <- m == k + 1
      res[done] <- lk1[done]
      if (!any(act)) break
    }
  }
  out[ok] <- res
  out
}

# Quantile of a discrete distribution on 0, 1, 2, ... for one probability p.
# `cdf_at(y)` must return the CDF at the integer vector y. The support is
# searched in blocks that double in size, so the cost is O(answer), p = 1
# returns Inf and p outside [0, 1] returns NaN.
.fcr_qsearch <- function(p, cdf_at, max_y = 1e7) {
  if (is.na(p)) return(NA_real_)
  if (p < 0 || p > 1) return(NaN)
  if (p == 1) return(Inf)
  # Tolerance for rounding in the summed CDF near 1. It is relative for small
  # p: an absolute 1e-12 would return 0 for every p below 1e-12.
  tol <- min(1e-12, 1e-8 * p)
  K <- 32
  repeat {
    cdf <- cdf_at(0:K)
    if (all(is.nan(cdf))) return(NaN)      # invalid parameters
    hit <- which(cdf >= p - tol)
    if (length(hit)) return(hit[1] - 1)
    if (K >= max_y || all(is.na(cdf))) return(NA_real_)
    K <- min(2 * K, max_y)
  }
}

# Vectorized wrapper around .fcr_qsearch for parameters recycled with p.
# `cdf_fun(y, ...)` gets scalar parameters.
.fcr_qcount <- function(p, cdf_fun, ...) {
  r <- .fcr_recycle(p, ...)
  pp <- r[[1]]
  pars <- r[-1]
  vapply(seq_along(pp), function(i) {
    args_i <- lapply(pars, `[[`, i)
    .fcr_qsearch(pp[i], function(y) do.call(cdf_fun, c(list(y), args_i)))
  }, numeric(1))
}

# Number of draws, following the base-R convention that length(n) > 1 means
# length(n) draws.
.fcr_rn <- function(n) {
  if (length(n) > 1L) return(length(n))
  if (length(n) != 1L || is.na(n) || !is.finite(n) || n < 0) {
    stop("`n` must be a non-negative number.", call. = FALSE)
  }
  as.integer(floor(n))
}

# Apply a scalar-parameter sampler `scalar_r(k, ...)` once per distinct
# parameter combination. `pars` is a named list of vectors of length n;
# `ok` flags valid rows (others are NA).
.fcr_r_grouped <- function(n, pars, ok, scalar_r) {
  out <- rep(NA_real_, n)
  if (!any(ok)) return(out)
  key <- do.call(paste, c(lapply(pars, function(v) format(v[ok], digits = 17)),
                          sep = "\r"))
  idx_ok <- which(ok)
  for (k in unique(key)) {
    idx <- idx_ok[key == k]
    args <- lapply(pars, `[`, idx[1])
    out[idx] <- do.call(scalar_r, c(list(length(idx)), args))
  }
  out
}

# Draws from the Lindley(theta) distribution as a gamma mixture:
# Gamma(1, theta) with probability theta / (theta + 1), else Gamma(2, theta).
.fcr_rlindley_mix <- function(n, theta) {
  shape <- 1 + (stats::runif(n) > theta / (theta + 1))
  stats::rgamma(n, shape = shape, rate = theta)
}

# The base R convention for invalid distribution parameters: NaN for that
# element and one warning.
.fcr_warn_nan <- function(bad) {
  if (any(bad, na.rm = TRUE)) warning("NaNs produced", call. = FALSE)
}

.fcr_warn_na <- function(bad, what) {
  if (any(bad)) {
    warning("NAs produced: ", what, call. = FALSE)
  }
}

# Shared contract of the count-distribution p-functions: arguments are
# recycled; NA q (or NaN) or a missing parameter gives NA; for valid
# parameters q < 0 (including -Inf) gives P(Y <= q) = 0 and q = Inf gives 1,
# while invalid parameters give the core's NA/NaN at every q; the upper tail
# and log scale follow. Only
# finite q >= 0 reach `core`, which gets the matching elements of the
# vector parameters in `pars` (NULL entries are dropped: NULL is each
# core's default) and the remaining arguments in `...`. With `logdens`
# (function(x, <pars>) returning log P(Y = x)), small upper tails are
# summed directly instead of taken as 1 - CDF.
.fcr_p_contract <- function(core, q, pars, lower.tail, log.p, ...,
                            logdens = NULL) {
  for (nm in c("lower.tail", "log.p")) {
    v <- get(nm)
    if (!is.logical(v) || length(v) != 1L || is.na(v)) {
      stop("'", nm, "' must be TRUE or FALSE", call. = FALSE)
    }
  }
  lens <- c(length(q), vapply(Filter(Negate(is.null), pars), length, 1L))
  if (any(lens == 0L)) return(numeric(0))
  n <- max(lens)
  q <- rep_len(as.numeric(q), n)
  sub <- function(p, ok) if (length(p) == 1L) p else rep_len(p, n)[ok]
  out <- rep(NA_real_, n)
  pars_nn <- Filter(Negate(is.null), pars)
  # Rows with a missing parameter are NA whatever q is
  par_na <- Reduce(`|`, lapply(pars_nn, function(p) is.na(rep_len(p, n))),
                   rep(FALSE, n))
  low <- !is.na(q) & q < 0 & !par_na
  high <- !is.na(q) & q == Inf & !par_na
  # The endpoints are 0 or 1 only for valid parameters: probe the core at
  # q = 0, which gives NA or NaN for parameters outside its domain
  ends <- which(low | high)
  if (length(ends)) {
    probe <- suppressWarnings(do.call(core, c(
      list(rep(0, length(ends))),
      lapply(pars_nn, function(p) if (length(p) == 1L) p else
        rep_len(p, n)[ends]),
      list(lower.tail = TRUE, log.p = FALSE), list(...))))
    bad <- ends[is.na(probe)]
    out[bad] <- probe[is.na(probe)]
    low[bad] <- FALSE
    high[bad] <- FALSE
  }
  out[low] <- if (lower.tail) 0 else 1
  out[high] <- if (lower.tail) 1 else 0
  if (log.p) out[low | high] <- log(out[low | high])
  ok <- !is.na(q) & is.finite(q) & q >= 0 & !par_na
  if (any(ok)) {
    pars <- Filter(Negate(is.null), pars)   # NULL is each core's default
    out[ok] <- do.call(core, c(list(q[ok]), lapply(pars, sub, ok = ok),
                               list(lower.tail = lower.tail, log.p = log.p),
                               list(...)))
    # Small upper tails: 1 - CDF has no relative precision below ~1e-16, so
    # sum the log-PMF beyond q directly
    if (!lower.tail && !is.null(logdens)) {
      small <- which(ok & !is.na(out) &
                       (if (log.p) out < log(1e-8) else out < 1e-8))
      for (i in small) {
        pi <- lapply(pars, function(p) if (length(p) == 1L) p else
          rep_len(p, n)[i])
        lt <- .fcr_log_upper_tail(floor(q[i]), function(x)
          do.call(logdens, c(list(x), pi)))
        if (!is.na(lt)) out[i] <- if (log.p) lt else exp(lt)
      }
    }
  }
  out
}

# log P(Y > q) by summing log-probabilities logd(x) for x = q + 1, q + 2, ...
# in blocks until the terms are negligible (NA if they do not decrease).
.fcr_log_upper_tail <- function(q, logd, block = 200L, max_terms = 1e6) {
  total <- -Inf
  start <- q + 1
  while (start - q <= max_terms) {
    x <- start + 0:(block - 1L)
    lp <- logd(x)
    if (anyNA(lp)) return(NA_real_)
    m <- max(lp, total)
    if (!is.finite(m)) return(total)
    total <- m + log(exp(total - m) + sum(exp(lp - m)))
    # Stop when the block's last term is negligible and the terms decrease
    if (lp[block] < total + log(1e-17) && lp[block] <= lp[1]) return(total)
    start <- start + block
  }
  NA_real_
}
