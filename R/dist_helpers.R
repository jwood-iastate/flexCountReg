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
  K <- 32
  repeat {
    cdf <- cdf_at(0:K)
    hit <- which(cdf >= p - 1e-12)
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

.fcr_warn_na <- function(bad, what) {
  if (any(bad)) {
    warning("NAs produced: ", what, call. = FALSE)
  }
}
