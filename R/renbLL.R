renb_ll <- function(y, mu, a, b, panels) { 
  # Random Effects Negative Binomial with Beta Distributed Random Effects (NB1)
  
  
  if (any(!is.finite(mu)) || any(mu <= 0) ||
      !is.finite(a) || a <= 0 ||
      !is.finite(b) || b <= 0) {
    return(rep(-Inf, length(unique(panels))))
  }
  
  # Panel-level sums (one row per panel, in sorted panel order)
  sums <- rowsum(cbind(mu, y, .fcr_lgamma_diff(mu, y) - lgamma(y + 1)),
                 group = panels, reorder = TRUE)
  df_model <- list(sum_mu = sums[, 1], sum_y = sums[, 2],
                   log_prod_element = sums[, 3])
  
  # log P = lbeta(a + M, b + Y) - lbeta(a, b) + sum(log terms), with
  # M = sum(mu) and Y = sum(y), computed in two parts so that no two huge
  # terms cancel (single lgamma terms can be ~1e300 when a, b or M is large):
  #  (1) lbeta(a + M, c) - lbeta(a, c), c = b + Y, which equals both
  #      D(a, M) - D(a + c, M) and D(a, c) - D(a + M, c), with
  #      D(x, d) = lgamma(x + d) - lgamma(x); the form with the smaller
  #      increment min(M, c) is used;
  #  (2) lbeta(a, b + Y) - lbeta(a, b) = -sum_{k < Y} log1p(a / (b + k)),
  #      exact for whole Y and accurate for every a and b.
  M <- df_model$sum_mu
  Y <- df_model$sum_y
  # A panel mean that overflows has no finite likelihood
  if (any(!is.finite(M))) return(rep(-Inf, length(M)))
  cc <- b + Y
  small_M <- M <= cc
  t1 <- numeric(length(M))
  t1[small_M] <- .fcr_lgamma_diff(a, M[small_M]) -
    .fcr_lgamma_diff(a + cc[small_M], M[small_M])
  t1[!small_M] <- .fcr_lgamma_diff(a, cc[!small_M]) -
    .fcr_lgamma_diff(a + M[!small_M], cc[!small_M])
  Yi <- round(Y)
  k <- sequence(Yi) - 1
  pid <- rep(seq_along(Yi), Yi)
  t2 <- numeric(length(Yi))
  if (length(k)) {
    s2 <- rowsum(log1p(a / (b + k)), pid, reorder = TRUE)
    t2[as.integer(rownames(s2))] <- -s2[, 1]
  }
  log_P <- t1 + t2 + df_model$log_prod_element

  # Sum across panels for total log-likelihood
  total_loglik <- log_P
  
  return(total_loglik)
}

# lgamma(x + d) - lgamma(x) for x > 0, d >= 0, accurate when x is much larger
# than d (where the direct difference cancels to zero in floating point).
.fcr_lgamma_diff <- function(x, d) {
  if (!length(x) || !length(d)) return(numeric(0))
  n <- max(length(x), length(d))
  x <- rep_len(x, n); d <- rep_len(d, n)
  out <- lgamma(x + d) - lgamma(x)
  big <- is.finite(x) & x > 1e6 * (1 + d^2)
  if (any(big)) {
    xb <- x[big]; db <- d[big]
    # Stirling series: d log x + d (d - 1) / (2 x) + O(d^3 / x^2)
    out[big] <- db * log(xb) + db * (db - 1) / (2 * xb)
  }
  out
}

# lgamma(x + d) - lgamma(x) - d * log(x): tends to 0 as x -> Inf (the
# Poisson limit of negative binomial panel terms with x = 1 / alpha), where
# the direct expression loses all significant digits. For whole d (counts)
# it is the finite sum of log1p(k / x), k = 0, ..., d - 1, which is accurate
# for every x; otherwise a Stirling series is used for large x.
.fcr_lgamma_diff_rel <- function(x, d) {
  if (!length(x) || !length(d)) return(numeric(0))
  n <- max(length(x), length(d))
  x <- rep_len(x, n); d <- rep_len(d, n)
  out <- lgamma(x + d) - lgamma(x) - d * log(x)
  int <- is.finite(x) & x > 0 & is.finite(d) & d >= 0 & d == round(d) &
    d <= 1e5
  if (any(int)) {
    out[int] <- vapply(which(int), function(i) {
      if (d[i] == 0) 0 else sum(log1p((0:(d[i] - 1)) / x[i]))
    }, numeric(1))
  }
  big <- !int & is.finite(x) & x > 1e4 * (1 + d^2)
  if (any(big)) {
    xb <- x[big]; db <- d[big]
    out[big] <- db * (db - 1) / (2 * xb) -
      db * (db - 1) * (2 * db - 1) / (12 * xb^2) +
      db^2 * (db - 1)^2 / (12 * xb^3)
  }
  out
}

# digamma(x + d) - digamma(x): the finite sum of 1 / (x + k) for whole d,
# otherwise a series for large x
.fcr_digamma_diff <- function(x, d) {
  n <- max(length(x), length(d))
  x <- rep_len(x, n); d <- rep_len(d, n)
  out <- digamma(x + d) - digamma(x)
  int <- is.finite(x) & x > 0 & is.finite(d) & d >= 0 & d == round(d) &
    d <= 1e5
  if (any(int)) {
    out[int] <- vapply(which(int), function(i) {
      if (d[i] == 0) 0 else sum(1 / (x[i] + 0:(d[i] - 1)))
    }, numeric(1))
  }
  big <- !int & is.finite(x) & x > 1e4 * (1 + d^2)
  if (any(big)) {
    xb <- x[big]; db <- d[big]
    out[big] <- db / xb - db * (db - 1) / (2 * xb^2) +
      db * (db - 1) * (2 * db - 1) / (6 * xb^3)
  }
  out
}
