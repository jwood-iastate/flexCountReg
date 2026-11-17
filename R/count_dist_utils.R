# Shared machinery for the count distributions added in flexCountReg 0.1.2
# (discrete Weibull, double Poisson, Delaporte, generalized Hermite,
# hyper-Poisson, gamma renewal count, and Polya-Aeppli).
#
# Each family supplies
#   valid(pars)                 TRUE where the (non-missing) parameters are in
#                               the parameter space;
#   logpmf(y, pars)             log P(Y = y) for whole numbers y >= 0;
#   logcdf(y, pars, lower)      log P(Y <= y) (lower = TRUE) or log P(Y > y),
#                               each computed directly (never as 1 - other);
#   start(pars, llow)           a rough guess of the quantile at log lower
#                               probability llow (only used to bracket the
#                               search; any value is safe);
#   rand(n, pars)               draws (optional; defaults to inversion).
# `pars` is a named list of parameter vectors of equal length. The public
# d/p/q/r functions are thin wrappers around .fcr_d(), .fcr_p(), .fcr_q() and
# .fcr_r() below, which implement the argument contract shared with base R:
# recycling, zero-length inputs, NA/NaN, invalid parameters (NaN with one
# warning), non-integer x (probability 0 with a warning), lower.tail and
# log.p.

# .fcr_log1mexp(a) = log(1 - exp(a)) for a <= 0 is defined in helpers.R.

# Log of the sum of exp(x)
.fcr_logsumexp <- function(x) {
  m <- max(x)
  if (!is.finite(m)) return(m)
  m + log(sum(exp(x - m)))
}

.fcr_check_flag <- function(v, nm) {
  if (!is.logical(v) || length(v) != 1L || is.na(v)) {
    stop("'", nm, "' must be TRUE or FALSE", call. = FALSE)
  }
}

# Recycle x and the parameters; returns list(x = , pars = , n = ).
.fcr_recycle_pars <- function(x, pars) {
  r <- do.call(.fcr_recycle, c(list(x), unname(pars)))
  list(x = as.numeric(r[[1]]),
       pars = stats::setNames(lapply(r[-1], as.numeric), names(pars)),
       n = length(r[[1]]))
}

.fcr_any_na <- function(pars, n) {
  Reduce(`|`, lapply(pars, is.na), rep(FALSE, n))
}

.fcr_subset_pars <- function(pars, idx) lapply(pars, `[`, idx)

# Rows with non-missing, valid parameters
.fcr_valid_rows <- function(valid, pars, na) {
  ok <- rep(FALSE, length(na))
  if (any(!na)) {
    v <- valid(.fcr_subset_pars(pars, !na))
    ok[!na] <- !is.na(v) & v
  }
  ok
}

# ---- density ----------------------------------------------------------------

.fcr_d <- function(x, pars, valid, logpmf, log) {
  .fcr_check_flag(log, "log")
  r <- .fcr_recycle_pars(x, pars)
  n <- r$n
  if (n == 0L) return(numeric(0))
  x <- r$x
  pars <- r$pars
  na <- is.na(x) | .fcr_any_na(pars, n)
  ok <- .fcr_valid_rows(valid, pars, na)
  bad <- !na & !ok
  lp <- rep(-Inf, n)
  nonint <- ok & is.finite(x) & x >= 0 & !.fcr_is_count(x)
  if (any(nonint)) {
    warning(sprintf("non-integer x = %f", x[which(nonint)[1]]),
            call. = FALSE)
  }
  cnt <- ok & .fcr_is_count(x)
  if (any(cnt)) {
    lp[cnt] <- logpmf(round(x[cnt]), .fcr_subset_pars(pars, cnt))
  }
  lp[bad] <- NaN
  lp[na] <- NA
  .fcr_warn_nan(bad)
  if (log) lp else exp(lp)
}

# ---- distribution function --------------------------------------------------

# Core used through .fcr_p_contract(), which handles recycling, NA, q < 0,
# q = Inf and the argument checks; the core gets finite q >= 0 and
# non-missing parameters.
.fcr_make_pcore <- function(valid, logcdf) {
  function(q, ..., lower.tail = TRUE, log.p = FALSE) {
    r <- .fcr_recycle_pars(q, list(...))
    out <- rep(NaN, r$n)
    if (r$n == 0L) return(numeric(0))
    ok <- valid(r$pars)
    ok <- !is.na(ok) & ok
    if (any(ok)) {
      out[ok] <- logcdf(floor(r$x[ok] + 1e-7), .fcr_subset_pars(r$pars, ok),
                        lower.tail)
    }
    .fcr_warn_nan(!ok)
    if (log.p) out else exp(out)
  }
}

.fcr_p <- function(q, pars, valid, logcdf, lower.tail, log.p) {
  .fcr_p_contract(.fcr_make_pcore(valid, logcdf), q, pars, lower.tail, log.p)
}

# ---- quantile function ------------------------------------------------------

# Smallest whole number y >= 0 with P(Y <= y) >= p. Vectorized bisection: the
# condition is evaluated for all unresolved elements at once. To keep tiny
# tail probabilities meaningful the comparison uses the lower tail when
# P(Y <= y) <= 1/2 is targeted and the upper tail otherwise (on the log
# scale, with a relative fuzz of 1e-12 as in base R's discrete quantiles).
.fcr_q <- function(p, pars, valid, logcdf, start, lower.tail, log.p,
                   max_y = 2^52) {
  .fcr_check_flag(lower.tail, "lower.tail")
  .fcr_check_flag(log.p, "log.p")
  r <- .fcr_recycle_pars(p, pars)
  n <- r$n
  if (n == 0L) return(numeric(0))
  p <- r$x
  pars <- r$pars
  na <- is.na(p) | .fcr_any_na(pars, n)
  ok <- .fcr_valid_rows(valid, pars, na)
  pbad <- !na & (if (log.p) p > 0 else (p < 0 | p > 1))
  out <- rep(NA_real_, n)
  out[!na & (!ok | pbad)] <- NaN
  .fcr_warn_nan(!na & (!ok | pbad))
  act <- which(ok & !pbad)
  if (!length(act)) return(out)
  lp <- if (log.p) p[act] else log(p[act])
  # log P(Y <= y) and log P(Y > y) targets
  if (lower.tail) {
    llow <- lp
    lupp <- .fcr_log1mexp(lp)
  } else {
    lupp <- lp
    llow <- .fcr_log1mexp(lp)
  }
  res <- rep(NA_real_, length(act))
  res[llow == -Inf] <- 0
  res[lupp == -Inf] <- Inf
  todo <- which(is.na(res))
  if (length(todo)) {
    sp <- .fcr_subset_pars(pars, act[todo])
    res[todo] <- .fcr_qbisect(llow[todo], lupp[todo], sp, logcdf,
                              start(sp, llow[todo]), max_y)
  }
  out[act] <- res
  out
}

.fcr_qbisect <- function(llow, lupp, pars, logcdf, s, max_y) {
  use_low <- llow <= log(0.5)
  cond <- function(y, idx) {
    res <- logical(length(idx))
    lo <- use_low[idx]
    pp <- .fcr_subset_pars(pars, idx)
    if (any(lo)) {
      res[lo] <- logcdf(y[lo], .fcr_subset_pars(pp, lo), TRUE) >=
        llow[idx][lo] - 1e-12
    }
    if (any(!lo)) {
      res[!lo] <- logcdf(y[!lo], .fcr_subset_pars(pp, !lo), FALSE) <=
        lupp[idx][!lo] + 1e-12
    }
    res[is.na(res)] <- FALSE
    res
  }
  m <- length(llow)
  s <- floor(s)
  s[!is.finite(s) | s < 0] <- 0
  s <- pmin(s, max_y)
  ans <- rep(NA_real_, m)
  hi <- s
  lo <- rep(-1, m)
  c_hi <- cond(hi, seq_len(m))
  # Quick exit: the guess is the answer
  chk <- which(c_hi & hi >= 1)
  if (length(chk)) {
    c_below <- cond(hi[chk] - 1, chk)
    ans[chk[!c_below]] <- hi[chk[!c_below]]
    lo[chk[c_below]] <- -1
  }
  ans[c_hi & hi == 0] <- 0
  # Expand the bracket upward where the guess is too small
  up <- which(!c_hi)
  while (length(up)) {
    lo[up] <- hi[up]
    hi[up] <- pmin(2 * hi[up] + 1, max_y)
    c_up <- cond(hi[up], up)
    stuck <- !c_up & hi[up] >= max_y
    if (any(stuck)) {
      ans[up[stuck]] <- NaN
      warning("quantile search reached the largest supported count",
              call. = FALSE)
    }
    up <- up[!c_up & !stuck]
  }
  # Bisection on (lo, hi]: cond(lo) is FALSE (or lo = -1), cond(hi) TRUE
  act <- which(is.na(ans))
  while (length(act)) {
    done <- hi[act] - lo[act] <= 1
    ans[act[done]] <- hi[act[done]]
    act <- act[!done]
    if (!length(act)) break
    mid <- floor((lo[act] + hi[act]) / 2)
    cm <- cond(mid, act)
    hi[act[cm]] <- mid[cm]
    lo[act[!cm]] <- mid[!cm]
  }
  ans
}

# Normal-approximation starting guess from a mean and variance
.fcr_qstart_normal <- function(mean, var, llow) {
  z <- stats::qnorm(pmax(llow, -700), log.p = TRUE)
  floor(pmax(mean + z * sqrt(pmax(var, 0)), 0))
}

# ---- random generation ------------------------------------------------------

.fcr_r <- function(n, pars, valid, rand) {
  n <- .fcr_rn(n)
  pars <- lapply(pars, function(v) rep_len(as.numeric(v), n))
  out <- rep(NA_real_, n)
  if (n == 0L) return(numeric(0))
  na <- .fcr_any_na(pars, n)
  ok <- .fcr_valid_rows(valid, pars, na)
  .fcr_warn_na(!ok, "invalid or missing parameter values.")
  if (any(ok)) out[ok] <- rand(sum(ok), .fcr_subset_pars(pars, ok))
  out
}

# ---- normalized weight tables -----------------------------------------------
#
# For distributions given by unnormalized log weights logw(y, <scalar pars>)
# on y = 0, 1, 2, ... (double Poisson, hyper-Poisson), the weights are
# computed on 0..K with K large enough that
#   * the last weight is at least 50 (log units) below the largest weight,
#     and
#   * the log weights are decreasing with non-increasing successive ratios at
#     K (both families have eventually decreasing ratios w(y+1)/w(y)), so the
#     omitted tail is at most w(K) r / (1 - r) with r = w(K) / w(K - 1).
# That geometric bound is added to the normalizing constant. K doubles until
# the conditions hold or the work budget `max_K` is exhausted (the result is
# then NA with a warning). Upper tails beyond the table are summed locally
# from the log weights, so arbitrarily large counts need no larger table.

.fcr_logaddexp <- function(a, b) {
  m <- pmax(a, b)
  out <- m + log1p(exp(-abs(a - b)))
  out[m == -Inf] <- -Inf
  out
}

.fcr_wtable <- function(logw, pars1, k0, max_K = 1e7) {
  K <- max(ceiling(k0), 16)
  repeat {
    if (K > max_K) return(NULL)
    lw <- do.call(logw, c(list(0:K), pars1))
    L <- length(lw)
    top <- max(lw)
    if (!is.finite(top)) return(NULL)
    d1 <- lw[L] - lw[L - 1]
    d2 <- lw[L - 1] - lw[L - 2]
    if (lw[L] < top - 50 && d1 < 0 && d1 <= d2) break
    K <- 2 * K
  }
  r <- exp(d1)
  ltail <- lw[L] + log(r) - log1p(-r)
  cum <- cumlogsumexp_cpp(lw)
  logZ <- .fcr_logaddexp(cum[L], ltail)
  # log sum_{j >= k} w_j for k = 0..K, including the omitted tail
  rcum <- .fcr_logaddexp(rev(cumlogsumexp_cpp(rev(lw))), ltail)
  list(lw = lw, logZ = logZ, cum = cum, rcum = rcum, K = K)
}

# Apply fun(idx, table, pars1) once per distinct parameter combination.
.fcr_by_table <- function(pars, logw, k0fun, fun) {
  n <- length(pars[[1]])
  out <- rep(NA_real_, n)
  key <- do.call(paste, c(lapply(pars, function(v) format(v, digits = 17)),
                          sep = "
"))
  failed <- FALSE
  for (k in unique(key)) {
    idx <- which(key == k)
    pars1 <- lapply(pars, `[`, idx[1])
    tab <- .fcr_wtable(logw, pars1, do.call(k0fun, pars1))
    if (is.null(tab)) {
      failed <- TRUE
      next
    }
    out[idx] <- fun(idx, tab, pars1)
  }
  if (failed) {
    warning("the normalizing constant could not be computed within the ",
            "work budget for some parameter values; NA returned",
            call. = FALSE)
  }
  out
}

.fcr_wtable_logpmf <- function(y, pars, logw, k0fun) {
  .fcr_by_table(pars, logw, k0fun, function(idx, tab, pars1) {
    do.call(logw, c(list(y[idx]), pars1)) - tab$logZ
  })
}

.fcr_wtable_logcdf <- function(y, pars, lower, logw, k0fun) {
  .fcr_by_table(pars, logw, k0fun, function(idx, tab, pars1) {
    yi <- y[idx]
    out <- numeric(length(yi))
    inside <- yi <= tab$K - 1
    if (any(inside)) {
      out[inside] <- if (lower) tab$cum[yi[inside] + 1] - tab$logZ else
        tab$rcum[yi[inside] + 2] - tab$logZ
    }
    for (i in which(!inside)) {
      # Past the table: sum the (decreasing) upper tail directly
      ls <- .fcr_log_upper_tail(yi[i], function(x)
        do.call(logw, c(list(x), pars1)) - tab$logZ)
      out[i] <- if (lower) .fcr_log1mexp(ls) else ls
    }
    out
  })
}

.fcr_wtable_moments <- function(pars, logw, k0fun) {
  m1 <- .fcr_by_table(pars, logw, k0fun, function(idx, tab, p1) {
    rep(sum((0:tab$K) * exp(tab$lw - tab$logZ)), length(idx))
  })
  m2 <- .fcr_by_table(pars, logw, k0fun, function(idx, tab, p1) {
    rep(sum((0:tab$K)^2 * exp(tab$lw - tab$logZ)), length(idx))
  })
  list(mean = m1, var = m2 - m1^2)
}
