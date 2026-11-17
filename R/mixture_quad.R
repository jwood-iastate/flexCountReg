# Mixed Poisson distributions computed by adaptive quadrature
# (src/mixtureQuad.cpp): Poisson-lognormal (PLN), Poisson-generalized
# exponential (PGE), Poisson-Weibull (PW) and Poisson-Lindley-lognormal
# (PLL).
#
# The log PMF is computed directly on the log scale: the log integrand is
# maximized over the (log) mixing variable and integrated with a
# double-exponential rule centered at its mode, so counts far in the tail
# keep their relative accuracy. The rule is deterministic, so the
# log-likelihood is a smooth function of the parameters.
#
# Earlier versions averaged the conditional probabilities over a fixed set
# of Halton draws (engine = "halton", kept to reproduce earlier fits). That
# average has no draws where the integrand has its mass for counts far in
# the tail, so it is badly biased there even when it is finite: for PGE with
# shape 1 (a geometric distribution with mean 1) and y = 1000, 1500 draws
# give log P = -3847 instead of the exact -693.8.

.fcr_mix_kind <- c(PLN = 0L, PGE = 1L, PW = 2L, PLL = 3L)

# log P(Y = x) of a mixture family, recycled over the arguments. `p1` is
# sigma (PLN), shape (PGE), alpha (PW) or theta (PLL); `p2` is sigma (PLL).
# Invalid parameters give NaN with a warning; missing values give NA.
.fcr_mix_logpmf <- function(x, mean, p1, p2 = 0, family) {
  r <- .fcr_recycle(x, mean, p1, p2)
  if (length(r[[1]]) == 0L) return(numeric(0))
  out <- mix_logpmf_cpp(as.numeric(r[[1]]), as.numeric(r[[2]]),
                        as.numeric(r[[3]]), as.numeric(r[[4]]),
                        .fcr_mix_kind[[family]])
  .fcr_warn_nan(is.nan(out) & !is.na(r[[1]]))
  out
}

# P(Y <= q) (or P(Y > q)) for finite q >= 0, by summing the PMF over
# 0, ..., floor(q) once per distinct parameter set. The p-function contract
# (.fcr_p_contract()) handles the endpoints, NA and small upper tails.
.fcr_mix_cdf <- function(q, mean, p1, p2 = 0, family, lower.tail = TRUE,
                         log.p = FALSE) {
  r <- .fcr_recycle(q, mean, p1, p2)
  q <- floor(r[[1]]); mean <- r[[2]]; p1 <- r[[3]]; p2 <- r[[4]]
  n <- length(q)
  out <- rep(NA_real_, n)
  key <- paste(format(mean, digits = 17), format(p1, digits = 17),
               format(p2, digits = 17), sep = "\r")
  for (k in unique(key)) {
    idx <- which(key == k)
    lp <- suppressWarnings(.fcr_mix_logpmf(0:max(q[idx]), mean[idx[1]],
                                           p1[idx[1]], p2[idx[1]], family))
    if (anyNA(lp)) {
      out[idx] <- lp[1]   # NaN for invalid parameters
      next
    }
    m <- max(lp)
    lcdf <- m + log(cumsum(exp(lp - m)))
    lcdf <- pmin(lcdf[q[idx] + 1], 0)
    out[idx] <- if (lower.tail) {
      if (log.p) lcdf else exp(lcdf)
    } else {
      if (log.p) .fcr_log1mexp(lcdf) else -expm1(lcdf)
    }
  }
  .fcr_warn_nan(is.nan(out))
  out
}

# The numerical engine of a mixture family: adaptive quadrature (default),
# the earlier Halton-draw average, or the poilog package (PLN only).
.fcr_mix_engine <- function(engine, poilog = FALSE) {
  choices <- c("quadrature", "halton", if (poilog) "poilog")
  if (is.null(engine)) return("quadrature")
  match.arg(engine[1], choices)
}

# Engine for the families without a poilog option (PGE, PW, PLL): "halton"
# keeps the earlier draws, anything else is quadrature.
.fcr_mix_engine_nopl <- function(engine) {
  if (identical(engine, "halton")) "halton" else "quadrature"
}
