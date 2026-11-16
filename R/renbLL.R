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
  
  # Calculate log-likelihood using lgamma for numerical stability
  # lgamma(a + sum_mu) - lgamma(a + b + sum_mu + sum_y), computed without
  # cancellation when sum_mu is large
  log_P <- .fcr_lgamma_diff(a, b) + lgamma(b + df_model$sum_y) - lgamma(b) -
    .fcr_lgamma_diff(a + df_model$sum_mu, b + df_model$sum_y) +
    df_model$log_prod_element
  
  # Sum across panels for total log-likelihood
  total_loglik <- log_P
  
  return(total_loglik)
}

# lgamma(x + d) - lgamma(x) for x > 0, d >= 0, accurate when x is much larger
# than d (where the direct difference cancels to zero in floating point).
.fcr_lgamma_diff <- function(x, d) {
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
