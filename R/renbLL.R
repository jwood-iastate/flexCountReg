#' @importFrom dplyr group_by reframe %>%
#' @importFrom tibble tibble
#'
renb_ll <- function(y, mu, a, b, panels) { 
  # Random Effects Negative Binomial with Beta Distributed Random Effects (NB1)
  
  
  if (any(!is.finite(mu)) || any(mu <= 0) ||
      !is.finite(a) || a <= 0 ||
      !is.finite(b) || b <= 0) {
    return(rep(-Inf, length(unique(panels))))
  }
  
  df <- tibble(y = y, mu = mu, panels = panels)
  
  # Calculate panel-level statistics
  df_model <- df %>%
    group_by(panels) %>%
    reframe(
      sum_mu = sum(mu),
      sum_y = sum(y),
      # Using lgamma instead of gamma for numerical stability
      # Note: y+1 in gamma(y+1) accounts for factorial
      log_prod_element = sum(.fcr_lgamma_diff(mu, y) - lgamma(y + 1))
    )
  
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
