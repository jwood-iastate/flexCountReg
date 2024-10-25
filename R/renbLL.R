renb_ll <- function(y, mu, a, b, panels) { # Random Effects Negative Binomial with Beta Distributed Random Effects (NB1)
  require(dplyr)
  require(tibble)
  
  df <- tibble(y = y, mu = mu, panels = panels)
  
  # Calculate panel-level statistics
  df_model <- df %>%
    group_by(panels) %>%
    reframe(
      sum_mu = sum(mu),
      sum_y = sum(y),
      # Using lgamma instead of gamma for numerical stability
      # Note: y+1 in gamma(y+1) accounts for factorial
      log_prod_element = sum(lgamma(mu + y) - lgamma(mu) - lgamma(y + 1))
    )
  
  # Calculate log-likelihood using lgamma for numerical stability
  log_num <- lgamma(a + b) + lgamma(a + df_model$sum_mu) + lgamma(b + df_model$sum_y)
  
  log_denom <- lgamma(a) + lgamma(b) + 
    lgamma(a + b + df_model$sum_mu + df_model$sum_y)
  
  # Final log-likelihood calculation
  log_P <- log_num - log_denom + df_model$log_prod_element
  
  # Sum across panels for total log-likelihood
  total_loglik <- sum(log_P)
  
  return(total_loglik)
}