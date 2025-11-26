# Helper functions for the flexCountReg package

#' @importFrom stats model.frame model.matrix model.response dnbinom dpois qnorm qgamma dnorm
#' @importFrom randtoolbox halton
NULL

# create function to clean data and run countreg for bootstrapping
# create function to clean data and run countreg for bootstrapping
mod.boot <- function(data, formula, family, offset, weights, 
                     dis_param_formula_1, dis_param_formula_2,
                     underreport_formula, underreport_family,
                     ndraws, method, max.iters, start.vals){
  
  # data comes in as a resample object from modelr, convert to df
  df <- as.data.frame(data)
  
  # Recursively call countreg with the bootstrap sample
  # We set stderr = "none" to prevent infinite recursion
  # We use the provided start.vals (from the original fit) to speed up convergence
  fit <- tryCatch({
    countreg(formula = formula, 
             data = df, 
             family = family, 
             offset = offset, 
             weights = weights, 
             verbose = FALSE, 
             dis_param_formula_1 = dis_param_formula_1, 
             dis_param_formula_2 = dis_param_formula_2, 
             underreport_formula = underreport_formula,
             underreport_family = underreport_family,
             ndraws = ndraws, 
             method = method, 
             max.iters = max.iters, 
             start.vals = start.vals, 
             stderr = "none", 
             bootstraps = NULL)
  }, error = function(e) return(NULL)) # Handle failures in bootstrap samples
  
  return(fit)
}


# Determine model to estimate, probability distribution to use, and parameters
get_params <- function(family) {
  switch(
    family,
    "POISSON" = list(NULL, NULL),
    "NB1" = list("ln(alpha)", NULL),
    "NB2" = list("ln(alpha)", NULL),
    "NBP" = list("ln(alpha)", "ln(p)"),
    "PLN" = list("ln(sigma)", NULL),
    "PGE" = list("ln(shape)", "ln(scale)"),
    "PIG1" = list("ln(eta)", NULL),
    "PIG2" = list("ln(eta)", NULL),
    "PIG" = list("ln(eta)", NULL),
    "PL" = list("ln(theta)", NULL),
    "PLG" = list("ln(theta)", "ln(alpha)"),
    "PLL" = list("ln(theta)", "ln(sigma)"),
    "PW" = list("ln(alpha)", "ln(sigma)"),
    "SI" = list("gamma", "ln(sigma)"),
    "GW" = list("ln(k)", "ln(rho)"),
    "COM" = list("ln(nu)", NULL),
  )
}

get_probFunc <- function(family){
  switch( # get the probability function for the specified distribution
    family,
    "POISSON" = function(y, predicted, alpha, sigma, haltons, normed_haltons) {
      return(stats::dpois(y, predicted))
    },
    "NB1" = function(y, predicted, alpha, sigma, haltons, normed_haltons) {
      return(stats::dnbinom(y, size = predicted/alpha, mu = predicted))
    },
    "NB2" = function(y, predicted, alpha, sigma, haltons, normed_haltons) {
      return(stats::dnbinom(y, size = 1/alpha, mu = predicted))
    },
    "NBP" = function(y, predicted, alpha, sigma, haltons, normed_haltons){
      return(stats::dnbinom(y, size = (predicted^(2-sigma))/alpha, mu = predicted))
    },
    "PLN" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpLnorm_cpp(x=y, mean=predicted, sigma=alpha, h=normed_haltons),
    "PGE" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpge(y, mean=predicted, shape=alpha, scale=sigma, haltons=haltons),
    "PIG1" = function(y, predicted, alpha, sigma, ...) dpinvgaus(y, mu=predicted, eta=alpha),
    "PIG2" = function(y, predicted, alpha, sigma, ...) dpinvgaus(y, mu=predicted, eta=alpha, form="Type 2"),
    "PIG" = function(y, predicted, alpha, sigma, ...) dpinvgamma(y, mu=predicted, eta=alpha),
    "PL" = function(y, predicted, alpha, sigma, ...) dplind(y, mean=predicted, theta=alpha),
    "PLG" = function(y, predicted, alpha, sigma, ...) dplindGamma(x=y, mean=predicted, theta=alpha, alpha=sigma),
    "PLL" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dplindLnorm(x=y, mean=predicted, theta=alpha, sigma=sigma, hdraws=normed_haltons),
    "PW" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpWeib_cpp(y, mean=predicted, alpha=alpha, sigma=sigma, h=haltons),
    "SI" = function(y, predicted, alpha, sigma, ...) dsichel(x=y, mu= predicted, sigma=sigma, gamma=log(alpha)),
    "GW" = function(y, predicted, alpha, sigma, ...) dgwar(y, mu= predicted, k=alpha, rho=sigma),
    "COM" = function(y, predicted, alpha, sigma, ...) dcom(x=y, mu=predicted, nu=alpha)
  )
}

# --- Internal Likelihood Functions ---

nb_prob <- function(y, mu, alpha, p = NULL, form="nb2") {
  if (form == 'nb2') {
    return(stats::dnbinom(y, size = 1/alpha, mu = mu))
  } else if (form == 'nb1') {
    return(stats::dnbinom(y, size = mu / alpha, mu = mu))
  } else if (form == 'nbp' && !is.null(p)) {
    return(stats::dnbinom(y, size = (mu^(2 - p)) / alpha, mu = mu))
  } else {
    stop("Invalid form or missing 'p'.")
  }
}

p_nb_rp <- function(p, y, X_Fixed, X_rand, ndraws, rpar, correlated, form, 
                    rpardists, hdraws, data, weights, X_offset = NULL, offset = NULL,
                    X_het_mean = NULL, X_het_var = NULL) {
  
  coef_info <- extract_coefficients(p, ncol(X_Fixed), length(rpar), form, 
                                    if(is.null(X_het_mean)) 0 else ncol(X_het_mean),
                                    if(is.null(X_het_var)) 0 else ncol(X_het_var),
                                    correlated)
  
  mu_fixed <- compute_fixed_effects(X_Fixed, coef_info$fixed_coefs, X_offset, offset)
  
  draws_info <- generate_random_draws(hdraws = hdraws, 
                                      random_coefs_means = coef_info$random_coefs_means, 
                                      rand_var_params = coef_info$rand_var_params, 
                                      rpardists = rpardists, 
                                      rpar = rpar, 
                                      X_rand = X_rand,
                                      het_mean_coefs = coef_info$het_mean_coefs, 
                                      X_het_mean = X_het_mean,
                                      het_var_coefs = coef_info$het_var_coefs, 
                                      X_het_var = X_het_var,
                                      correlated = correlated)
  
  log_probs <- compute_log_likelihoods(draws_info, mu_fixed, y, coef_info$alpha, 
                                       coef_info$p, form, weights, data$panel_id)
  
  return(log_probs)
}

extract_coefficients <- function(p, N_fixed, N_rand, form, N_het_mean, N_het_var, correlated) {
  coefs <- as.array(p)
  current_idx <- 0
  
  fixed_coefs <- coefs[(current_idx + 1):(current_idx + N_fixed)]
  current_idx <- current_idx + N_fixed
  
  random_coefs_means <- coefs[(current_idx + 1):(current_idx + N_rand)]
  current_idx <- current_idx + N_rand
  
  n_var_params <- if(correlated) N_rand * (N_rand + 1) / 2 else N_rand
  rand_var_params <- coefs[(current_idx + 1):(current_idx + n_var_params)]
  current_idx <- current_idx + n_var_params
  
  het_mean_coefs <- NULL
  if (N_het_mean > 0) {
    het_mean_coefs <- coefs[(current_idx + 1):(current_idx + N_het_mean)]
    current_idx <- current_idx + N_het_mean
  }
  
  het_var_coefs <- NULL
  if (N_het_var > 0) {
    het_var_coefs <- coefs[(current_idx + 1):(current_idx + N_het_var)]
    current_idx <- current_idx + N_het_var
  }
  
  log_alpha <- coefs[current_idx + 1]
  p_param <- if(form == 'nbp') coefs[current_idx + 2] else NULL
  
  return(list(
    fixed_coefs = fixed_coefs,
    random_coefs_means = random_coefs_means,
    rand_var_params = rand_var_params,
    het_mean_coefs = het_mean_coefs,
    het_var_coefs = het_var_coefs,
    alpha = exp(log_alpha),
    p = p_param
  ))
}

compute_fixed_effects <- function(X_Fixed, fixed_coefs, X_offset, offset) {
  linear_pred <- X_Fixed %*% fixed_coefs
  if (!is.null(offset)) {
    if (length(offset) > 1) linear_pred <- linear_pred + rowSums(X_offset)
    else linear_pred <- linear_pred + offset
  }
  return(exp(linear_pred))
}

generate_random_draws <- function(hdraws, random_coefs_means, rand_var_params, 
                                  rpardists, rpar, X_rand,
                                  het_mean_coefs = NULL, X_het_mean = NULL,
                                  het_var_coefs = NULL, X_het_var = NULL,
                                  correlated = FALSE) {
  
  n_obs <- nrow(X_rand)
  n_draws <- nrow(hdraws)
  n_params <- length(rpar)
  
  # 1. Calculate Heterogeneity Multipliers
  mean_multiplier <- rep(1, n_obs)
  var_multiplier <- rep(1, n_obs)
  
  if (!is.null(X_het_mean) && !is.null(het_mean_coefs)) {
    mean_multiplier <- exp(as.vector(X_het_mean %*% het_mean_coefs))
  }
  if (!is.null(X_het_var) && !is.null(het_var_coefs)) {
    var_multiplier <- exp(as.vector(X_het_var %*% het_var_coefs))
  }
  
  # 2. Calculate Random Effects
  if (correlated) {
    # --- Correlated (Normal Only) ---
    L <- matrix(0, n_params, n_params)
    L[lower.tri(L, diag = TRUE)] <- rand_var_params
    
    z_draws <- stats::qnorm(hdraws) 
    corr_draws <- z_draws %*% t(L) 
    
    mu_matrix <- matrix(random_coefs_means, nrow=n_obs, ncol=n_params, byrow=TRUE)
    xb_deterministic <- rowSums(X_rand * mu_matrix * mean_multiplier)
    
    W_variance <- X_rand * var_multiplier
    xb_stochastic <- W_variance %*% t(corr_draws)
    
    xb_rand_mat <- xb_stochastic + xb_deterministic
    return(list(xb_rand_mat = xb_rand_mat))
    
  } else {
    # --- Uncorrelated (Supports Distributions) ---
    xb_rand_mat <- matrix(0, nrow = n_obs, ncol = n_draws)
    
    for (k in 1:n_params) {
      mu_base <- random_coefs_means[k]
      sd_base <- abs(rand_var_params[k])
      
      mu_i <- mu_base * mean_multiplier
      sd_i <- sd_base * var_multiplier
      
      dist_type <- if (is.null(rpardists)) "n" else rpardists[k]
      
      if (dist_type == "g") {
        # Gamma: Recalculate shape/rate per observation
        shape_vec <- (mu_i^2) / (sd_i^2)
        rate_vec <- mu_i / (sd_i^2)
        h_col <- hdraws[, k]
        
        param_draws <- matrix(0, nrow = n_obs, ncol = n_draws)
        for (d in 1:n_draws) {
          param_draws[, d] <- stats::qgamma(h_col[d], shape = shape_vec, rate = rate_vec)
        }
        xb_rand_mat <- xb_rand_mat + (X_rand[, k] * param_draws)
        
      } else {
        # Location-Scale Families
        if (dist_type == "ln") {
          z_scores <- stats::qnorm(hdraws[, k])
          # exp(mu_i + sd_i * Z)
          exponent_mat <- outer(mu_i, rep(1, n_draws)) + outer(sd_i, z_scores)
          param_draws <- exp(exponent_mat)
          xb_rand_mat <- xb_rand_mat + (X_rand[, k] * param_draws)
        } else {
          # Normal, Triangle, Uniform
          std_draws <- switch(dist_type,
                              "n" = stats::qnorm(hdraws[, k]),
                              "t" = qtri(hdraws[, k], 0, 1),
                              "u" = hdraws[, k] - 0.5,
                              stats::qnorm(hdraws[, k]))
          
          term1 <- X_rand[, k] * mu_i 
          term2_scale <- X_rand[, k] * sd_i 
          
          stochastic_part <- outer(term2_scale, std_draws)
          xb_rand_mat <- xb_rand_mat + stochastic_part + term1 
        }
      }
    }
    return(list(xb_rand_mat = xb_rand_mat))
  }
}

compute_log_likelihoods <- function(draws_info, mu_fixed, y, alpha, p, form, 
                                    weights, panel_id) {
  rpar_mat <- exp(draws_info$xb_rand_mat)
  pred_mat <- sweep(rpar_mat, 1, mu_fixed, "*")
  prob_mat <- apply(pred_mat, 2, nb_prob, y = y, alpha = alpha, p = p, form = form)
  prob_mat <- prob_mat^weights
  
  log_prob_mat <- log(prob_mat)
  log_prob_df <- as.data.frame(log_prob_mat)
  log_prob_df$panel_id <- panel_id
  
  log_probs <- aggregate(. ~ panel_id, data = log_prob_df, FUN = sum)
  log_probs$panel_id <- NULL
  
  log_probs_mat <- as.matrix(log_probs)
  probs <- rowMeans(exp(log_probs_mat))
  
  return(log(probs))
}