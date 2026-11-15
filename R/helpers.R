# Helper functions for the flexCountReg package

#' @importFrom stats model.frame model.matrix model.response
#' @importFrom stats dnbinom dpois qnorm qgamma dnorm
#' @importFrom randtoolbox halton
NULL

# ---- Internal helpers for estimators ----------------------------------------

# Families supported by countreg() (names used by get_params/get_probFunc)
.fcr_families <- c("POISSON", "GP1", "GP2", "NB1", "NB2", "NBP", "PLN", "PGE",
                   "PIG1", "PIG2", "PIG", "PL", "PLG", "PLL", "PW", "SI", "GW",
                   "COM")

# Validate and normalize the optimizer name. An unknown method is an error.
.fcr_check_method <- function(method) {
  valid <- c("SANN", "NM", "BFGS", "BFGSR", "CG", "NR", "BHHH")
  m <- toupper(gsub("[^[:alnum:]]", "", method))
  if (length(m) != 1L || !m %in% valid) {
    stop("`method` must be one of: ", toString(paste0('"', valid, '"')), ".",
         call. = FALSE)
  }
  m
}

# Validate user starting values: right length, matched by name when named.
.fcr_check_start <- function(start.vals, x_names) {
  start.vals <- unlist(start.vals)
  if (!is.numeric(start.vals)) {
    stop("`start.vals` must be numeric.", call. = FALSE)
  }
  if (length(start.vals) != length(x_names)) {
    stop("`start.vals` has ", length(start.vals), " value(s), but the model ",
         "has ", length(x_names), " parameter(s): ", toString(x_names), ".",
         call. = FALSE)
  }
  nm <- names(start.vals)
  if (!is.null(nm) && all(nzchar(nm))) {
    if (!setequal(nm, x_names)) {
      stop("Names of `start.vals` do not match the model parameters. ",
           "Expected: ", toString(x_names), ".", call. = FALSE)
    }
    start.vals <- start.vals[x_names]
  }
  names(start.vals) <- x_names
  start.vals
}

# Probability of an event being reported, given the underreporting linear
# predictor. Both links use the same sign convention: a positive coefficient
# means the event is more likely to be reported.
.report_prob <- function(lin, family = c("logit", "probit")) {
  switch(match.arg(family),
         logit = stats::plogis(lin),
         probit = stats::pnorm(lin))
}

# Did maxLik converge? Return codes differ by optimizer.
.fcr_converged <- function(fit, method) {
  code <- fit$code
  if (is.null(code)) return(NA)
  if (method %in% c("NR", "BHHH", "BFGSR")) {
    code %in% c(1L, 2L, 8L)
  } else {
    code == 0L
  }
}

# Covariance from the Hessian, with a generalized-inverse fallback.
.fcr_hessian_vcov <- function(fit) {
  H <- fit$hessian
  V <- tryCatch(solve(-H), error = function(e) NULL)
  if (is.null(V)) {
    warning("The Hessian is singular; standard errors use a generalized ",
            "inverse and may be unreliable. The model may not be identified.",
            call. = FALSE)
    V <- MASS::ginv(-H)
  }
  V
}

# Nonparametric bootstrap of a set of estimates.
#   n          number of rows (or clusters) to resample
#   B          number of replicates
#   coef_names names of the estimates, in order
#   refit      function(idx) returning a named numeric vector of estimates
# Returns a B x length(coef_names) matrix; failed replicates are NA rows.
.fcr_bootstrap <- function(n, B, coef_names, refit) {
  draws <- matrix(NA_real_, nrow = B, ncol = length(coef_names),
                  dimnames = list(NULL, coef_names))
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    est <- tryCatch(suppressWarnings(suppressMessages(refit(idx))),
                    error = function(e) NULL)
    if (!is.null(est) && all(coef_names %in% names(est))) {
      draws[b, ] <- est[coef_names]
    }
  }
  n_ok <- sum(stats::complete.cases(draws))
  if (n_ok < 2) {
    stop("Fewer than 2 bootstrap replicates succeeded; bootstrapped ",
         "standard errors cannot be computed.", call. = FALSE)
  }
  if (n_ok < 0.9 * B) {
    warning(B - n_ok, " of ", B, " bootstrap replicates failed; the ",
            "bootstrapped standard errors use the ", n_ok, " that succeeded.",
            call. = FALSE)
  }
  draws
}

#' Tidy a flexCountReg object
#'
#' @param x An object of class flexCountReg
#' @param ... Additional arguments
#' @returns A tibble object with columns: term, estimate
#' @importFrom tibble tibble
#' @export
tidy.flexCountReg <- function(x, ...) {
  # Attempt to locate coefficients based on standard countreg structure
  # Check if coefficients are at the top level or inside a 'model' slot
  if (!is.null(x$coefficients)) {
    ests <- x$coefficients
  } else if (!is.null(x$model) && !is.null(x$model$coefficients)) {
    ests <- x$model$coefficients
  } else {
    stop("Could not find coefficients in flexCountReg object for tidying.")
  }
  
  # Return the tibble structure required by broom
  tibble::tibble(
    term = names(ests),
    estimate = as.numeric(ests)
  )
}


# Determine model to estimate, probability distribution to use, and parameters
get_params <- function(family) {
  switch( # pattern: list("alpha", "sigma")
    family,
    "POISSON" = list(NULL, NULL),
    "GP1" = list("ln(phi+1)", NULL),
    "GP2" = list("ln(alpha+1)", NULL),
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
    "GW" = list("ln(k)", "ln(rho-1)"),
    "COM" = list("ln(nu)", NULL)
  )
}

get_probFunc <- function(family){
  switch( # get the probability function for the specified distribution
    family,
    "POISSON" = function(y, predicted, alpha, sigma, haltons, normed_haltons, engine) {
      return(stats::dpois(y, predicted))
    },
    "GP1"=function(y, predicted, alpha, sigma, haltons, normed_haltons, engine){
      return(dgp1(y, mu=predicted, phi=(alpha-1)))
    },
    "GP2"=function(y, predicted, alpha, sigma, haltons, normed_haltons, engine){
      return(dgp2(y, mu=predicted, alpha=(alpha-1)))
    },
    "NB1" = function(y, predicted, alpha, sigma, haltons, normed_haltons, engine) {
      return(stats::dnbinom(y, size = predicted/alpha, mu = predicted))
    },
    "NB2" = function(y, predicted, alpha, sigma, haltons, normed_haltons, engine) {
      return(stats::dnbinom(y, size = 1/alpha, mu = predicted))
    },
    "NBP" = function(y, predicted, alpha, sigma, haltons, normed_haltons, engine){
      return(
        stats::dnbinom(y, size = (predicted^(2-sigma))/alpha, mu = predicted))
    },
    "PLN" = function(y, predicted, alpha, sigma, haltons, normed_haltons, engine) 
      dpLnorm(x = y, mean = predicted, sigma = alpha, hdraws = haltons, engine=engine),
    "PGE" = function(y, predicted, alpha, sigma, haltons, normed_haltons, engine) 
      dpge(y, mean = predicted, shape = alpha, scale = sigma,haltons = haltons),
    "PIG1" = function(y, predicted, alpha, sigma, ...) 
      dpinvgaus(y, mu = predicted, eta = alpha),
    "PIG2" = function(y, predicted, alpha, sigma, ...) 
      dpinvgaus(y, mu=predicted, eta = alpha, form = "Type 2"),
    "PIG" = function(y, predicted, alpha, sigma, ...) 
      dpinvgamma(y, mu = predicted, eta = alpha),
    "PL" = function(y, predicted, alpha, sigma, ...) 
      dplind(y, mean = predicted, theta = alpha),
    "PLG" = function(y, predicted, alpha, sigma, ...) 
      dplindGamma(x = y, mean = predicted, theta = alpha, alpha = sigma),
    "PLL" = function(y, predicted, alpha, sigma, haltons, normed_haltons, engine) 
      dplindLnorm(x = y, mean = predicted, theta = alpha, sigma = sigma, 
                  hdraws = haltons),
    "PW" = function(y, predicted, alpha, sigma, haltons, normed_haltons, engine) 
      dpWeib_cpp(y, mean = predicted, alpha = alpha, sigma = sigma, h =haltons),
    "SI" = function(y, predicted, alpha, sigma, ...) 
      dsichel(x = y, mu = predicted, sigma = sigma, gamma = log(alpha)),
    "GW" = function(y, predicted, alpha, sigma, ..., log = FALSE)
      dgwar(y, mu = predicted, k = alpha, rho = 1 + sigma,
            log = log),
    "COM" = function(y, predicted, alpha, sigma, ...) 
      dcom(x = y, mu = predicted, nu = alpha)
  )
}

# --- Internal Likelihood Functions ---

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
    
    mu_matrix <- 
      matrix(random_coefs_means, nrow = n_obs, ncol = n_params, byrow = TRUE)
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
          param_draws[, d] <- 
            stats::qgamma(h_col[d], shape = shape_vec, rate = rate_vec)
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