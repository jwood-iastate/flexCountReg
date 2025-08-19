#' Random Parameters Negative Binomial Model
#'
#' @description Estimates a random parameters count model for any of the distributions also used in `countreg`
#'
#' @name rpcount
#' @param formula An R formula for the fixed part of the model.
#' @param rpar_formula A one-sided formula for the random parameters (e.g., `~ var1 + var2`).
#' @param data A data frame containing all model variables.
#' @param family The model distribution type (e.g., "NB1", "NB2", "NBP").
#' @param rpardists A named character vector specifying the distribution for each random
#'   parameter. Options: "n" (normal), "ln" (log-normal), "t" (triangular), "u" (uniform), "g" (gamma).
#'   Defaults to normal for all.
#' @param ndraws The number of Halton draws for simulation.
#' @param scrambled Logical. If `TRUE`, uses scrambled Halton draws.
#' @param correlated Logical. If `TRUE`, estimates correlated random parameters.
#'   Only the normal distribution is supported for correlated parameters.
#' @param panel A character vector of variable names defining the panel structure.
#' @param offset A character vector of variable names to be used as an offset.
#' @param weights A string specifying the variable name for frequency weights.
#' @param method The optimization method for `maxLik`. See `?maxLik::maxLik`.
#' @param max.iters Maximum iterations for the optimizer.
#' @param start.vals An optional named vector of starting values.
#' @param print.level Verbosity level for the optimizer (0, 1, or 2).
#'
#' @import randtoolbox stats modelr
#' @importFrom maxLik maxLik
#' @importFrom MASS glm.nb
#' @include tri.R get_chol.R helpers.R
#'
#' @examples
#' \donttest{
#' data("washington_roads")
#'
#' ## Uncorrelated Random Parameters Poisson-Lognormal Model
#' nb2.rp <- rpcount(
#'   Total_crashes ~ lnlength + lnaadt,
#'   rpar_formula = ~ -1 + speed50,
#'   data = washington_roads,
#'   ndraws = 50, # Using fewer draws for a quick example
#'   rpardists = c(speed50 = "n"),
#'   family = 'PLN',
#'   method = "bfgs"
#' )
#' summary(nb2.rp)
#'
#' ## Correlated Random Parameters NB1 Model
#' nb1.rp.corr <- rpcount(
#'   Total_crashes ~ lnlength,
#'   rpar_formula = ~ lnaadt + speed50,
#'   data = washington_roads,
#'   ndraws = 50,
#'   correlated = TRUE,
#'   family = 'nb1',
#'   method = "bfgs"
#' )
#' summary(nb1.rp.corr)
#' }
#' @export
rpcount <- function(formula, rpar_formula, data, family = 'nb2',
                 rpardists = NULL,
                 ndraws = 1500, scrambled = FALSE,
                 correlated = FALSE, panel = NULL,
                 weights = NULL, offset = NULL,
                 method = 'BHHH', max.iters = 1000,
                 start.vals = NULL, print.level = 0) {
  
  # --- 1. Initial Setup and Data Preparation ---
  data <- as.data.frame(data)
  
  # Generate a unique panel ID for each observation or group
  if (is.null(panel)) {
    # If no panel structure is specified, each observation is its own panel
    panel_id <- 1:nrow(data)
  } else {
    # Combine multiple panel variables into a single identifier
    panel_id <- apply(data[, panel, drop = FALSE], 1, paste, collapse = "_")
  }
  # Store panel_id in the data frame for use in the likelihood function
  data$panel_id <- as.factor(panel_id)
  
  # Prepare model matrices
  X_Fixed <- stats::model.matrix(formula, data)
  X_rand <- stats::model.matrix(rpar_formula, data)
  y <- stats::model.response(stats::model.frame(formula, data))
  
  # Validate that intercept is not in both fixed and random parts
  if ("(Intercept)" %in% colnames(X_Fixed) && "(Intercept)" %in% colnames(X_rand)) {
    stop("Intercept cannot be in both `formula` and `rpar_formula`. Please remove one.")
  }
  
  # Prepare weights and offset vectors
  weights_vec <- if (is.null(weights)) rep(1, nrow(data)) else data[[weights]]
  offset_vec <- if (!is.null(offset)) rowSums(data[, offset, drop = FALSE]) else NULL
  
  # Validate and prepare random parameter distributions
  rand_terms <- colnames(X_rand)
  if (correlated) {
    if (!is.null(rpardists) && any(rpardists != "n")) {
      warning("For correlated parameters, only the normal distribution is used. Ignoring `rpardists`.")
    }
    rpardists <- setNames(rep("n", length(rand_terms)), rand_terms)
  } else {
    if (is.null(rpardists)) {
      rpardists <- setNames(rep("n", length(rand_terms)), rand_terms)
    } else {
      # Ensure all random terms have a specified distribution
      if (!all(rand_terms %in% names(rpardists))) {
        stop("Not all random parameters in `rpar_formula` have a distribution specified in `rpardists`.")
      }
    }
  }
  # Ensure order is correct
  rpardists <- rpardists[rand_terms]
  
  # Generate Halton draws
  hdraws <- randtoolbox::halton(ndraws, length(rand_terms), mixed = scrambled)
  haltons <- randtoolbox::halton(nrow(hdraws), (ncol(hdraws)+1), mixed = TRUE)[,(ncol(hdraws)+1)]
  normed_haltons <- randtoolbox::halton(nrow(hdraws), (ncol(hdraws)+1), mixed = TRUE, normal=TRUE)[,(ncol(hdraws)+1)]
  
  # --- 2. Starting Values ---
  if (is.null(start.vals)) {
    # Combine formulas for initial model fitting
    full_formula <- update(formula, paste(". ~ . +", paste(setdiff(colnames(X_rand), "(Intercept)"), collapse = " + ")))
    
    # Get starting values from a standard NB model, with fallback to Poisson
    nb_model <- tryCatch({
      MASS::glm.nb(full_formula, data = data)
    }, error = function(e) {
      warning("glm.nb failed for starting values, falling back to poisson model. Error: ", e$message)
      glm(full_formula, data = data, family = "poisson")
    })
    
    # Extract coefficients and align them
    model_coefs <- coef(nb_model)
    fixed_coef_names <- colnames(X_Fixed)
    rand_coef_names <- colnames(X_rand)
    
    start_fixed <- model_coefs[names(model_coefs) %in% fixed_coef_names]
    start_rand_means <- model_coefs[names(model_coefs) %in% rand_coef_names]
    
    # Handle cases where some coefficients might not be estimated (e.g., collinearity)
    if(length(start_fixed) < length(fixed_coef_names)) start_fixed <- setNames(rep(0, length(fixed_coef_names)), fixed_coef_names)
    if(length(start_rand_means) < length(rand_coef_names)) start_rand_means <- setNames(rep(0, length(rand_coef_names)), rand_coef_names)
    
    start <- c(start_fixed, start_rand_means)
    
    # Starting values for standard deviations or Cholesky
    if (correlated) {
      # For correlated, start with a diagonal covariance matrix
      chol_start_mat <- diag(0.1, length(rand_terms))
      chol_start_vec <- chol_start_mat[upper.tri(chol_start_mat, diag = TRUE)]
      start <- c(start, chol_start_vec)
    } else {
      start <- c(start, rep(0.1, length(rand_terms))) # St. Devs
    }
    
    # Add ancillary parameters 
    # Get number of distribution params based on family
    if (toupper(family) == 'NB2') {
      start <- c(start, 0) # ln(alpha)
      parnames <- "ln(alpha)"
    } else if (toupper(family) == 'NB1') {
      start <- c(start, 0) # ln(alpha)
      parnames <- "ln(alpha)"
    } else if (toupper(family) == 'NBP') {
    start <- c(start, 0, log(1.5)) # ln(alpha)+
    parnames <- c("ln(alpha)", "ln(P)")
    } else if (toupper(family) == 'PLN') {
      start <- c(start, 0) # ln(alpha)+
      parnames <- "ln(sigma)"
    }else if (toupper(family) == 'PGE') {
      start <- c(start, 0, 0) # 
      parnames <- c("ln(shape)", "ln(scale)")
    }else if (toupper(family) == 'PIG1') {
      start <- c(start, 0) # 
      parnames <- "ln(eta)"
    }else if (toupper(family) == 'PIG2') {
      start <- c(start, 0) # 
      parnames <- "ln(eta)"
    }else if (toupper(family) == 'PIG') {
      start <- c(start, 0) # 
      parnames <- "ln(eta)"
    }else if (toupper(family) == 'PL') {
      start <- c(start, 0) # 
      parnames <- "ln(theta)"
    }else if (toupper(family) == 'PLL') {
      start <- c(start, 0, 0) # 
      parnames <- c("ln(theta)","ln(sigma)")
    }else if (toupper(family) == 'PLG') {
      start <- c(start, 0, 0) # 
      parnames <- c("ln(theta)","ln(alpha)")
    }else if (toupper(family) == 'PW') {
      start <- c(start, 0, 0) # 
      parnames <- c("ln(alpha)","ln(sigma)")
    }else if (toupper(family) == 'SI') {
      start <- c(start, 0, 0) # 
      parnames <- c("ln(gamma)","ln(sigma)")
    }else if (toupper(family) == 'GW') {
      start <- c(start, 0, 0) # 
      parnames <- c("ln(k)","ln(rho)")
    }else if (toupper(family) == 'COM') {
      start <- c(start, 0) # 
      parnames <- "ln(nu)"
    }
    
    # Create descriptive names for the parameters
    param_names <- c(
      colnames(X_Fixed),
      paste0("Mean(", colnames(X_rand), ")")
    )
    if(correlated) {
      chol_names <- outer(colnames(X_rand), colnames(X_rand), function(r, c) paste0("chol_", r, ":", c))
      param_names <- c(param_names, chol_names[upper.tri(chol_names, diag=TRUE)])
    } else {
      param_names <- c(param_names, paste0("SD(", colnames(X_rand), ")"))
    }
    param_names <- c(param_names, parnames)
    names(start) <- param_names
  }
    
  # --- 3. Maximum Likelihood Estimation ---
  fit <- maxLik::maxLik(
    logLik = p_nb_rp,
    start = start,
    method = method,
    control = list(iterlim = max.iters, printLevel = print.level),
    # Pass all necessary data to the likelihood function
    y = y,
    X_Fixed = X_Fixed,
    X_rand = X_rand,
    hdraws = hdraws,
    haltons = haltons,
    normed_haltons = normed_haltons,
    panel_id = data$panel_id,
    weights_vec = weights_vec,
    offset_vec = offset_vec,
    correlated = correlated,
    family = toupper(family),
    rpardists = rpardists
  )
  
  # --- 4. Post-Estimation Processing ---
  N_fixed <- ncol(X_Fixed)
  N_rand <- ncol(X_rand)
  coefs <- fit$estimate
  
  # Create descriptive names for the parameters
  names(fit$estimate) <- param_names
  
  # Extract and store key results
  fit$coefficients <- coefs[1:(N_fixed + N_rand)]
  
  if (correlated) {
    chol_vals <- coefs[grepl("chol_", names(coefs))]
    fit$Cholesky <- get_chol_cpp(chol_vals, N_rand)
    colnames(fit$Cholesky) <- rownames(fit$Cholesky) <- rand_terms
    fit$Covariance <- t(fit$Cholesky) %*% fit$Cholesky
    fit$sd <- sqrt(diag(fit$Covariance))
  } else {
    # For uncorrelated, ensure SDs are positive
    sd_vals <- coefs[grepl("SD\\(", names(coefs))]
    fit$sd <- abs(sd_vals)
    names(fit$sd) <- rand_terms
    # Update estimates to store the absolute value
    fit$estimate[grepl("SD\\(", names(coefs))] <- fit$sd
  }
  
  fit$formula <- formula
  fit$rpar_formula <- rpar_formula
  fit$numdraws <- ndraws
  fit$correlated <- correlated
  fit$family <- toupper(family)
  fit$rpardists <- rpardists
  fit$modelType <- "rpnb"
  fit$se <- sqrt(diag(solve(-fit$hessian)))
  
  # Create final object
  obj <- .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}

#' Optimized Log-Likelihood Function for Random Parameters NB Model
#'
#' @noRd
p_nb_rp <- function(p, y, X_Fixed, X_rand, hdraws, panel_id, weights_vec,
                    offset_vec, correlated, family, rpardists, haltons,
                    normed_haltons) {
  
  # --- 1. Unpack Parameters ---
  N_fixed <- ncol(X_Fixed)
  N_rand <- ncol(X_rand)
  
  fixed_coefs <- p[1:N_fixed]
  rand_means <- p[(N_fixed + 1):(N_fixed + N_rand)]
  
  # Ancillary parameters (alpha, and P for NBP)
  if (family == 'NBP' | 
      family == 'PLG' |
      family == 'PLL' |
      family == 'PGE' |
      family == 'PW' |
      family == 'SI' |
      family == 'GW') {
    ancillary_params <- utils::tail(p, N_rand * (N_rand + 1) / 2 * correlated + N_rand * (!correlated) + 2)
    alpha <- exp(ancillary_params[length(ancillary_params) - 1])
    P_val <- exp(ancillary_params[length(ancillary_params)])
    Ndistparams <- 2
  } else {
    ancillary_params <- utils::tail(p, N_rand * (N_rand + 1) / 2 * correlated + N_rand * (!correlated) + 1)
    alpha <- exp(ancillary_params[length(ancillary_params)])
    P_val <- NULL
    Ndistparams <- 1
  }
  
  # --- 2. Generate Random Parameter Draws ---
  if (correlated) {
    chol_params <- ancillary_params[1:(length(ancillary_params) - Ndistparams)]
    Ch <- get_chol_cpp(chol_params, N_rand)
    # Transform uniform Halton draws to standard normal, then to correlated normal
    norm_draws <- t(qnorm(hdraws))
    rand_pars <- t(rand_means + Ch %*% norm_draws)
  } else {
    # Uncorrelated case: use the efficient generate_draws helper
    rand_sds <- ancillary_params[1:(length(ancillary_params) - Ndistparams)]
    rand_pars <- generate_draws(hdraws, rand_means, rand_sds, rpardists)
  }
  
  # --- 3. Vectorized Likelihood Calculation ---
  # Fixed part of the mean
  mu_fixed_linear <- X_Fixed %*% fixed_coefs
  if (!is.null(offset_vec)) {
    mu_fixed_linear <- mu_fixed_linear + offset_vec
  }
  mu_fixed <- exp(mu_fixed_linear)
  
  # Random part of the mean, broadcast across all draws
  mu_rand_linear <- X_rand %*% t(rand_pars)
  
  # Combine fixed and random parts to get the full predicted mean matrix
  # Each row is an observation, each column is a draw
  pred_mat <- as.vector(mu_fixed) * exp(mu_rand_linear)
  
  
  
  # --- 4. Calculate Probabilities and Aggregate by Panel ---
  # Calculate probabilities for all observations and draws in one vectorized call
  prob_mat <- switch(family,
                     'NB2' = stats::dnbinom(y, size = alpha, mu = pred_mat),
                     'NB1' = stats::dnbinom(y, size = pred_mat / alpha, mu = pred_mat),
                     'NBP' = stats::dnbinom(y, size = (pred_mat^(2 - P_val)) / alpha, mu = pred_mat),
                     'PLN' = dpLnorm_cpp(x = y, mean = pred_mat, sigma = alpha, h = normed_haltons),
                     "PL"  = dplind_cpp(x = y, mean = pred_mat, theta = alpha),
                     "PLG" = dplindgamma_cpp(x = y, mean = pred_mat, theta = alpha, alpha = P_val, h = haltons),
                     "PLL" = dplindlogn_cpp(x = y, mean = pred_mat, theta = alpha, sigma = P_val, h = normed_haltons),
                     "PW"  = dpWeib_cpp(x = y, mean = pred_mat, alpha = alpha, sigma = P_val, h = haltons),
                     "GW"  = genWaring_cpp(x = y, mean = pred_mat, k = alpha, p = P_val),
                     "PGE" = dpge(y, mean = pred_mat, shape = alpha, scale = P_val, haltons = haltons),
                     "PIG1" = dpinvgaus(y, mu = pred_mat, eta = alpha, form = "Type 1"),
                     "PIG2" = dpinvgaus(y, mu = pred_mat, eta = alpha, form = "Type 2"),
                     "PIG" = dpinvgamma(y, mu = pred_mat, eta = alpha),
                     "SI" = dsichel(x = y, mu = pred_mat, sigma = alpha, gamma = P_val),
                     "COM" = dcom(x = y, mu = pred_mat, nu = alpha),
                     
                     # Default case for unrecognized families
                     stop(paste("Family", family, "not recognized."))
  )
  
  # Apply weights
  if (!is.null(weights_vec) && !all(weights_vec == 1)) {
    prob_mat <- prob_mat ^ weights_vec
  }
  
  # Use logs for numerical stability
  log_prob_mat <- log(prob_mat)
  log_prob_mat[prob_mat == 0] <- -1e20 # Replace -Inf with a large negative number
  
  # Sum log-probabilities by panel_id for each draw (the fastest way in R)
  panel_log_probs_sum <- rowsum(log_prob_mat, panel_id, reorder = FALSE)
  
  # --- 5. Numerically Stable Averaging of Panel Probabilities ---
  # Use the log-sum-exp trick to average probabilities without underflow/overflow
  max_log_prob <- apply(panel_log_probs_sum, 1, max)
  # Subtract max, exponentiate, get mean, log, then add max back
  log_avg_prob <- max_log_prob + log(rowMeans(exp(panel_log_probs_sum - max_log_prob)))
  
  return(log_avg_prob)
}