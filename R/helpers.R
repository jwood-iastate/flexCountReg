# Helper functions for the flexCountReg package

# create function to clean data and run maxLik for bootstrapping
mod.boot <- function(data){
  # Prepare model matrices
  mod1_frame <- stats::model.frame(formula, data)
  X <- stats::model.matrix(formula, data)
  
  # If an offset is specified, create a vector for the offset
  if (!is.null(offset)){
    X_offset <- data[, offset]
  }
  
  if (!is.null(dis_param_formula_1)) {
    mod_alpha_frame <- as.matrix(modelr::model_matrix(data, dis_param_formula_1))
  } 
  
  if (!is.null(dis_param_formula_2)) {
    mod_sigma_frame <- as.matrix(modelr::model_matrix(data, dis_param_formula_2))
    
  } 
  
  if (!is.null(underreport_formula)){
    X_underreport <- as.matrix(modelr::model_matrix(data, underreport_formula))
  }
  y <- stats::model.response(mod1_frame)
  
  int_res <- maxLik::maxLik(logLikFunc, 
                            start = fit$estimate,
                            method = method, 
                            control = list(iterlim = max.iters, 
                                           printLevel = print.level))
  return(int_res)
}


# Determine model to estimate, probability distribution to use, and parameters
get_params <- function(family) {
  switch(
    family,
    "Poisson" = list(NULL, NULL),
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
  "Poisson" = function(y, predicted, alpha, sigma, haltons, normed_haltons) {
    return(stats::dpois(y, predicted))
  },
  "NB1" = function(y, predicted, alpha, sigma, haltons, normed_haltons) {
    return(stats::dnbinom(y, size = predicted/alpha, mu = predicted))
  },
  "NB2" = function(y, predicted, alpha, sigma, haltons, normed_haltons) {
    return(stats::dnbinom(y, size = alpha, mu = predicted))
  },
  "NBP" = function(y, predicted, alpha, sigma, haltons, normed_haltons){
    return(stats::dnbinom(y, size = (predicted^(2-sigma))/alpha, mu = predicted))
  },
  "PLN" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpLnorm_cpp(x=y, mean=predicted, sigma=alpha, h=normed_haltons),
  "PGE" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpge(y, mean=predicted, shape=alpha, scale=sigma, haltons=haltons),
  "PIG1" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpinvgaus(y, mu=predicted, eta=alpha),
  "PIG2" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpinvgaus(y, mu=predicted, eta=alpha, form="Type 2"),
  "PIG" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpinvgamma(y, mu=predicted, eta=alpha),
  "PL" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dplind(y, mean=predicted, theta=alpha),
  "PLG" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dplindGamma(x=y, mean=predicted, theta=alpha, alpha=sigma, hdraws=haltons),
  "PLL" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dplindLnorm(x=y, mean=predicted, theta=alpha, sigma=sigma, hdraws=normed_haltons),
  "PW" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dpWeib_cpp(y, mean=predicted, alpha=alpha, sigma=sigma, h=haltons),
  "SI" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dsichel(x=y, mu= predicted, sigma=sigma, gamma=log(alpha)),
  "GW" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dgwar(y, mu= predicted, k=alpha, rho=sigma),
  "COM" = function(y, predicted, alpha, sigma, haltons, normed_haltons) dcom(x=y, mu=predicted, nu=alpha)
  )
}


# Define the function for generating and scaling random draws using map2()
generate_draws <- function(hdraws, random_coefs_means, rand_sdevs, rpardists) {
  
  # Create a function to handle the random draws based on the distribution type
  generate_column <- function(draws_col, mean_val, sd_val, dist_type) {
    switch(dist_type,
           "n"  = stats::qnorm(draws_col, mean = mean_val, sd = abs(sd_val)),
           "ln" = stats::qlnorm(draws_col, meanlog = mean_val, sdlog = abs(sd_val)),
           "t"  = qtri(draws_col, mean_val, abs(sd_val)),  # Assuming qtri() is defined
           "u"  = mean_val + (draws_col - 0.5) * abs(sd_val),
           "g"  = {
             shape_param <- (mean_val^2) / (sd_val^2)
             rate_param <- mean_val / (sd_val^2)
             stats::qgamma(draws_col, shape = shape_param, rate = rate_param)
           },
           stats::qnorm(draws_col, mean = mean_val, sd = abs(sd_val))  # Default to normal distribution
    )
  }
  
  # Apply the generate_column function column-wise using map2
  draws <- purrr::map2(asplit(hdraws, 2), seq_along(random_coefs_means), 
                ~ generate_column(.x, random_coefs_means[.y], rand_sdevs[.y], rpardists[.y]))
  
  # Combine the list of columns back into a matrix
  draws <- do.call(cbind, draws)
  
  return(draws)
}

parse_complex_formula <- function(formula) {
  # Convert formula to string
  formula_str <- deparse(formula)
  
  # Separate the left-hand side (response variable) from the right-hand side
  rhs <- strsplit(formula_str, "~")[[1]][2]
  
  # Extract random effects terms (e.g., (x|group) or (x|group|cluster))
  random_effect_matches <- gregexpr("\\((.*?)\\)", rhs)
  random_effect_terms <- regmatches(rhs, random_effect_matches)[[1]]
  
  # Extract fixed effects (everything outside parentheses)
  fixed_effect_str <- gsub("\\((.*?)\\)", "", rhs) # Remove random effect terms
  fixed_effects <- strsplit(fixed_effect_str, "\\+|\\-")[[1]] # Split by operators
  fixed_effects <- trimws(fixed_effects) # Remove whitespace
  fixed_effects <- fixed_effects[fixed_effects != ""] # Remove empty terms
  
  # Parse random effects and build hierarchical structure
  random_effects <- list()
  for (term in random_effect_terms) {
    # Remove parentheses and split by "|"
    term_clean <- gsub("[()]", "", term)
    parts <- strsplit(term_clean, "\\|")[[1]]
    variables <- trimws(parts[1]) # Random effect variables
    levels <- trimws(parts[-1])  # Grouping levels
    random_effects[[length(random_effects) + 1]] <- list(variables = variables, levels = levels)
  }
  
  # Aggregate variables by hierarchical level
  hierarchy <- list()
  for (effect in random_effects) {
    levels_key <- paste(effect$levels, collapse = " | ") # Key by levels
    if (!levels_key %in% names(hierarchy)) {
      hierarchy[[levels_key]] <- c()
    }
    hierarchy[[levels_key]] <- unique(c(hierarchy[[levels_key]], effect$variables))
  }
  
  # Return parsed components
  list(
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    hierarchy = hierarchy
  )
}