#' Function for estimating a random parameter negative binomial with the ability to specify if the NB-1, NB-2, or NB-P should be used
#'
#' @name rpnb
#' @param formula an R formula.. This formula should specify the outcome and the independent variables that have fixed parameters.
#' @param rpar_formula a symbolic description of the model related specifically to the random parameters. This should not include an outcome variable. If the intercept is random, include it in this formula. If the intercept is fixed, include it in \code{formula} but not in \code{rpar_formula}. To remove the intercept, use \code{0 + vars} or \code{-1 + vars},
#' @param data a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula},
#' @param form the version of the negative binomial to estimate (\code{"nb2"} estimates the NB-2, \code{"nb1"} estimates the NB-1, \code{"nbp"} estimates the NB-P)
#' @param rpardists an optional named vector whose names are the random parameters and values the distribution. The distribution options include normal ("n"), lognormal ("ln"), triangular ("t"), uniform ("u"), and gamma ("g"). If this is not provided, normal distributions are used for all random coefficients,
#' @param ndraws the number of Halton draws to use for estimating the random parameters,
#' @param scrambled if the Halton draws should be scrambled or not. \code{scrambled = FALSE} results in standard Halton draws while \code{scrambled = TRUE} results in scrambled Halton draws,
#' @param correlated if the random parameters should be correlated (\code{correlated = FALSE} results in uncorrelated random coefficients, \code{correlated = TRUE} results in correlated random coefficients). If the random parameters are correlated, only the normal distribution is used for the random coefficients,
#' @param panel an optional variable or vector of variables that can be used to specify a panel structure in the data. If this is specified, the function will estimate the random parameters using a panel structure,
#' @param offset offset the name of a variable, or vector of variable names, in the data frame that should be used 
#'        as an offset (i.e., included but forced to have a coefficient of 1).
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param max.iters the maximum number of iterations to allow the optimization method to perform,
#' @param weights the name of a variable in the data frame that should be used
#'        as a frequency weight.
#' @param start.vals an optional vector of starting values for the regression coefficients
#' @param print.level determines the level of verbosity for printing details of the optimization as it is computed. A value of 0 does not print out any information, a value of 1 prints minimal information, and a value of 2 prints the most information.
#' @import randtoolbox stats modelr
#' @importFrom MASS glm.nb
#' @importFrom utils head  tail
#' @importFrom dplyr mutate %>% row_number group_by across all_of summarize ungroup reframe pull
#' @importFrom tibble as_tibble
#' @importFrom tidyr unite
#' @importFrom purrr map2
#' @importFrom bbmle mle2
#' @importFrom maxLik maxLik
#' @include tri.R get_chol.R helpers.R createFlexCountReg.R
#' 
#' @examples
#' \donttest{
#'
#' ## Random Parameters Negative Binomial model (NB-1)
#' data("washington_roads")
#' nb1.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                rpar_formula = ~ speed50,
#'                data = washington_roads,
#'                ndraws = 100,
#'                correlated = FALSE,
#'                rpardists = c(intercept="u", speed50="t"),
#'                form = 'nb1',
#'                method = "nm",
#'                print.level = 2)
#'
#' summary(nb1.rp)
#'
#' ## Random Parameters Negative Binomial model (NB-2)
#' nb2.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                rpar_formula = ~ speed50,
#'                data = washington_roads,
#'                ndraws = 100,
#'                correlated = FALSE,
#'                form = 'nb2',
#'                method = "bfgs",
#'                print.level = 1)
#'
#' summary(nb2.rp)
#'
#' ## Random Parameters Negative Binomial model (NB-P)
#' nbp.rp <- rpnb(Total_crashes ~ - 1 + lnlength + lnaadt,
#'                rpar_formula = ~ speed50,
#'                data = washington_roads,
#'                ndraws = 100,
#'                correlated = FALSE,
#'                form = 'nbp',
#'                method = "bfgs",
#'                print.level = 1)
#'
#' summary(nbp.rp)}
#' @export
rpnb <- function(formula, rpar_formula, data, form = 'nb2',
                 rpardists = NULL,
                 ndraws = 1500, scrambled = FALSE,
                 correlated = FALSE, panel=NULL, 
                 weights=NULL, offset = NULL,
                 method = 'BHHH', max.iters = 1000,
                 start.vals = NULL, print.level = 0) {
  
  # Generating model matrices
  mod1_frame <- stats::model.frame(formula, data)
  y_name <- all.vars(formula)[1]
  y <- stats::model.response(mod1_frame)
  # start.vals can be a vector or a named vector with the starting values for the parameters
  # print.level is used to determine the level of details for the optimization to print (for the maxLik function call)
  
  sd.start <- 0.1 # starting value for each of the standard deviations of random parameters
  
  # ensure data is a tibble
  data <- as_tibble(data, .name_repair="minimal")
  
  # Generate a panel ID for the model - using the row number if no panel is specified
  if (!"panel_id" %in% names(data)) {
    if (is.null(panel)) {
      data <- data %>% mutate(panel_id = row_number())
    } else {
      # Check if panel columns exist
      missing_cols <- setdiff(panel, names(data))
      if (length(missing_cols) > 0) {
        stop("Panel columns not found in data: ", paste(missing_cols, collapse = ", "))
      }
      
      if (length(panel) > 1) {
        data <- data %>% unite("panel_id", all_of(panel), sep = "_", remove = FALSE)
      } else {
        data <- data %>% mutate(panel_id = as.character(data[[panel]]))
      }
    }
  }
  
  # Now safely convert panel_id to a factor
  data <- data %>% mutate(panel_id = as.factor(panel_id))
  
  # Handling Weights
  if (is.null(weights)){
    weights.df <- rep(1, length(y))
  }else{
    weights.df <- data %>% pull(weights)
  }
  
  
  # Check and correct the rpardists if the random parameters are correlated
  if(correlated && !is.null(rpardists) && any(rpardists != "n")){
    warning("When the random parameters are correlated, only the normal distribution is used.")
    rpardists <- NULL
  }
  
  # Function to generate Halton draws
  halton_draws <- function(ndraws, rpar, scrambled) {
    num_params <- length(rpar)
    halton_draws_matrix <- randtoolbox::halton(ndraws, num_params, mixed = scrambled)
    return(halton_draws_matrix)
  }
  
  # Generating model matrices
  X_Fixed <- modelr::model_matrix(data, formula)
  X_rand <- modelr::model_matrix(data, rpar_formula)
  #X_rand <-stats::model.matrix(rpar_formula, data)
  y_name <- all.vars(formula)[1]
  
  X_expanded <- cbind(X_Fixed, X_rand, X_rand) # for use in the gradient
  
  # Create named vectors
  fixed_terms <- names(X_Fixed)
  rand_terms <- names(X_rand)
  predictor_terms <- c(fixed_terms, rand_terms)
  nb_vars <- predictor_terms[!grepl("ntercept", predictor_terms)] # remove the intercept
  
  X_Fixed <- as.matrix(X_Fixed)
  X_rand <- as.matrix(X_rand)
  
  # check if the intercept is included in both the fixed and random parameters
  if("\\(Intercept\\)" %in% colnames(X_Fixed) && "\\(Intercept\\)" %in% colnames(X_rand)){
    stop("Do not include the intercept in both the fixed parameters (in `formula`) and random parameters (in `rpar_formula`). Use `- 1` in the formula you want to remove the intercept from.")
  }
  
  
  if(!is.null(rpardists)){ # check the random parameter distributions for the intercept (and correct, if needed)
    names(rpardists) <- gsub("intercept", "\\(Intercept\\)", names(rpardists)) # correct intercept name
    names(rpardists) <- gsub("constant", "\\(Intercept\\)", names(rpardists)) # correct intercept name
  }
  
  # Check if the specified distributions include all specified random parameters and no extras
  if (!correlated && !is.null(rpardists)){
    prob = FALSE
    for (i in rand_terms){
      if(i %in% names(rpardists)){
        next
      }
      else{
        print(paste("Variable ", i, " is specified as random but is not included in the `rpardists`. Please include it or do not specify the distributions."))
        prob = TRUE
      }
    }
    
    for (i in names(rpardists)){
      if (i %in% rand_terms){
        next
      }
      else{
        print(paste("Variable ", i, " is included in `rpardists` but is not in the specified random parameters. Please include it as a random parameter or remove it from `rpardists`."))
        prob = TRUE
      }
    }
    if (prob) stop() # if a problem was identified, the instructions were printed and the function will now stop executing
  }
  
  nb_formula <- reformulate(nb_vars, response = y_name)
  
  x_fixed_names <- colnames(X_Fixed)
  rpar <- colnames(X_rand)
  
  # If an offset is specified, create a vector for the offset
  if (!is.null(offset)){
    X_offset <- data %>% select(all_of(offset))
  } else{
    X_offset <- NULL
  }
  
  hdraws <- halton_draws(ndraws, rpar, scrambled)
  
  if(!is.null(start.vals)){
    params <- unname(start.vals)
    Lparams <- length(params)
    Lrpar = length(rpar)
    start <- params
    x_names <- names(start.vals)
  }
  else{
    nb_model <- tryCatch({
      MASS::glm.nb(nb_formula, data)
    }, error = function(e) {
      stop("Error fitting initial NB model for starting values: ", e$message)
    })
    params <- coef(nb_model)
    Lparams <- length(params)
    Lrpar = length(rpar)
    
    if (length(rpar)<2){
      correlated = FALSE
    }
    
    start <- params
    # Make sure the start parameters are ordered correctly, including the location of the intercept
    varnameorder <- c(x_fixed_names, rpar)
    match_indices <- match(varnameorder, names(start))
    start <- start[match_indices] # correctly ordered
    
    x_rand_names_mean <- paste0(rpar, ':Mean')
    x_names = c(x_fixed_names, x_rand_names_mean)
    
    if (correlated){
      randparam_means = tail(start, Lrpar)
      rparam_var <- rep(0.1, length(randparam_means))
      rparam_var <- diag(rparam_var)
      Chl <- chol(rparam_var)
      
      for (i in 1:length(rpar)){
        for (j in 1:length(rpar)){
          if (i >= j){
            start <- append(start, Chl[j,i])
            x_names <- append(x_names, paste('Cholesky Value for' ,paste(rpar[j], rpar[i], sep=":")))
          }
        }
      }
    }
    else{
      start <- append(start, rep(sd.start, length(rpar)))
      x_rand_names_sd <- paste0(rpar, ':St. Dev.')
      x_names <- c(x_names, x_rand_names_sd)
    }
    start <- append(start, log(0.1)) # initial log of overdispersion parameter
    x_names <- append(x_names, 'ln(alpha)')
    
    if (form=='nbp'){
      start <- append(start, 1.5) # initial value for parameter P
      x_names <- append(x_names, 'P')
    }
    names(start) <- x_names
  }
  
  # fit <- mle2(
  #   minuslogl = neg_LL,
  #   start = start,
  #   data = list(
  #     y = y,
  #     X_Fixed = X_Fixed,
  #     X_rand = X_rand,
  #     ndraws = ndraws,
  #     rpar = rpar,
  #     correlated = correlated,
  #     form = form,
  #     rpardists = rpardists,
  #     hdraws = hdraws,
  #     data = data,
  #     weights = weights.df,
  #     X_offset = X_offset,
  #     offset = offset
  #   ),
  #   control = list(maxit = max.iters, trace = print.level)
  # )
  
  fit <- maxLik::maxLik(p_nb_rp,
                        start = start,
                        y = y,
                        X_Fixed = X_Fixed,
                        X_rand = X_rand,
                        ndraws = ndraws,
                        hdraws=hdraws,
                        rpar = rpar,
                        correlated = correlated,
                        form=form,
                        rpardists=rpardists,
                        data=data,
                        method = method,
                        weights=weights.df,
                        offset=offset,
                        X_offset = X_offset,
                        control = list(iterlim = max.iters, printLevel = print.level))

  N_fixed = ncol(X_Fixed)
  N_rand = length(rpar)
  
  # names(fit$estimate) <- x_names
  
  param.splits <- as.factor(ifelse((grepl("St. Dev", x_names) + grepl("Cholesky", x_names)==1), "rpr", "coef"))
  
  coefs <- as.array(fit$estimate)
  
  # split.coefs <- split(coefs, param.splits)
  
  betas <- coefs[1:(N_fixed+N_rand)] # includes means of random parameters
  t <- coefs[(N_fixed+N_rand+1):length(coefs)]
  
  if (correlated){
    if(form=='nbp'){
      chol_vals <- coefs[(N_fixed+N_rand+1):(length(coefs)-2)]
    }
    else{
      chol_vals <-coefs[(N_fixed+N_rand+1):(length(coefs)-1)]
    }
   
    
    Cholesky <- matrix(0, N_rand, N_rand)
    chol_idx <- 1
    for (i in 1:N_rand) {
      for (j in 1:i) {  # Only lower triangle
        Cholesky[i, j] <- chol_vals[chol_idx]
        chol_idx <- chol_idx + 1
      }
    }
    
    Covariance <- t(Cholesky) %*% Cholesky
    Correlation <- cov2cor(Covariance)
  }
  else{
    if(form=='nbp'){
      sd <- abs(coefs[(N_fixed+N_rand+1):(length(coefs)-2)])
      fit$estimate[(N_fixed+N_rand+1):(length(coefs)-2)] <- sd
    }
    else{
      sd <- abs(coefs[(N_fixed+N_rand+1):(length(coefs)-1)])
      fit$estimate[(N_fixed+N_rand+1):(length(coefs)-1)] <- sd
    }
  }
  
  
  if(form =='nbp'){
    t <- tail(coefs,2)
    alpha <- exp(t[1])
    P <- t[2]
  }
  else{
    alpha <- exp(tail(fit$estimate,1))
    P <- NULL
  }
  #names(fit$estimate) <- x_names
  if(correlated){
    fit$Cholesky <- Cholesky
    fit$Covariance <- Covariance
    fit$Correlation <- Correlation
    fit$sd <- as.vector(sqrt(diag(Covariance)))
  }
  else{
    fit$sd <- abs(sd)
  }
  fit$coefs <- betas
  fit$alpha = alpha
  fit$P <- P
  #fit$estimate <- ifelse(grepl("St.", names(fit$estimate)), abs(fit$estimate ), fit$estimate)
  names(fit$estimate) <- x_names # ensure it retains the names
  fit$formula <- formula
  fit$x_names <- x_names
  fit$rpar_formula <- rpar_formula
  fit$scrambled <- scrambled
  fit$numdraws <- ndraws
  fit$correlated <- correlated
  fit$bootstraps <- NULL
  fit$se <- sqrt(diag(-1/(fit$hessian)))
  fit$form = form
  if (!correlated){
    fit$rpardists = rpardists
    names(fit$sd) <- colnames(X_rand)
  }
  else{
    dst <- as.list(rep("n", length(colnames(X_rand))))
    names(dst) <- colnames(X_rand)
    fit$rpardists = dst
  }
  
  fit$modelType <- "rpnb"
  
  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}

# function to compute probabilities
nb_prob <- function(y, mu, alpha, p = NULL, form="nb2") {
  if (form == 'nb2') {
    return(stats::dnbinom(y, size = 1/alpha, mu = mu))
  } else if (form == 'nb1') {
    return(stats::dnbinom(y, size = mu / alpha, mu = mu))
  } else if (form == 'nbp' && !is.null(p)) {
    return(stats::dnbinom(y, size = (mu^(2 - p)) / alpha, mu = mu))
  } else {
    stop("Invalid form or missing 'p' for NB-P model.")
  }
}


# Negative LL
neg_LL <- function(p, y, X_Fixed, X_rand, ndraws, rpar, correlated, form, rpardists, hdraws, data, weights, X_offset, offset){
  -sum(p_nb_rp(p, y, X_Fixed, X_rand, ndraws, rpar, correlated, form, rpardists, hdraws, data, weights, X_offset, offset))
}

# Improved main function for estimating log-likelihoods
p_nb_rp <- function(p, y, X_Fixed, X_rand, ndraws, rpar, correlated, form, 
                    rpardists, hdraws, data, weights, X_offset = NULL, offset = NULL) {
  
  # Validate inputs
  # validate_inputs(p, y, X_Fixed, X_rand, rpar, form, data, weights)
  
  # Set gradient computation method based on correlation
  exact.gradient <- !correlated
  
  # Extract coefficients
  coef_info <- extract_coefficients(p, ncol(X_Fixed), length(rpar), form)
  
  # Compute fixed effects
  mu_fixed <- compute_fixed_effects(X_Fixed, coef_info$fixed_coefs, X_offset, offset)
  
  # Generate random draws
  draws_info <- generate_random_draws(hdraws, coef_info$random_coefs_means, 
                                      coef_info$rand_sdevs, rpardists, rpar, X_rand)
  
  # Compute probabilities and log-likelihoods
  log_probs <- compute_log_likelihoods(draws_info, mu_fixed, y, coef_info$alpha, 
                                       coef_info$p, form, weights, data$panel_id)
  
  return(log_probs)
}

# Helper function to validate inputs
# validate_inputs <- function(p, y, X_Fixed, X_rand, rpar, form, data, weights) {
#   if (is.null(p) || length(p) == 0) stop("Parameter vector p cannot be empty")
#   if (is.null(y) || length(y) == 0) stop("Response vector y cannot be empty")
#   if (is.null(X_Fixed) || ncol(X_Fixed) == 0) stop("Fixed effects matrix cannot be empty")
#   if (!form %in% c("nb", "nbp")) stop("Form must be 'nb' or 'nbp'")
#   if (!"panel_id" %in% names(data)) stop("Data must contain panel_id column")
#   if (length(weights) != length(y)) stop("Weights must match response length")
# }

# Helper function to extract coefficients from parameter vector
extract_coefficients <- function(p, N_fixed, N_rand, form) {
  coefs <- as.array(p)
  
  # Extract fixed coefficients
  fixed_coefs <- head(coefs, N_fixed)
  
  # Extract random coefficients means
  random_coefs_means <- coefs[(N_fixed + 1):(N_fixed + N_rand)]
  
  # Extract distribution parameters
  if (form == 'nbp') {
    dist_params <- tail(coefs, 2)
    log_alpha <- dist_params[1]
    p_param <- dist_params[2]
    sdev_end_offset <- 2
  } else {
    dist_params <- tail(coefs, 1)
    log_alpha <- dist_params[1]
    p_param <- NULL
    sdev_end_offset <- 1
  }
  
  alpha <- exp(log_alpha)
  
  # Extract random coefficients standard deviations
  if (N_rand == 1) {
    rand_sdevs <- coefs[length(coefs) - sdev_end_offset]
  } else {
    start_idx <- N_fixed + N_rand + 1
    end_idx <- length(coefs) - sdev_end_offset
    rand_sdevs <- coefs[start_idx:end_idx]
  }
  
  return(list(
    fixed_coefs = fixed_coefs,
    random_coefs_means = random_coefs_means,
    rand_sdevs = rand_sdevs,
    alpha = alpha,
    p = p_param
  ))
}

# Helper function to compute fixed effects
compute_fixed_effects <- function(X_Fixed, fixed_coefs, X_offset, offset) {
  linear_pred <- X_Fixed %*% fixed_coefs
  
  if (!is.null(offset)) {
    if (length(offset) > 1) {
      X_offset_i <- rowsum(X_offset)
      linear_pred <- linear_pred + X_offset_i
    } else {
      linear_pred <- linear_pred + offset
    }
  }
  
  return(exp(linear_pred))
}

# Helper function to generate random draws
generate_random_draws <- function(hdraws, random_coefs_means, rand_sdevs, 
                                  rpardists, rpar, X_rand) {
  if (length(rpar) > 1) {
    draws <- generate_scaled_draws(hdraws, random_coefs_means, rand_sdevs, 
                                   rpardists, rpar)
    xb_rand_mat <- crossprod(t(X_rand), draws)
  } else {
    draws <- generate_single_param_draws(hdraws, random_coefs_means, rand_sdevs, 
                                         rpardists)
    xb_rand_mat <- sapply(draws, function(x) X_rand * x)
  }
  
  return(list(
    draws = draws,
    xb_rand_mat = xb_rand_mat
  ))
}

# Helper function for single parameter draws
generate_single_param_draws <- function(hdraws, random_coefs_means, rand_sdevs, rpardists) {
  if (is.null(rpardists)) {
    return(hdraws * rand_sdevs + random_coefs_means[1])
  }
  
  dist_type <- rpardists[1]
  mean_val <- random_coefs_means[1]
  sd_val <- abs(rand_sdevs)
  
  switch(dist_type,
         "ln" = stats::qlnorm(hdraws, mean_val, sd_val),
         "t" = qtri(hdraws, mean_val, sd_val),
         "u" = mean_val + (hdraws - 0.5) * sd_val,
         "g" = stats::qgamma(hdraws, shape = mean_val^2 / sd_val^2, 
                             rate = mean_val / sd_val^2),
         "n" = stats::qnorm(hdraws, mean_val, abs(sd_val)),
         stats::qnorm(hdraws, mean_val, abs(sd_val))  # default to normal
  )
}

# Helper function to compute log-likelihoods
compute_log_likelihoods <- function(draws_info, mu_fixed, y, alpha, p, form, 
                                    weights, panel_id) {
  # Compute random parameter matrix
  rpar_mat <- exp(draws_info$xb_rand_mat)
  
  # Compute predictions
  pred_mat <- sweep(rpar_mat, 1, mu_fixed, "*")  # More efficient than apply
  
  # Compute probabilities
  prob_mat <- apply(pred_mat, 2, nb_prob, y = y, alpha = alpha, p = p, form = form)
  
  # Apply weights
  prob_mat <- prob_mat^weights
  
  # Compute log probabilities
  log_prob_mat <- log(prob_mat)
  
  # Create data frame for grouping
  log_prob_df <- as.data.frame(log_prob_mat)
  colnames(log_prob_df) <- paste0("draw_", seq_len(ncol(log_prob_mat)))
  log_prob_df$panel_id <- panel_id
  
  # Aggregate by panel
  log_probs <- aggregate(. ~ panel_id, data = log_prob_df, FUN = sum)
  log_probs$panel_id <- NULL  # Remove panel_id column
  
  # Convert to matrix for efficiency
  log_probs_mat <- as.matrix(log_probs)
  
  # Compute final probabilities
  probs <- rowMeans(exp(log_probs_mat))
  
  return(log(probs))
}

# Improved function for generating Halton draws
generate_draws <- function(ndraws, num_params, scrambled = FALSE) {
  if (ndraws <= 0 || num_params <= 0) {
    stop("Number of draws and parameters must be positive")
  }
  
  randtoolbox::halton(ndraws, num_params, mixed = scrambled)
}

# Improved function for generating scaled draws
generate_scaled_draws <- function(hdraws, random_coefs_means, rand_sdevs, 
                                  rpardists, rpar) {
  
  if (length(rpar) == 1) {
    return(generate_single_param_draws(hdraws, random_coefs_means, rand_sdevs, rpardists))
  }
  
  # Multiple parameters case
  if (is.null(rpardists)) {
    # Default normal distribution
    return(apply(hdraws, 1, function(x) stats::qnorm(x, random_coefs_means, abs(rand_sdevs))))
  }
  
  # Custom distributions
  draws <- hdraws
  n_params <- length(rpar)
  
  for (i in seq_len(n_params)) {
    mean_val <- random_coefs_means[i]
    sd_val <- abs(rand_sdevs[i])
    
    draws[, i] <- switch(rpardists[i],
                         "ln" = stats::qlnorm(hdraws[, i], mean_val, sd_val),
                         "t" = qtri(hdraws[, i], mean_val, sd_val),
                         "u" = mean_val + (hdraws[, i] - 0.5) * sd_val,
                         "g" = stats::qgamma(hdraws[, i], shape = mean_val^2 / sd_val^2, 
                                             rate = mean_val / sd_val^2),
                         "n" = stats::qnorm(hdraws[, i], mean_val, abs(sd_val)),
                         stats::qnorm(hdraws[, i], mean_val, abs(sd_val))  # default to normal
    )
  }
  
  return(t(draws))
}

# Additional utility function for parameter extraction validation
validate_parameter_vector <- function(p, expected_length) {
  if (length(p) != expected_length) {
    stop(paste("Parameter vector length mismatch. Expected:", expected_length, 
               "Got:", length(p)))
  }
}

# Function to create coefficient summary
summarize_coefficients <- function(coef_info) {
  list(
    n_fixed = length(coef_info$fixed_coefs),
    n_random = length(coef_info$random_coefs_means),
    alpha = coef_info$alpha,
    has_p_param = !is.null(coef_info$p)
  )
}