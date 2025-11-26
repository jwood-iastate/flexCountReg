#' Predictions for flexCountReg models
#'
#' @name predict.flexCountReg
#' @param object a model object estimated using this R package.
#' @param ... optional arguments passed to the function. These include `data` 
#'        and `method`.
#' 
#' @note optional parameter `data`: a dataframe that has all of the variables in 
#'        the \code{formula} and \code{rpar_formula}.
#' @note optional parameter `method`: Only valid for random parameters models. 
#'        Options include \code{Simulated} (default) or \code{Individual}. 
#'
#' @import randtoolbox stats lamW modelr rlang 
#' @importFrom utils head tail
#' @include corr_haltons.R halton_dists.R get_chol.R helpers.R
#' 
#' #' @examples
#' \donttest{
#' # Load data and create a dummy variable
#' data("washington_roads")
#' washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
#' 
#' # =========================================================================
#' # 1. Fixed Parameter Model (Standard Countreg)
#' # =========================================================================
#' nb2_fixed <- countreg(Total_crashes ~ lnaadt + lnlength + speed50,
#'                       data = washington_roads, 
#'                       family = "NB2")
#' 
#' # Standard prediction (Expected value: mu)
#' pred_fixed <- predict(nb2_fixed, data = washington_roads)
#' head(pred_fixed)
#' 
#' # =========================================================================
#' # 2. Random Parameters Model (Countreg.rp)
#' # =========================================================================
#' # Estimate a Random Parameters NB2 model
#' rp_nb2 <- countreg.rp(Total_crashes ~ lnaadt + lnlength,
#'                       rpar_formula = ~ -1 + speed50,
#'                       data = washington_roads,
#'                       family = "NB2",
#'                       rpardists = c(speed50 = "n"),
#'                       ndraws = 100)
#' 
#' # --- Method A: Simulated Prediction (Population Average) ---
#' # Calculates E[y] = Mean over draws of exp(X*beta + random_draws)
#' # This is the standard prediction for forecasting on new data.
#' pred_sim <- predict(rp_nb2, data = washington_roads, method = "Simulated")
#' head(pred_sim)
#' 
#' # --- Method B: Individual Prediction (Bayesian Posterior) ---
#' # Calculates E[y_new | y_observed]
#' # This conditions the random parameters on the specific observed counts 
#' # for that site/observation. Useful for identifying high-risk sites.
#' # Note: The 'data' object MUST contain the outcome variable.
#' pred_ind <- predict(rp_nb2, data = washington_roads, method = "Individual")
#' head(pred_ind)
#' 
#' # =========================================================================
#' # 3. Generalized Random Parameters (Complex Structure)
#' # =========================================================================
#' # Model where dispersion (alpha) is a function of 'lnlength'
#' rp_gen <- countreg.rp(Total_crashes ~ lnaadt,
#'                       rpar_formula = ~ -1 + speed50,
#'                       dis_param_formula_1 = ~ lnlength,
#'                       data = washington_roads,
#'                       family = "NB2",
#'                       ndraws = 100)
#' 
#' # Predict takes into account the varying dispersion parameter
#' pred_gen <- predict(rp_gen, data = washington_roads, method = "Simulated")
#' head(pred_gen)
#' }
#' @export
predict.flexCountReg <- function(object, ...){
  
  # Extract optional parameters from '...'
  additional_args  <- list(...)
  
  # Handle data argument
  if (is.null(additional_args$data)) {
    data <- object$data 
  } else {
    data <- as.data.frame(additional_args$data)
  }
  
  model <- object$model
  
  if (!is.null(object$modelType)) {
    modtype <- object$modelType
  } else if (!is.null(model$modelType)) {
    modtype <- model$modelType
  } else {
    modtype <- "countreg" 
  }
  
  # Handle method argument
  if (is.null(additional_args$method)) {
    if(modtype %in% c("rppois", "countreg.rp")) {
      method <- "Simulated"
    } else {
      method <- "Exact"
    }
  } else {
    method <- additional_args$method
  }
  
  # =========================================================================
  # RANDOM PARAMETERS MODELS (RPPOIS / COUNTREG.RP)
  # =========================================================================
  if(modtype %in% c("rppois", "countreg.rp")){
    
    # --- 1. Setup & Matrix Generation ---
    form <- model$form
    family <- if(!is.null(model$family)) model$family else "NB2" 
    rpar_formula <- model$rpar_formula
    rpardists <- model$rpardists
    correlated <- model$correlated
    scrambled <- model$scrambled
    ndraws <- max(model$ndraws, 500)
    
    formula_fixed <- delete.response(terms(model$formula))
    data <- as.data.frame(data)
    
    if (!is.null(model$offset)){
      if(length(model$offset) > 1 || !model$offset %in% names(data)) {
        X_offset <- rep(0, nrow(data))
      } else {
        X_offset <- data[[model$offset]]
      }
    } else {
      X_offset <- rep(0, nrow(data))
    }
    
    # Matrices
    X_Fixed <- as.matrix(modelr::model_matrix(data, formula_fixed))
    X_rand <- as.matrix(modelr::model_matrix(data, rpar_formula))
    
    # INTERCEPT CLEANUP (Match estimation logic)
    if("(Intercept)" %in% colnames(X_rand) && !is.null(rpardists)){
      if(!any(grepl("intercept", names(rpardists), ignore.case=TRUE))){
        X_rand <- X_rand[ , colnames(X_rand) != "(Intercept)", drop=FALSE]
      }
    }
    
    # Heterogeneity Matrices
    if (!is.null(model$het_mean_formula)) {
      X_het_mean <- model.matrix(model$het_mean_formula, data)
      if ("(Intercept)" %in% colnames(X_het_mean)) X_het_mean <- X_het_mean[, -1, drop=FALSE]
    } else { X_het_mean <- NULL }
    
    if (!is.null(model$het_var_formula)) {
      X_het_var <- model.matrix(model$het_var_formula, data)
      if ("(Intercept)" %in% colnames(X_het_var)) X_het_var <- X_het_var[, -1, drop=FALSE]
    } else { X_het_var <- NULL }
    
    # Distribution Parameter Matrices
    X_dis_1 <- NULL
    X_dis_2 <- NULL
    vec_1 <- NULL
    vec_2 <- NULL
    
    # --- 2. Extract Coefficients ---
    coefs <- unlist(model$estimate, recursive = TRUE, use.names = FALSE)
    N_fixed <- ncol(X_Fixed)
    N_rand <- ncol(X_rand)
    
    current_idx <- 0
    
    # Fixed
    fixed_coefs <- coefs[(current_idx + 1):(current_idx + N_fixed)]
    current_idx <- current_idx + N_fixed
    
    # Random Means
    random_coefs_means <- coefs[(current_idx + 1):(current_idx + N_rand)]
    current_idx <- current_idx + N_rand
    
    # Random Vars
    n_var_params <- if(correlated) N_rand * (N_rand + 1) / 2 else N_rand
    rand_var_params <- coefs[(current_idx + 1):(current_idx + n_var_params)]
    current_idx <- current_idx + n_var_params
    
    # Heterogeneity
    het_mean_coefs <- NULL
    if (!is.null(X_het_mean)) {
      N_het_mean <- ncol(X_het_mean)
      het_mean_coefs <- coefs[(current_idx + 1):(current_idx + N_het_mean)]
      current_idx <- current_idx + N_het_mean
    }
    
    het_var_coefs <- NULL
    if (!is.null(X_het_var)) {
      N_het_var <- ncol(X_het_var)
      het_var_coefs <- coefs[(current_idx + 1):(current_idx + N_het_var)]
      current_idx <- current_idx + N_het_var
    }
    
    if(modtype == "countreg.rp"){
      params <- get_params(family)
      
      # Param 1
      if (!is.null(model$dis_param_formula_1)) {
        X_dis_1 <- as.matrix(modelr::model_matrix(data, model$dis_param_formula_1))
        N_dis_1 <- ncol(X_dis_1)
        p_dis_1 <- coefs[(current_idx + 1):(current_idx + N_dis_1)]
        current_idx <- current_idx + N_dis_1
        vec_1 <- exp(as.vector(X_dis_1 %*% p_dis_1))
      } else if (!is.null(params[[1]])) {
        p_dis_1 <- coefs[(current_idx + 1):(current_idx + 1)]
        current_idx <- current_idx + 1
        vec_1 <- rep(exp(p_dis_1), nrow(data))
      }
      
      # Param 2
      if (!is.null(model$dis_param_formula_2)) {
        X_dis_2 <- as.matrix(modelr::model_matrix(data, model$dis_param_formula_2))
        N_dis_2 <- ncol(X_dis_2)
        p_dis_2 <- coefs[(current_idx + 1):(current_idx + N_dis_2)]
        current_idx <- current_idx + N_dis_2
        vec_2 <- exp(as.vector(X_dis_2 %*% p_dis_2))
      } else if (!is.null(params[[2]])) {
        p_dis_2 <- coefs[(current_idx + 1):(current_idx + 1)]
        current_idx <- current_idx + 1
        vec_2 <- rep(exp(p_dis_2), nrow(data))
      }
      
    } else {
      if(!is.null(model$alpha)) vec_1 <- rep(model$alpha, nrow(data))
      if(!is.null(model$P)) vec_2 <- rep(model$P, nrow(data))
    }
    
    # --- 3. Draw Generation ---
    hdraws <- as.matrix(randtoolbox::halton(ndraws, N_rand, mixed = scrambled))
    
    draws_info <- generate_random_draws(hdraws = hdraws, 
                                        random_coefs_means = random_coefs_means, 
                                        rand_var_params = rand_var_params, 
                                        rpardists = rpardists, 
                                        rpar = colnames(X_rand), 
                                        X_rand = X_rand,
                                        het_mean_coefs = het_mean_coefs, 
                                        X_het_mean = X_het_mean,
                                        het_var_coefs = het_var_coefs, 
                                        X_het_var = X_het_var,
                                        correlated = correlated)
    
    rpar_mat <- exp(draws_info$xb_rand_mat)
    mu_fixed <- exp(as.vector(X_Fixed %*% fixed_coefs) + X_offset)
    pred_mat <- sweep(rpar_mat, 1, mu_fixed, "*")
    
    if(family == "PLN" && !is.null(vec_1)){
      adj_factor <- exp((vec_1^2)/2)
      pred_mat <- sweep(pred_mat, 1, adj_factor, "*")
    }
    
    # --- 4. Method Implementation ---
    if (method == 'Simulated') {
      return(rowMeans(pred_mat))
    }
    else if (method == 'Individual') {
      y_name <- all.vars(model$formula)[1]
      if(!y_name %in% names(data)) stop("Method 'Individual' requires outcome variable in data.")
      y_obs <- data[[y_name]]
      
      probFunc <- get_probFunc(family) 
      
      N <- nrow(pred_mat)
      D <- ncol(pred_mat)
      
      flat_y <- rep(y_obs, times = D)
      flat_pred <- as.vector(pred_mat)
      flat_alpha <- if(!is.null(vec_1)) rep(vec_1, times = D) else NULL
      flat_sigma <- if(!is.null(vec_2)) rep(vec_2, times = D) else NULL
      
      dist_haltons <- randtoolbox::halton(max(model$ndraws, 50), 1)
      normed_dist_haltons <- stats::qnorm(dist_haltons)
      
      if(family == "PLN" && !is.null(vec_1)){
        adj_factor <- exp((vec_1^2)/2)
        flat_lambda <- flat_pred / rep(adj_factor, times = D)
      } else {
        flat_lambda <- flat_pred
      }
      
      probs_flat <- probFunc(y = flat_y, 
                             predicted = flat_lambda, 
                             alpha = flat_alpha, 
                             sigma = flat_sigma,
                             haltons = dist_haltons,
                             normed_haltons = normed_dist_haltons)
      
      prob_mat <- matrix(probs_flat, nrow = N, ncol = D)
      
      numerator <- rowSums(pred_mat * prob_mat)
      denominator <- rowSums(prob_mat)
      pred_i <- ifelse(denominator == 0, rowMeans(pred_mat), numerator / denominator)
      return(pred_i)
    }
    else if (method == 'Exact') {
      stop("Exact method not supported for random parameter models. Use 'Simulated'.")
    }
  }
  
  # =========================================================================
  # FIXED PARAMETER MODELS (CountReg)
  # =========================================================================
  else {
    X <- as.matrix(modelr::model_matrix(data, model$formula))
    coefs <- unlist(model$estimate, recursive = TRUE, use.names = FALSE)
    N_x <- ncol(X)
    
    beta_pred <- as.vector(coefs[1:ncol(X)])
    
    if (!is.null(model$dis_param_formula_1)){
      alpha_X <- as.matrix(modelr::model_matrix(data, model$dis_param_formula_1))
      N_alpha <- ncol(alpha_X)
      alpha_pars <- coefs[(N_x+1):(N_x+N_alpha)]
      alpha <- exp(alpha_X %*% alpha_pars)
    } else {
      params <- get_params(model$family)
      if(!is.null(params[[1]])){
        alpha <- exp(coefs[N_x+1])
      } else {
        alpha <- NULL
      }
    }
    
    if(!is.null(model$underreport_formula)){
      X_underreport <- as.matrix(modelr::model_matrix(data, model$underreport_formula))
      N_underreport <- ncol(X_underreport)
      N_pars <- length(coefs)
      underrep_pars <- coefs[(N_pars-N_underreport+1):N_pars]
      underrep_lin <- X_underreport %*% underrep_pars
      
      if(model$underreport_family=="logit"){
        mu_adj <- 1/(1+exp(-underrep_lin))
      } else {
        mu_adj <- stats::pnorm(underrep_lin, lower.tail = FALSE)
      }
    } else {
      mu_adj <- rep(1, nrow(X))
    }
    
    pred_base <- exp(X %*% beta_pred)
    
    if (modtype=="countreg" && model$family=="PLN"){
      predictions <- pred_base * exp((alpha^2)/2) * mu_adj
    } else {
      predictions <- pred_base * mu_adj
    }
    
    return(as.vector(predictions))
  }
}