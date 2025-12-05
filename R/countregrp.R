#' Random Parameters Count Regression Models
#'
#' @name countreg.rp
#' @param formula an R formula. This formula should specify the outcome and the
#'   independent variables that have fixed parameters.
#' @param rpar_formula a symbolic description of the model related specifically
#'   to the random parameters. This should not include an outcome variable. If
#'   the intercept is random, include it in this formula. If the intercept is
#'   fixed, include it in \code{formula} but not in \code{rpar_formula}.
#' @param data a dataframe that has all of the variables in the \code{formula}
#'   and \code{rpar_formula}.
#' @param family the name of the distribution/model type to estimate. Default is
#'   "NB2". Options include "Poisson", "NB1", "NB2", "NBP", "PIG", "Sichel",
#'   etc. (See \code{\link{countreg}} for full list).
#' @param rpardists an optional named vector whose names are the random
#'   parameters and values the distribution. The distribution options include
#'   normal ("n"), lognormal ("ln"), triangular ("t"), uniform ("u"), and gamma
#'   ("g"). If not provided, normal is used.
#' @param dis_param_formula_1 a symbolic description of the model for the 
#'        first parameter of the count distribution (e.g., ln(alpha) for NB2).
#' @param dis_param_formula_2 a symbolic description of the model for the 
#'        second parameter of the count distribution (if applicable).
#' @param het_mean_formula an optional symbolic description of the model 
#'        for heterogeneity in the means of the random parameters.
#' @param het_var_formula an optional symbolic description of the model 
#'        for heterogeneity in the variances of the random parameters.
#' @param ndraws the number of Halton draws to use for estimating the random 
#'        parameters.
#' @param scrambled if the Halton draws should be scrambled.
#' @param correlated if the random parameters should be correlated. 
#'        If TRUE, only normal distributions are used.
#' @param panel_id an optional variable name (string) or vector defining the
#'   panel structure (repeated measures). If provided, the standard errors and
#'   likelihood are estimated accounting for the panel structure.
#' @param offset variable name to be used as an offset.
#' @param weights variable name to be used as frequency weights.
#' @param method optimization method (e.g., "BHHH", "BFGS", "NM").
#' @param max.iters maximum number of iterations.
#' @param start.vals optional vector of starting values.
#' @param verbose logical.
#' 
#' @examples
#' \donttest{
#' # Load data
#' data("washington_roads")
#' washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
#' 
#' # 1. Basic Random Parameters Negative Binomial (NB2)
#' rp_nb2 <- countreg.rp(Total_crashes ~ lnaadt + lnlength,
#'                       rpar_formula = ~ -1 + speed50,
#'                       data = washington_roads,
#'                       family = "NB2",
#'                       rpardists = c(speed50 = "n"),
#'                       ndraws = 100,
#'                       method = "BHHH")
#' summary(rp_nb2)
#' 
#' # 2. Random Parameters with Panel Structure (if 'site_id' exists)
#' # rp_panel <- countreg.rp(Total_crashes ~ -1 + lnaadt,
#' #                        rpar_formula = ~ speed50,
#' #                        data = washington_roads,
#' #                        panel_id = "site_id",
#' #                        family = "NB2",
#' #                        ndraws = 100)
#' 
#' # 3. Generalized Random Parameters Model with Heterogeneity
#' rp_gen <- countreg.rp(Total_crashes ~ lnaadt,
#'                       rpar_formula = ~ -1 + speed50,
#'                       dis_param_formula_1 = ~ lnlength, 
#'                       het_mean_formula = ~ AADT10kplus,
#'                       data = washington_roads,
#'                       family = "NB2",
#'                       rpardists = c(speed50 = "n"),
#'                       ndraws = 100)
#' summary(rp_gen)
#' 
#' # 4. Random Parameters Poisson Model with panel specification
#' rp_poisson <- countreg.rp(Total_crashes ~ lnaadt,
#'                       rpar_formula = ~ -1 + speed50,
#'                       dis_param_formula_1 = ~ lnlength, 
#'                       het_mean_formula = ~ AADT10kplus,
#'                       data = washington_roads,
#'                       family = "POISSON",
#'                       rpardists = c(speed50 = "n"),
#'                       ndraws = 100,
#'                       panel_id = "ID")
#' summary(rp_poisson)
#' }
#' 
#' @importFrom maxLik maxLik
#' @importFrom randtoolbox halton
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom stats terms reformulate coef
#' @importFrom dplyr mutate pull select all_of
#' @importFrom tibble as_tibble
#' @importFrom tidyr unite
#' @include helpers.R tri.R  createFlexCountReg.R
#' @export
countreg.rp <- function(formula, rpar_formula, data, family = "NB2",
                         rpardists = NULL,
                         dis_param_formula_1 = NULL,
                         dis_param_formula_2 = NULL,
                         het_mean_formula = NULL, 
                         het_var_formula = NULL,
                         ndraws = 500, scrambled = FALSE,
                         correlated = FALSE, panel_id = NULL, 
                         weights=NULL, offset = NULL,
                         method = 'BHHH', max.iters = 1000,
                         start.vals = NULL, verbose = FALSE) {
  
  print.level <- ifelse(verbose, 2, 0)
  
  # --- 1. Family and Parameter Setup ---
  family <- toupper(gsub("[^[:alnum:]]", "", family))
  params <- get_params(family) # From helpers.R
  # NOTE: probFunc is now retrieved inside the likelihood function to prevent
  # scoping issues.
  
  # --- 2. Data Preparation ---
  mod1_frame <- stats::model.frame(formula, data)
  y_name <- all.vars(formula)[1]
  y <- stats::model.response(mod1_frame)
  
  data <- as_tibble(data, .name_repair="minimal")
  
  # Panel ID Generation
  if (is.null(panel_id)) {
    nn <- nrow(data)
    data$panel_id_internal <- factor(1:nn)
  } else {
    cond1 <- (length(panel_id) == 1) && (is.character(panel_id)) && 
      (panel_id %in% names(data))
    if (cond1) {
      data$panel_id_internal <- factor(data[[panel_id]])
    } else if (length(panel_id) == nrow(data)) {
      data$panel_id_internal <- factor(panel_id)
    } else {
      msg <- paste0("panel_id must be a column name in 'data' or a vector of ", 
                    "the same length as the data.")
      warning(msg)
    }
  }
  
  panel_group <- data$panel_id_internal 
  
  # Weights
  if (is.null(weights)){
    weights.df <- rep(1, length(y))
  }else{
    weights.df <- data %>% pull(weights)
  }
  
  # Offset
  if (!is.null(offset)){
    
    X_offset <- data %>% select(all_of(offset))
    
    if (ncol(X_offset) > 1){
      X_offset <- rowSums(X_offset) 
    } else {
      X_offset <- X_offset[[1]]
      }
  } else {
    
    X_offset <- rep(0, nrow(data))
  }
  
  # Handle Correlated Flag
  if(correlated && !is.null(rpardists) && any(rpardists != "n")){
    msg <- paste0("When correlated=TRUE, only normal distribution is used. ", 
                  "Resetting rpardists to 'n'.")
    warning(msg)
    rpardists <- NULL 
  }
  
  # --- 3. Matrix Generation ---
  
  # Fixed Matrix
  X_Fixed <- as.matrix(modelr::model_matrix(data, formula))
  
  # Random Matrix
  X_rand <- as.matrix(modelr::model_matrix(data, rpar_formula))
  
  # FIX: Auto-remove Intercept if present but not in rpardists
  # This prevents mismatch errors when formula is ~var but rpardists =
  # c(var="n")
  if("(Intercept)" %in% colnames(X_rand) && !is.null(rpardists)){
    # Check if user explicitly named "Intercept" or "(Intercept)" in rpardists
    if(!any(grepl("intercept", names(rpardists), ignore.case = TRUE))){
      # If not, remove the intercept column
      X_rand <- X_rand[ , colnames(X_rand) != "(Intercept)", drop = FALSE]
    }
  }
  
  # Distribution Parameter Matrices
  if (is.null(params[[1]])) {
    N_dis_1 <- 0; X_dis_1 <- NULL
  } else if (!is.null(dis_param_formula_1)) {
    X_dis_1 <- as.matrix(modelr::model_matrix(data, dis_param_formula_1))
    N_dis_1 <- ncol(X_dis_1)
  } else {
    X_dis_1 <- NULL; N_dis_1 <- 1 
  }
  
  if (is.null(params[[2]])) {
    N_dis_2 <- 0; X_dis_2 <- NULL
  } else if (!is.null(dis_param_formula_2)) {
    X_dis_2 <- as.matrix(modelr::model_matrix(data, dis_param_formula_2))
    N_dis_2 <- ncol(X_dis_2)
  } else {
    X_dis_2 <- NULL; N_dis_2 <- 1 
  }
  
  # Heterogeneity Matrices
  if (!is.null(het_mean_formula)) {
    X_het_mean <- model.matrix(het_mean_formula, data)
    if ("(Intercept)" %in% colnames(X_het_mean)) {
      X_het_mean <- X_het_mean[, -which(colnames(X_het_mean) == "(Intercept)"), 
                               drop = FALSE]
    }
    N_het_mean <- ncol(X_het_mean)
  } else {
    X_het_mean <- NULL; N_het_mean <- 0
  }
  
  if (!is.null(het_var_formula)) {
    X_het_var <- model.matrix(het_var_formula, data)
    if ("(Intercept)" %in% colnames(X_het_var)) {
      X_het_var <- X_het_var[, -which(colnames(X_het_var) == "(Intercept)"), 
                             drop = FALSE]
    }
    N_het_var <- ncol(X_het_var)
  } else {
    X_het_var <- NULL; N_het_var <- 0
  }
  
  cond <- 
    "\\(Intercept\\)" %in% colnames(X_Fixed) && 
    "\\(Intercept\\)" %in% colnames(X_rand) 
  if (cond){
    warning("Do not include Intercept in both fixed and random formulas.")
  }
  
  rpar <- colnames(X_rand)
  N_fixed <- ncol(X_Fixed)
  N_rand <- length(rpar)
  
  # Draws
  hdraws <- 
    as.matrix(randtoolbox::halton(ndraws, length(rpar), mixed = scrambled))
  
  # Haltons for the Distribution
  dist_haltons <- randtoolbox::halton(ndraws, 1) 
  normed_dist_haltons <- stats::qnorm(dist_haltons)
  
  # --- 4. Starting Values ---
  if(!is.null(start.vals)){
    start <- unname(start.vals)
    x_names <- names(start.vals)
  } else {
    # 1. Fixed & Random Means (via GLM Poisson)
    all_vars <- c(colnames(X_Fixed), colnames(X_rand))
    all_vars <- all_vars[!grepl("ntercept", all_vars)] 
    init_form <- reformulate(all_vars, response = y_name)
    
    glm_mod <- tryCatch(glm(init_form, data = data, family = poisson), 
                        error = function(e) NULL)
    
    if(is.null(glm_mod)) {
      start <- rep(0, N_fixed + N_rand)
    } else {
      glm_coef <- coef(glm_mod)
      s_fixed <- glm_coef[match(colnames(X_Fixed), names(glm_coef))]
      s_fixed[is.na(s_fixed)] <- 0
      
      s_rand <- glm_coef[match(colnames(X_rand), names(glm_coef))]
      s_rand[is.na(s_rand)] <- 0
      
      start <- c(s_fixed, s_rand)
    }
    x_names <- c(colnames(X_Fixed), paste0(rpar, ":Mean"))
    
    # 2. Random Variances
    if (correlated){
      chol_starts <- numeric(0)
      chol_names <- numeric(0)
      for (i in seq_along(rpar)){
        for (j in 1:i){
          val <- if(i==j) 0.1 else 0
          chol_starts <- c(chol_starts, val)
          chol_names <- c(chol_names, paste0('Chol:', rpar[i], '_', rpar[j]))
        }
      }
      start <- c(start, chol_starts)
      x_names <- c(x_names, chol_names)
    } else {
      start <- c(start, rep(0.1, N_rand))
      x_names <- c(x_names, paste0(rpar, ':St.Dev'))
    }
    
    # 3. Heterogeneity
    if(N_het_mean > 0) {
      start <- c(start, rep(0, N_het_mean))
      x_names <- c(x_names, paste0("HetMean:", colnames(X_het_mean)))
    }
    if(N_het_var > 0) {
      start <- c(start, rep(0, N_het_var))
      x_names <- c(x_names, paste0("HetVar:", colnames(X_het_var)))
    }
    
    # 4. Distribution Parameters
    if(N_dis_1 > 0) {
      if(is.null(X_dis_1)) {
        start <- c(start, -2) 
        x_names <- c(x_names, params[[1]])
      } else {
        start <- c(start, rep(0, N_dis_1)) 
        x_names <- c(x_names, paste0(params[[1]], ":", colnames(X_dis_1)))
      }
    }
    
    if(N_dis_2 > 0) {
      if(is.null(X_dis_2)) {
        start <- c(start, -2) 
        x_names <- c(x_names, params[[2]])
      } else {
        start <- c(start, rep(0, N_dis_2)) 
        x_names <- c(x_names, paste0(params[[2]], ":", colnames(X_dis_2)))
      }
    }
    
    names(start) <- x_names
  }
  
  # --- 5. Likelihood Function ---
  
  ll_countreg_rp <- function(p) {
    # DEFINE probFunc INSIDE THE LIKELIHOOD TO FIX SCOPING
    local_probFunc <- get_probFunc(family)
    
    if(is.null(local_probFunc)) {
      warning(paste0("Probability function not found for family: ", family))
    }
    
    current_idx <- 0
    
    fixed_coefs <- p[(current_idx + 1):(current_idx + N_fixed)]
    current_idx <- current_idx + N_fixed
    
    rand_means <- p[(current_idx + 1):(current_idx + N_rand)]
    current_idx <- current_idx + N_rand
    
    n_var_params <- if(correlated) N_rand * (N_rand + 1) / 2 else N_rand
    rand_vars <- p[(current_idx + 1):(current_idx + n_var_params)]
    current_idx <- current_idx + n_var_params
    
    het_mean_coefs <- NULL
    if (N_het_mean > 0) {
      het_mean_coefs <- p[(current_idx + 1):(current_idx + N_het_mean)]
      current_idx <- current_idx + N_het_mean
    }
    
    het_var_coefs <- NULL
    if (N_het_var > 0) {
      het_var_coefs <- p[(current_idx + 1):(current_idx + N_het_var)]
      current_idx <- current_idx + N_het_var
    }
    
    p_dis_1 <- NULL
    if (N_dis_1 > 0) {
      p_dis_1 <- p[(current_idx + 1):(current_idx + N_dis_1)]
      current_idx <- current_idx + N_dis_1
    }
    
    p_dis_2 <- NULL
    if (N_dis_2 > 0) {
      p_dis_2 <- p[(current_idx + 1):(current_idx + N_dis_2)]
      current_idx <- current_idx + N_dis_2
    }
    
    # Construct Distribution Vectors
    if (N_dis_1 > 0) {
      if (is.null(X_dis_1)) {
        val_1 <- exp(p_dis_1) 
        vec_1 <- rep(val_1, nrow(data))
      } else {
        vec_1 <- exp(as.vector(X_dis_1 %*% p_dis_1))
      }
    } else {
      vec_1 <- NULL
    }
    
    if (N_dis_2 > 0) {
      if (is.null(X_dis_2)) {
        val_2 <- exp(p_dis_2)
        vec_2 <- rep(val_2, nrow(data))
      } else {
        vec_2 <- exp(as.vector(X_dis_2 %*% p_dis_2))
      }
    } else {
      vec_2 <- NULL
    }
    
    # Calculate Mu (Obs x Draws)
    mu_fixed <- exp(as.vector(X_Fixed %*% fixed_coefs) + X_offset)
    
    draws_info <- generate_random_draws(hdraws = hdraws, 
                                        random_coefs_means = rand_means, 
                                        rand_var_params = rand_vars, 
                                        rpardists = rpardists, 
                                        rpar = rpar, 
                                        X_rand = X_rand,
                                        het_mean_coefs = het_mean_coefs, 
                                        X_het_mean = X_het_mean,
                                        het_var_coefs = het_var_coefs, 
                                        X_het_var = X_het_var,
                                        correlated = correlated)
    
    rpar_mult <- exp(draws_info$xb_rand_mat)
    mu_mat <- sweep(rpar_mult, 1, mu_fixed, "*") 
    
    # Calculate Probabilities
    N <- nrow(mu_mat)
    D <- ncol(mu_mat)
    
    flat_mu <- as.vector(mu_mat) 
    flat_y <- rep(y, times = D)
    
    flat_alpha <- if(!is.null(vec_1)) rep(vec_1, times = D) else NULL
    flat_sigma <- if(!is.null(vec_2)) rep(vec_2, times = D) else NULL
    
    # Use the local_probFunc variable
    probs_flat <- local_probFunc(y = flat_y, 
                                 predicted = flat_mu, 
                                 alpha = flat_alpha, 
                                 sigma = flat_sigma,
                                 haltons = dist_haltons, 
                                 normed_haltons = normed_dist_haltons)
    
    probs_mat <- matrix(probs_flat, nrow = N, ncol = D)
    probs_mat <- probs_mat^weights.df 
    
    # Panel Aggregation
    log_probs_mat <- log(probs_mat)
    sum_log_probs <- rowsum(log_probs_mat, group = panel_group, reorder = FALSE)
    lik_panel_draws <- exp(sum_log_probs)
    lik_panel <- rowMeans(lik_panel_draws)
    
    ll_total <- log(lik_panel)
    
    return(ll_total)
  }
  
  # --- 6. Optimization ---
  fit <- maxLik::maxLik(
    ll_countreg_rp,
    start = start,
    method = method,
    control = list(iterlim = max.iters, printLevel = print.level))
  
  # --- 7. Result Packaging ---
  fit$x_names <- names(start)
  names(fit$estimate) <- names(start)
  
  fit$formula <- formula
  fit$rpar_formula <- rpar_formula
  fit$het_mean_formula <- het_mean_formula
  fit$het_var_formula <- het_var_formula
  fit$dis_param_formula_1 <- dis_param_formula_1
  fit$dis_param_formula_2 <- dis_param_formula_2
  fit$family <- family
  fit$modelType <- "countreg.rp"
  fit$ndraws <- ndraws
  fit$scrambled <- scrambled
  fit$correlated <- correlated
  fit$panel_id <- panel_id 
  fit$rpardists <- if(is.null(rpardists) || correlated) {
    structure(rep("n", length(rpar)), names=rpar) 
  } else {
    rpardists
  }
  
  fit$se <- sqrt(diag(-1/(fit$hessian)))
  
  obj <- .createFlexCountReg(model = fit, 
                             data = data, 
                             call = match.call(), 
                             formula = formula)
  return(obj)
}