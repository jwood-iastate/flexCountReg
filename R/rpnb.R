#' Function for estimating a random parameter negative binomial with the ability to specify if the NB-1, NB-2, or NB-P should be used
#'
#' @name rpnb
#' @param formula an R formula. This formula should specify the outcome and the 
#'        independent variables that have fixed parameters.
#' @param rpar_formula a symbolic description of the model related specifically 
#'        to the random parameters. This should not include an outcome variable. 
#'        If the intercept is random, include it in this formula. If the 
#'        intercept is fixed, include it in \code{formula} but not in 
#'        \code{rpar_formula}. To remove the intercept, use \code{0 + vars} or 
#'        \code{-1 + vars}.
#' @param data a dataframe that has all of the variables in the \code{formula} 
#'        and \code{rpar_formula}.
#' @param form the version of the negative binomial to estimate (\code{"nb2"} 
#'        estimates the NB-2, \code{"nb1"} estimates the NB-1, \code{"nbp"} 
#'        estimates the NB-P).
#' @param rpardists an optional named vector whose names are the random 
#'        parameters and values the distribution. The distribution options 
#'        include normal ("n"), lognormal ("ln"), triangular ("t"), uniform 
#'        ("u"), and gamma ("g"). If this is not provided, normal distributions 
#'        are used for all random coefficients.
#' @param het_mean_formula an optional symbolic description of the model 
#'        (formula) for the heterogeneity in the means of the random parameters. 
#'        Variables included here act as multipliers on the mean of the random 
#'        parameters.
#' @param het_var_formula an optional symbolic description of the model 
#'        (formula) for the heterogeneity in the variances (standard deviations) 
#'        of the random parameters. Variables included here act as multipliers 
#'        on the standard deviation of the random parameters.
#' @param ndraws the number of Halton draws to use for estimating the random 
#'        parameters.
#' @param scrambled if the Halton draws should be scrambled or not. 
#'        \code{scrambled = FALSE} results in standard Halton draws while 
#'        \code{scrambled = TRUE} results in scrambled Halton draws.
#' @param correlated if the random parameters should be correlated 
#'        (\code{correlated = FALSE} results in uncorrelated random 
#'        coefficients, \code{correlated = TRUE} results in correlated random 
#'        coefficients). If the random parameters are correlated, only the 
#'        normal distribution is used for the random coefficients.
#' @param panel an optional variable or vector of variables that can be used to 
#'        specify a panel structure in the data. If this is specified, the 
#'        function will estimate the random parameters using a panel structure.
#' @param offset offset the name of a variable, or vector of variable names, in 
#'        the data frame that should be used as an offset (i.e., included but 
#'        forced to have a coefficient of 1).
#' @param method a method to use for optimization in the maximum likelihood 
#'        estimation. For options, see \code{\link[maxLik]{maxLik}}.
#' @param max.iters the maximum number of iterations to allow the optimization 
#'        method to perform.
#' @param weights the name of a variable in the data frame that should be used
#'        as a frequency weight.
#' @param start.vals an optional vector of starting values for the regression 
#'        coefficients.
#' @param verbose determines the level of verbosity for printing details of the 
#'        optimization as it is computed. A value of `FALSE` indicates no 
#'        intermediate output while `TRUE` indicates full output. Default is 
#'        `FALSE`.
#' @include helpers.R
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
#'                verbose = FALSE)
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
#'                verbose = TRUE)
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
#'                verbose = TRUE)
#'
#' summary(nbp.rp)
#' 
#' ## Random Parameters NB-2 with Heterogeneity in Means and Variances
#' # In this example, the mean of the random parameter (lnaadt) is shifted 
#' # by 'speed50', and the variance of the random parameter is scaled by 'speed50'.
#' nb2.het <- rpnb(Total_crashes ~ lnlength,
#'                 rpar_formula = ~ -1 + lnaadt,
#'                 het_mean_formula = ~ speed50,
#'                 het_var_formula = ~ speed50,
#'                 data = washington_roads,
#'                 ndraws = 50, 
#'                 form = 'nb2',
#'                 verbose = TRUE)
#'                 
#' summary(nb2.het)
#' }
#' @export
rpnb <- function(formula, rpar_formula, data, form = 'nb2',
                 rpardists = NULL,
                 het_mean_formula = NULL, 
                 het_var_formula = NULL,
                 ndraws = 1500, scrambled = FALSE,
                 correlated = FALSE, panel=NULL, 
                 weights=NULL, offset = NULL,
                 method = 'BHHH', max.iters = 1000,
                 start.vals = NULL, verbose = FALSE) {
  
  print.level = ifelse(verbose, 2, 0)
  
  # --- Data Preparation ---
  mod1_frame <- stats::model.frame(formula, data)
  y_name <- all.vars(formula)[1]
  y <- stats::model.response(mod1_frame)
  
  data <- as_tibble(data, .name_repair="minimal")
  
  # Generate panel ID
  if (!"panel_id" %in% names(data)) {
    if (is.null(panel)) {
      data <- data %>% mutate(panel_id = row_number())
    } else {
      missing_cols <- setdiff(panel, names(data))
      if (length(missing_cols) > 0) stop("Panel columns not found.")
      if (length(panel) > 1) {
        data <- data %>% unite("panel_id", all_of(panel), sep = "_", remove = FALSE)
      } else {
        data <- data %>% mutate(panel_id = as.character(data[[panel]]))
      }
    }
  }
  data <- data %>% mutate(panel_id = as.factor(panel_id))
  
  # Weights
  if (is.null(weights)){
    weights.df <- rep(1, length(y))
  }else{
    weights.df <- data %>% pull(weights)
  }
  
  if(correlated && !is.null(rpardists) && any(rpardists != "n")){
    warning("When correlated=TRUE, only normal distribution is used. Resetting rpardists to 'n'.")
    rpardists <- NULL 
  }
  
  # --- Matrix Generation ---
  X_Fixed <- modelr::model_matrix(data, formula)
  X_rand <- modelr::model_matrix(data, rpar_formula)
  
  # Heterogeneity Matrices
  if (!is.null(het_mean_formula)) {
    X_het_mean <- model.matrix(het_mean_formula, data)
    if ("(Intercept)" %in% colnames(X_het_mean)) {
      X_het_mean <- X_het_mean[, -which(colnames(X_het_mean) == "(Intercept)"), drop = FALSE]
    }
    N_het_mean <- ncol(X_het_mean)
  } else {
    X_het_mean <- NULL
    N_het_mean <- 0
  }
  
  if (!is.null(het_var_formula)) {
    X_het_var <- model.matrix(het_var_formula, data)
    if ("(Intercept)" %in% colnames(X_het_var)) {
      X_het_var <- X_het_var[, -which(colnames(X_het_var) == "(Intercept)"), drop = FALSE]
    }
    N_het_var <- ncol(X_het_var)
  } else {
    X_het_var <- NULL
    N_het_var <- 0
  }
  
  # --- Checks & Setup ---
  fix_col_names <- colnames(X_Fixed)
  rand_col_names <- colnames(X_rand)
  
  common_terms <- intersect(setdiff(fix_col_names, "(Intercept)"), 
                            setdiff(rand_col_names, "(Intercept)"))
  if(length(common_terms) > 0){
    stop("Terms cannot be in both fixed and random parts: ", paste(common_terms, collapse = ", "))
  }
  
  if("\\(Intercept\\)" %in% colnames(X_Fixed) && "\\(Intercept\\)" %in% colnames(X_rand)){
    stop("Do not include Intercept in both fixed and random formulas.")
  }
  
  if(!is.null(rpardists)){ 
    names(rpardists) <- gsub("intercept", "\\(Intercept\\)", names(rpardists))
    names(rpardists) <- gsub("constant", "\\(Intercept\\)", names(rpardists))
  }
  
  X_Fixed <- as.matrix(X_Fixed)
  X_rand <- as.matrix(X_rand)
  rpar <- colnames(X_rand)
  
  halton_draws <- function(ndraws, rpar, scrambled) {
    as.matrix(randtoolbox::halton(ndraws, length(rpar), mixed = scrambled))
  }
  hdraws <- halton_draws(ndraws, rpar, scrambled)
  
  if (!is.null(offset)){
    X_offset <- data %>% select(all_of(offset))
  } else{
    X_offset <- NULL
  }
  
  # --- Starting Values ---
  sd.start <- 0.1
  
  if(!is.null(start.vals)){
    start <- unname(start.vals)
    x_names <- names(start.vals)
  } else {
    # Initial GLM NB for fixed effects
    # BUG FIX: Use colnames() instead of names() because X_Fixed/X_rand are now matrices
    nb_vars <- c(colnames(X_Fixed), colnames(X_rand))
    nb_vars <- nb_vars[!grepl("ntercept", nb_vars)]
    nb_formula <- reformulate(nb_vars, response = y_name)
    
    nb_model <- tryCatch({
      MASS::glm.nb(nb_formula, data)
    }, error = function(e) {
      stop("Error fitting initial NB model: ", e$message)
    })
    
    params <- coef(nb_model)
    
    # 1. Fixed parameters
    start <- params[match(colnames(X_Fixed), names(params))]
    if(any(is.na(start))) start[is.na(start)] <- 0 
    x_names <- colnames(X_Fixed)
    
    # 2. Random Means
    rand_means <- params[match(rpar, names(params))]
    if(any(is.na(rand_means))) rand_means[is.na(rand_means)] <- 0
    start <- c(start, rand_means)
    x_names <- c(x_names, paste0(rpar, ':Mean'))
    
    # 3. Random Variances
    if (correlated){
      chol_starts <- numeric(0)
      for (i in 1:length(rpar)){
        for (j in 1:i){
          val <- if(i==j) 0.1 else 0
          chol_starts <- c(chol_starts, val)
          x_names <- c(x_names, paste0('Chol:', rpar[i], '_', rpar[j]))
        }
      }
      start <- c(start, chol_starts)
    } else {
      start <- c(start, rep(sd.start, length(rpar)))
      x_names <- c(x_names, paste0(rpar, ':St.Dev'))
    }
    
    # 4. Heterogeneity
    if (N_het_mean > 0) {
      start <- c(start, rep(0, N_het_mean))
      x_names <- c(x_names, paste0("HetMean:", colnames(X_het_mean)))
    }
    
    if (N_het_var > 0) {
      start <- c(start, rep(0, N_het_var))
      x_names <- c(x_names, paste0("HetVar:", colnames(X_het_var)))
    }
    
    # 5. Distribution Params
    start <- c(start, log(0.1)) 
    x_names <- c(x_names, 'ln(alpha)')
    
    if (form=='nbp'){
      start <- c(start, 1.5) 
      x_names <- c(x_names, 'P')
    }
    names(start) <- x_names
  }
  
  # --- Optimization ---
  fit <- maxLik::maxLik(p_nb_rp,
                        start = start,
                        y = y,
                        X_Fixed = X_Fixed,
                        X_rand = X_rand,
                        ndraws = ndraws,
                        hdraws = hdraws,
                        rpar = rpar,
                        correlated = correlated,
                        form = form,
                        rpardists = rpardists,
                        data = data,
                        method = method,
                        weights = weights.df,
                        offset = offset,
                        X_offset = X_offset,
                        X_het_mean = X_het_mean, 
                        X_het_var = X_het_var,   
                        control = list(iterlim = max.iters, printLevel = print.level))
  
  # --- Result Processing ---
  fit$x_names <- names(start)
  names(fit$estimate) <- names(start)
  
  coefs <- fit$estimate
  N_fixed <- ncol(X_Fixed)
  N_rand <- length(rpar)
  current_idx <- 0
  
  # 1. Fixed
  fit$coefs <- coefs[(current_idx + 1):(current_idx + N_fixed)]
  current_idx <- current_idx + N_fixed
  
  # 2. Random Means
  current_idx <- current_idx + N_rand
  
  # 3. Random Variances
  if(correlated){
    n_var <- N_rand * (N_rand + 1) / 2
    chol_vals <- coefs[(current_idx + 1):(current_idx + n_var)]
    
    Cholesky <- matrix(0, N_rand, N_rand)
    idx <- 1
    for (i in 1:N_rand) {
      for (j in 1:i) {
        Cholesky[i, j] <- chol_vals[idx]
        idx <- idx + 1
      }
    }
    fit$Cholesky <- Cholesky
    fit$Covariance <- Cholesky %*% t(Cholesky)
    fit$Correlation <- cov2cor(fit$Covariance)
    fit$sd <- sqrt(diag(fit$Covariance))
    current_idx <- current_idx + n_var
  } else {
    fit$sd <- abs(coefs[(current_idx + 1):(current_idx + N_rand)])
    current_idx <- current_idx + N_rand
  }
  
  # 4. Heterogeneity
  fit$het_mean_coefs <- if(N_het_mean > 0) coefs[(current_idx + 1):(current_idx + N_het_mean)] else NULL
  current_idx <- current_idx + N_het_mean
  
  fit$het_var_coefs <- if(N_het_var > 0) coefs[(current_idx + 1):(current_idx + N_het_var)] else NULL
  current_idx <- current_idx + N_het_var
  
  # 5. Alpha / P
  fit$alpha <- exp(coefs[current_idx + 1])
  fit$P <- if(form == 'nbp') coefs[current_idx + 2] else NULL
  
  # Meta info
  fit$formula <- formula
  fit$rpar_formula <- rpar_formula
  fit$het_mean_formula <- het_mean_formula
  fit$het_var_formula <- het_var_formula
  fit$scrambled <- scrambled
  fit$numdraws <- ndraws
  fit$correlated <- correlated
  fit$form <- form
  fit$modelType <- "rpnb"
  fit$se <- sqrt(diag(-1/(fit$hessian)))
  
  if(correlated || is.null(rpardists)) {
    fit$rpardists <- rep("n", N_rand)
    names(fit$rpardists) <- rpar
  } else {
    fit$rpardists <- rpardists
  }
  
  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}