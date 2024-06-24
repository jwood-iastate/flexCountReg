#' Function for estimating a random parameters Poisson regression model
#'
#' This function estimates a Poisson regression model with random parameters.
#'
#' @name rppois
#' @param formula An R formula specifying the outcome and the independent variables with fixed parameters.
#' @param rpar_formula A symbolic description of the model related to the random parameters. This should not include an outcome variable. If the intercept is random, include it in this formula. If the intercept is fixed, include it in \code{formula} but not in \code{rpar_formula}. To remove the intercept, use \code{0 + vars} or \code{-1 + vars}.
#' @param data A dataframe containing all of the variables in the \code{formula} and \code{rpar_formula}.
#' @param rpardists An optional named vector specifying the distributions of the random parameters. Options include normal ("n"), lognormal ("ln"), triangular ("t"), uniform ("u"), and gamma ("g"). If not provided, normal distributions are used for all random coefficients.
#' @param ndraws The number of Halton draws to use for estimating the random parameters.
#' @param scrambled Logical, whether the Halton draws should be scrambled. \code{FALSE} results in standard Halton draws, while \code{TRUE} results in scrambled Halton draws.
#' @param correlated Logical, whether the random parameters should be correlated. \code{FALSE} results in uncorrelated random coefficients, \code{TRUE} results in correlated random coefficients. If \code{TRUE}, only the normal distribution is used for the random coefficients.
#' @param method The optimization method to use for maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}}.
#' @param max.iters The maximum number of iterations for the optimization method.
#' @param start.vals An optional vector of starting values for the regression coefficients.
#' @param print.level Determines the verbosity level for printing optimization details. \code{0} prints no information, \code{1} prints minimal information, \code{2} prints detailed information.
#' @import MASS nlme randtoolbox maxLik stats modelr
#' @importFrom utils head tail
#' @include tri.R createFlexCountReg.R
#' @examples
#' \donttest{
#' ## Random Parameters Poisson Regression model
#' data("washington_roads")
#' rpp <- rppois(Total_crashes ~ -1 + lnlength + lnaadt,
#'               rpar_formula = ~ speed50,
#'               data = washington_roads,
#'               ndraws = 100,
#'               correlated = FALSE,
#'               rpardists = c(intercept = "u", speed50 = "t"),
#'               method = "bfgs",
#'               print.level = 2)
#' summary(rpp)
#' 
#' Generate predictions and create a CURE Plot
#' }
#' @export
rppois <- function(formula, rpar_formula, data,
                   rpardists = NULL, ndraws = 1500, scrambled = FALSE,
                   correlated = FALSE, method = 'BHHH', max.iters = 1000,
                   start.vals = NULL, print.level = 0) {
  
  mod1_frame <- stats::model.frame(formula, data)
  X_Fixed <- modelr::model_matrix(data, formula)
  X_rand <- modelr::model_matrix(data, rpar_formula)
  y <- stats::model.response(mod1_frame)
  
  X_Fixed <- as.matrix(X_Fixed)
  X_rand <- as.matrix(X_rand)
  rpar <- colnames(X_rand)
  x_fixed_names <- colnames(X_Fixed)
  
  if ("\\(Intercept\\)" %in% colnames(X_Fixed) && "\\(Intercept\\)" %in% rpar) {
    stop("Do not include the intercept in both the fixed parameters (in `formula`) and random parameters (in `rpar_formula`). Use `-1` in the formula you want to remove the intercept from.")
  }
  
  if (!is.null(rpardists)) { 
    names(rpardists) <- gsub("intercept|constant|Constant", "\\(Intercept\\)", names(rpardists))
  }
  
  if (!correlated && !is.null(rpardists)) {
    rand_term_check <- setdiff(names(rpardists), rpar)
    if (length(rand_term_check) > 0) {
      stop(paste("Variables ", paste(rand_term_check, collapse = ", "), " are included in `rpardists` but not in the specified random parameters."))
    }
  }
  
  y_name <- all.vars(formula)[1]
  fixed_terms <- colnames(X_Fixed)
  rand_terms <- colnames(X_rand)
  predictor_terms <- c(fixed_terms, rand_terms)
  nb_vars <- predictor_terms[!grepl("ntercept", predictor_terms)] # remove the intercept
  
  nb_formula <- reformulate(nb_vars, response = y_name)
  
  x_fixed_names <- colnames(X_Fixed)
  rpar <- colnames(X_rand)
  
  hdraws <- halton_draws(ndraws, rpar, scrambled)
  
  if (length(rpar)<2){
    correlated = FALSE
  }
  
  if (!is.null(start.vals)) {
    start <- unname(start.vals)
    x_names <- names(start.vals)
  } else {
    x_rand_names_mean <- paste0(rpar, ':Mean')
    x_names = c(x_fixed_names, x_rand_names_mean)
    
    nb_model <- MASS::glm.nb(nb_formula, data)
    params <- coef(nb_model)
    start <- params
    
    # Make sure the start parameters are ordered correctly, including the location of the intercept
    varnameorder <- c(x_fixed_names, rpar)
    match_indices <- match(varnameorder, names(start))
    start <- start[match_indices] # correctly ordered
    
    hdraws <- halton_draws(ndraws, length(rpar), scrambled)
    
    if (correlated && length(rpar) > 1) {
      randparam_means = tail(start, length(rpar))
      rparam_var <- abs(randparam_means)/2
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
    else {
      x_rand_names_sd <- paste0(rpar, ':St. Dev.')
      start <- c(start, rep(0.1, length(rpar)))
      if(length(rpar) > 1){
        x_names <- c(x_names, x_rand_names_sd)
      }
      else{
        x_names <- append(x_names, x_rand_names_sd)
      }
    }
  }
  names(start) <- x_names
  
  p_pois_rp <- function(p, y, X_Fixed, X_rand, ndraws, hdraws, correlated, rpardists) {
    N_fixed <- ncol(X_Fixed)
    N_rand <- ncol(X_rand)
    
    fixed_coefs <- p[1:N_fixed]
    rand_means <- p[(N_fixed + 1):(N_fixed + N_rand)]
    rand_sdevs <- if (correlated) NULL else p[(N_fixed + N_rand + 1):(N_fixed + 2 * N_rand)]
    
    mu_fixed <- exp(X_Fixed %*% fixed_coefs)
    
    if (correlated && N_rand > 1) {
      chol_vals <- p[(N_fixed + N_rand + 1):(length(p))]
      Ch <- matrix(0, N_rand, N_rand)
      Ch[lower.tri(Ch, diag = TRUE)] <- chol_vals
      scaled_draws <- qnorm(hdraws) %*% Ch
    } else {
      if (is.null(rpardists)) {
        scaled_draws <- hdraws * rep(rand_sdevs, each = ndraws)
      } else {
        scaled_draws <- matrix(NA, nrow = ndraws, ncol = N_rand)
        for (i in 1:N_rand) {
          dist_type <- rpardists[colnames(X_rand)[i]]
          scaled_draws[, i] <- switch(dist_type,
                                      ln = stats::qlnorm(hdraws[, i], rand_means[i], abs(rand_sdevs[i])),
                                      t = qtri(hdraws[, i], rand_means[i], abs(rand_sdevs[i])),
                                      u = rand_means[i] + (hdraws[, i] - 0.5) * abs(rand_sdevs[i]),
                                      g = stats::qgamma(hdraws[, i], shape = rand_means[i]^2 / (rand_sdevs[i]^2), rate = rand_means[i] / (rand_sdevs[i]^2)),
                                      stats::qnorm(hdraws[, i], rand_means[i], abs(rand_sdevs[i])))
        }
      }
    }
    
    xb_rand_mat <- if (N_rand > 1) {
      crossprod(t(X_rand), t(scaled_draws))
    } else {
      apply(scaled_draws, function(x) X_rand * x)
    }
    
    rpar_mat <- exp(xb_rand_mat)
    pred_mat <- apply(rpar_mat, 2, function(x) x * mu_fixed)

    prob_mat <- apply(pred_mat, 2, pois_prob, y = y)
    
    log_probs <- log(rowMeans(prob_mat))
    if (method == 'BHHH') log_probs else sum(log_probs)
  }
  
  gradFun <- function(params, y, X_Fixed, X_rand, ndraws, hdraws, correlated, rpardists) {
    N_fixed <- ncol(X_Fixed)
    N_rand <- ncol(X_rand)
    
    fixed_coefs <- params[1:N_fixed]
    rand_means <- params[(N_fixed + 1):(N_fixed + N_rand)]
    
    if (!correlated) {
      rand_sdevs <- params[(N_fixed + N_rand + 1):(N_fixed + 2 * N_rand)]
    }
    
    mu_fixed <- exp(X_Fixed %*% fixed_coefs)
    
    if (correlated && N_rand > 1) {
      chol_vals <- params[(N_fixed + N_rand + 1):length(params)]
      Ch <- matrix(0, N_rand, N_rand)
      Ch[lower.tri(Ch, diag = TRUE)] <- chol_vals
      scaled_draws <- qnorm(hdraws) %*% Ch
    } else {
      if (is.null(rpardists)) {
        scaled_draws <- hdraws * matrix(rep(rand_sdevs, each = ndraws), nrow = ndraws)
      } else {
        scaled_draws <- matrix(NA, nrow = ndraws, ncol = N_rand)
        for (i in 1:N_rand) {
          dist_type <- rpardists[colnames(X_rand)[i]]
          scaled_draws[, i] <- switch(dist_type,
                                      ln = stats::qlnorm(hdraws[, i], rand_means[i], abs(rand_sdevs[i])),
                                      t = qtri(hdraws[, i], rand_means[i], abs(rand_sdevs[i])),
                                      u = rand_means[i] + (hdraws[, i] - 0.5) * abs(rand_sdevs[i]),
                                      g = stats::qgamma(hdraws[, i], shape = rand_means[i]^2 / (rand_sdevs[i]^2), rate = rand_means[i] / (rand_sdevs[i]^2)),
                                      stats::qnorm(hdraws[, i], rand_means[i], abs(rand_sdevs[i])))
        }
      }
    }
    
    xb_rand_mat <- if (N_rand > 1) {
      crossprod(t(X_rand), t(scaled_draws))
    } else {
      apply(scaled_draws, function(x) X_rand * x)
    }
    
    rpar_mat <- exp(xb_rand_mat)
    pred_mat <- apply(rpar_mat, 2, function(x) x * mu_fixed)
    mu <- rowMeans(pred_mat)
    resid <- y - mu
    
    grad <- matrix(0, nrow = length(y), ncol = length(params))
    
    pois_prob <- function(y, mu) {
      dpois(y, mu)
    }
    
    # Gradient for fixed coefficients
    grad[, 1:N_fixed] <- sweep(X_Fixed, 1, resid, `*`)
    
    # Gradient for random means
    grad[, (N_fixed + 1):(N_fixed + N_rand)] <- sweep(X_rand, 1, resid, `*`)
    
    # Gradient for random standard deviations if not correlated
    if (!correlated) {
      if (N_rand == 1) {
        grad[, (N_fixed + N_rand + 1)] <- scaled_draws * resid
      } else {
        for (j in 1:N_rand) {
          grad[, (N_fixed + N_rand + j)] <- scaled_draws[, j] * resid
        }
      }
    } else {
      # Gradient for Cholesky values
      counter <- 1
      for (k in 1:N_rand) {
        for (j in 1:k) {
          grad[, (N_fixed + N_rand + counter)] <- rowSums(sweep(scaled_draws, 2, resid, `*`) * sweep(scaled_draws, 2, Ch[j, ], `*`))
          counter <- counter + 1
        }
      }
    }
    
    return(grad)
  }
  
  fit <- maxLik::maxLik(
    logLik = p_pois_rp,
    start = start,
    grad = if (method == 'BHHH') gradFun else function(...) colSums(gradFun(...)),
    hess = NULL,
    y = y,
    X_Fixed = X_Fixed,
    X_rand = X_rand,
    ndraws = ndraws,
    hdraws = hdraws,
    correlated = correlated,
    rpardists = rpardists,
    method = method,
    control = list(iterlim = max.iters, printLevel = print.level)
  )
  
  N_fixed <- ncol(X_Fixed)
  N_rand <- ncol(X_rand)
  
  param.splits <- as.factor(ifelse((grepl("St. Dev", x_names) + grepl("Cholesky", x_names)==1), "rpr", "coef"))
  
  coefs <- as.array(fit$estimate)
  
  split.coefs <- split(coefs, param.splits)
  
  if (correlated){
    chol_vals <- split.coefs$coef
    Cholesky <- matrix(0, N_rand, N_rand)
    counter = 1
    for (i in 1:N_rand){
      for (j in 1:N_rand){
        if (j<=i){
          Cholesky[j,i] <- chol_vals[counter]
          counter <- counter + 1
        }
      }
    }
    
    Covariance <- t(Cholesky) %*% Cholesky
    Correlation <- cov2cor(Covariance)
    
    fit$Cholesky <- Cholesky
    fit$Covariance <- Covariance
    fit$Correlation <- Correlation
    fit$sd <- as.vector(sqrt(diag(Covariance)))
  }
  else{
    fit$sd <- abs(split.coefs$rpr)
  }
  
  names(fit$estimate) <- x_names
  fit$coefs <- coefs[1:(N_fixed+N_rand)]
  fit$estimate <- ifelse(grepl("St.", names(fit$estimate)), abs(fit$estimate ), fit$estimate)
  fit$formula <- formula
  fit$x_names <- x_names
  fit$rpar_formula <- rpar_formula
  fit$scrambled <- scrambled
  fit$numdraws <- ndraws
  fit$correlated <- correlated
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

  
  fit$modelType <- "rppois"
  .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
}
