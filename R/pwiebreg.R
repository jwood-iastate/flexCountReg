#' Poisson-Weibull Regression
#'
#' Estimates a Poisson-Weibull regression model. 
#'
#' @name pwiebreg
#' @param formula A symbolic description of the model to be fitted, specifying the outcome
#'        and the independent variables with fixed effects.
#' @param alpha_formula A symbolic description of the model for the alpha parameter,
#'        specifying the distribution parameter as a function of predictor variables.
#' @param sigma_formula A symbolic description of the model for the sigma parameter,
#'        specifying the distribution parameter as a function of predictor variables.
#' @param data A dataframe containing all variables included in 'formula', 'alpha_formula', and 'sigma_formula'.
#' @param ndraws The number of Halton draws for integrating the Weibull distribution.
#' @param method Optimization method to be used for maximum likelihood estimation.
#'        See `maxLik` documentation for options.
#' @param max.iters Maximum number of iterations for the optimization method.
#' @param start.vals Optional vector of starting values for the optimization.
#' @param print.level Integer specifying the verbosity of output during optimization.
#' @param bootstraps Optional integer specifying the number of bootstrap samples to be used
#'        for estimating standard errors.
#' @import modelr        
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom purrr map map_df
#' @importFrom broom tidy
#' @importFrom dplyr group_by %>% summarise
#' @importFrom tibble deframe
#' @importFrom maxLik maxLik
#' @include poisWeib.R
#' 
#' @details
#' For the Poisson-Weibull Regression model, the expected values is:
#' \deqn{E[Y]=\lambda\sigma\Gamma\left(1+\frac{1}{\alpha}\right)}
#' Where \eqn{\lambda} is the mean of the Poisson distribution, \eqn{\alpha} is the shape parameter, and \eqn{\sigma} is the scale parameter.
#' 
#' To ensure that the regression model predicts the mean value, the regression utilizes:
#' \deqn{\mu=\exp{X\gamma}=\lambda\sigma\Gamma\left(1+\frac{1}{\alpha}\right)}
#' Where \eqn{X} is a matrix of independent variables and \eqn{\gamma} is a vector of coefficients. 
#' 
#' This leads to:
#' \deqn{\lambda=\frac{\mu}{\sigma\Gamma\left(1+\frac{1}{\alpha}\right)}}
#' 
#' The variance for the Poisson-Weibull regression is:
#' 
#' \deqn{V[Y]=\mu+\left(\frac{\Gamma\left(1+\frac{2}{\alpha}\right)}{\Gamma\left(1+\frac{1}{\alpha}\right)^2}-1\right)\mu^2}
#' 
#' @examples
#' data("washington_roads")
#' pw <- pwiebreg(Total_crashes ~ offset(lnaadt) + lnlength,
#'                 ndraws = 1500,
#'                 data = washington_roads,
#'                 alpha_formula = ~ -1 + lnaadt,
#'                 sigma_formula = ~ lnaadt,
#'                 method = 'NM',
#'                 bootstraps=30)
#' print(summary(pw))
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @export
pwiebreg <- function(formula, alpha_formula = NULL, sigma_formula = NULL, data,
                     ndraws = 1500,  method = 'NM', max.iters = 1000,
                     start.vals = NULL, print.level = 0, bootstraps = NULL) {
  
  # Prepare model matrices for fixed effects, alpha parameter, and sigma parameter
  mod1_frame <- stats::model.frame(formula, data)
  X_Fixed <- stats::model.matrix(formula, data)
  
  x_names <- colnames(X_Fixed)
  
  if (!is.null(alpha_formula)) {
    mod_alpha_frame <- stats::model.frame(alpha_formula, data)
    X_alpha <- stats::model.matrix(alpha_formula, data)
    x_names <- c(x_names, paste0("ln(alpha):",colnames(X_alpha)))
  } else {
    X_alpha <- matrix(1, nrow(data), 1)  # Use an intercept-only model if alpha_formula is not provided
    x_names <- append(x_names, "ln(alpha)")
  }
  
  if (!is.null(sigma_formula)) {
    mod_sigma_frame <- stats::model.frame(sigma_formula, data)
    X_sigma <- stats::model.matrix(sigma_formula, data)
    x_names <- c(x_names, paste0("ln(sigma):",colnames(X_sigma)))
  } else {
    X_sigma <- matrix(1, nrow(data), 1)  # Use an intercept-only model if sigma_formula is not provided
    x_names <- append(x_names, "ln(sigma)")
  }
  
  y <- stats::model.response(mod1_frame)
  
  # Generate Halton draws to use as quantile values
  h <- randtoolbox::halton(ndraws)
  
  # Efficient computatation of probabilities using C++ code
  dpoisweibull_cpp <- Vectorize(function(x, lambda, alpha , sigma, h) {
    p <- dpWeib_cpp(x, lambda, alpha, sigma, h)
    return(max(p,1e-10))})
  
  # Define the main function for computing log-likelihood
  p_poisweibull <- function(p, y, X_Fixed, X_alpha, X_sigma, ndraws, 
                               est_method) {
    N_fixed = ncol(X_Fixed)
    coefs <- as.array(p)
    fixed_coefs <- head(coefs, N_fixed)
    
    N_alpha = ncol(X_alpha)
    if (!is.null(alpha_formula)){
      alpha_coefs <- coefs[(N_fixed + 1):(N_fixed + N_alpha)]
      alpha <- exp(X_alpha %*% alpha_coefs)
    }
    else{
      alpha <- exp(coefs[(N_fixed + 1)])
    }
    
    
    N_sigma = ncol(X_sigma)
    if (!is.null(sigma_formula)){
      sigma_coefs <- coefs[(N_fixed + N_alpha + 1):(N_fixed + N_alpha + N_sigma)]
      sigma <- exp(X_sigma %*% sigma_coefs)
    }
    else{
      sigma <- exp(coefs[(N_fixed + N_alpha + 1)])
    }
    
    pred <- exp(X_Fixed %*% fixed_coefs)

    lambda <- pred / (sigma * gamma(1 + 1 / alpha))

    probs <- dpoisweibull_cpp(x = y, 
                          lambda = lambda, 
                          alpha = alpha,
                          sigma = sigma, 
                          h = h)
    
    ll <- sum(log(probs))
    
    if (est_method == 'bhhh' || method == 'BHHH') {
      return(log(probs))
    } else {
      return(ll)
    }
  }
  
  # Initialize starting values if not provided
  start <- if (is.null(start.vals)) {
    
    # Use the NB2 from MASS as starting values
    p_model <- glm.nb(formula, data = data)
    start <- unlist(p_model$coefficients)
    start <- c(start, rep(0, ncol(X_alpha) + ncol(X_sigma)))
  } else {
    start.vals
  }
  
  names(start) <- x_names
  
  # Run the maximum likelihood estimation
  fit <- maxLik::maxLik(p_poisweibull, start = start,
                        y = y, X_Fixed = X_Fixed, 
                        X_alpha = X_alpha, X_sigma = X_sigma,
                        ndraws = ndraws, est_method = method,
                        method = method, 
                        control = list(iterlim = max.iters, 
                                       printLevel = print.level))
  
  # Optionally, compute bootstrapped standard errors
  # create function to clean data and run maxLik
  weib.boot <- function(data){
    mod1_frame <- stats::model.frame(formula, data)
    X_Fixed <- stats::model.matrix(formula, data)
    y <- stats::model.response(mod1_frame)
    
    if (!is.null(alpha_formula)) {
      mod_alpha_frame <- stats::model.frame(alpha_formula, data)
    } else {
      X_alpha <- matrix(1, nrow(data), 1)  # Use an intercept-only model if alpha_formula is not provided
    }
    
    if (!is.null(sigma_formula)) {
      mod_sigma_frame <- stats::model.frame(sigma_formula, data)
      X_sigma <- stats::model.matrix(sigma_formula, data)
    } else {
      X_sigma <- matrix(1, nrow(data), 1)  # Use an intercept-only model if sigma_formula is not provided
    }
    
    int_res <- maxLik::maxLik(p_poisweibull, start = fit$estimate,
                              y = y, X_Fixed = X_Fixed,  X_alpha = X_alpha, X_sigma = X_sigma,
                              ndraws = ndraws, est_method = method,
                              method = method, control = list(iterlim = max.iters, printLevel = print.level))
    return(int_res)
  }
  
  
  if (!is.null(bootstraps) & is.numeric(bootstraps)) {
    bs.data <- modelr::bootstrap(data, n = bootstraps)
    
    
    mod1_frame <- stats::model.frame(formula, data)
    X_Fixed <- stats::model.matrix(formula, data)
    y <- stats::model.response(mod1_frame)
   
    models <- map(bs.data$strap, ~ weib.boot(data = .))
    tidied <- map_df(models, broom::tidy, .id = "id")
    
    SE <- tidied %>%
      group_by(term) %>%
      summarise(sd = sd(estimate)) %>% deframe()
    
    fit$bootstrapped_se <- SE
  }
  
  fit$coefficients = fit$estimate
  fit$se = if (!is.null(bootstraps) & is.numeric(bootstraps)) fit$bootstrapped_se else sqrt(diag(fit$hessian))
  fit$logLik = fit$maximum
  fit$converged = fit$convergence
  fit$model = "rpoisweibull_regression"
  fit$method = method
  fit$data = data
  fit$formula = formula
  fit$alpha_formula = alpha_formula
  fit$sigma_formula = sigma_formula
  fit$ndraws = ndraws
  fit$bootstraps = if (!is.null(bootstraps)) bootstraps else NULL
  
  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}
