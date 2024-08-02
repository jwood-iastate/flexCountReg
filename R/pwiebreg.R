#' Poisson-Weibull Regression
#'
#' Estimates a Poisson-Weibull regression model. 
#'
#' @name pwiebreg
#' @param formula A symbolic description of the model to be fitted, specifying the outcome
#'        and the independent variables with fixed effects.
#' @param alpha_formula A symbolic description of the model for the alpha parameter,
#'        specifying the distribution parameter as a function of predictor variables.
#' @param beta_formula A symbolic description of the model for the beta parameter,
#'        specifying the distribution parameter as a function of predictor variables.
#' @param data A dataframe containing all variables included in 'formula', 'alpha_formula', and 'beta_formula'.
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
#' @importFrom maxLik maxLik
#' @include poisWeib.R
#' 
#' @details
#' For the Poisson-Weibull Regression model, the expected values is:
#' \deqn{E[Y]=\lambda\beta\Gamma\left(1+\frac{1}{\alpha}\right)}
#' Where \eqn{\lambda} is the mean of the Poisson distribution, \eqn{\alpha} is the shape parameter, and \eqn{\beta} is the scale parameter.
#' 
#' To ensure that the regression model predicts the mean value, the regression utilizes:
#' \deqn{\mu=\exp{X\gamma}=\lambda\beta\Gamma\left(1+\frac{1}{\alpha}\right)}
#' Where \eqn{X} is a matrix of independent variables and \eqn{\gamma} is a vector of coefficients. 
#' 
#' This leads to:
#' \deqn{\lambda=\frac{\mu}{\beta\Gamma\left(1+\frac{1}{\alpha}\right)}}
#' 
#' The variance for the Poisson-Weibull regression is:
#' 
#' \deqn{V[Y]=\mu+\left(\frac{\Gamma\left(1+\frac{2}{\alpha}\right)}{\Gamma\left(1+\frac{1}{\alpha}\right)^2}-1\right)\mu^2}
#' 
#' @examples
#' data("washington_roads")
#' pw_rp <- pwiebreg(Total_crashes ~ lnlength + lnaadt,
#'                                  data = washington_roads,
#'                                  ndraws = 100)
#' print(summary(pw_rp))
#' 
#' @export
pwiebreg <- function(formula, alpha_formula = NULL, beta_formula = NULL, data,
                     ndraws = 1500,  method = 'BHHH', max.iters = 1000,
                     start.vals = NULL, print.level = 0, bootstraps = NULL) {
  
  # Prepare model matrices for fixed effects, alpha parameter, and beta parameter
  mod1_frame <- stats::model.frame(formula, data)
  X_Fixed <- stats::model.matrix(formula, data)
  
  if (!is.null(alpha_formula)) {
    mod_alpha_frame <- stats::model.frame(alpha_formula, data)
    X_alpha <- stats::model.matrix(alpha_formula, data)
  } else {
    X_alpha <- matrix(1, nrow(data), 1)  # Use an intercept-only model if alpha_formula is not provided
  }
  
  if (!is.null(beta_formula)) {
    mod_beta_frame <- stats::model.frame(beta_formula, data)
    X_beta <- stats::model.matrix(beta_formula, data)
  } else {
    X_beta <- matrix(1, nrow(data), 1)  # Use an intercept-only model if beta_formula is not provided
  }
  
  y <- stats::model.response(mod1_frame)
  
  # Define the main function for computing log-likelihood
  p_poisweibull_rp <- function(p, y, X_Fixed, X_alpha, X_beta, ndraws, 
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
    
    
    N_beta = ncol(X_beta)
    if (!is.null(beta_formula)){
      beta_coefs <- coefs[(N_fixed + N_alpha + 1):(N_fixed + N_alpha + N_beta)]
      beta <- exp(X_beta %*% beta_coefs)
    }
    else{
      beta <- exp(coefs[(N_fixed + N_alpha + 1)])
    }
    
    
    pred <- exp(X_Fixed %*% fixed_coefs)
    
    lambda <- pred / (beta * gamma(1 + 1 / alpha))
    
    # Call the provided Poisson-Weibull density function
    probs <- dpoisweibull(x = y, 
                          lambda = lambda, 
                          alpha = alpha,
                          beta = beta, 
                          log = FALSE)
    
    ll <- sum(log(probs))
    
    if (est_method == 'bhhh' || method == 'BHHH') {
      return(log(probs))
    } else {
      return(ll)
    }
  }
  
  # Initialize starting values if not provided
  start <- if (is.null(start.vals)) {
    rep(0, ncol(X_Fixed) + ncol(X_alpha) + ncol(X_beta))
  } else {
    start.vals
  }
  
  # Run the maximum likelihood estimation
  fit <- maxLik::maxLik(p_poisweibull_rp, start = start,
                        y = y, X_Fixed = X_Fixed,  X_alpha = X_alpha, X_beta = X_beta,
                        ndraws = ndraws, est_method = method,
                        method = method, control = list(iterlim = max.iters, printLevel = print.level))
  
  # Optionally, compute bootstrapped standard errors
  if (!is.null(bootstraps) & is.numeric(bootstraps)) {
    boot_err <- matrix(0, ncol = bootstraps, nrow = length(start))
    for (i in 1:bootstraps) {
      # Generate random sample of row indices with replacement - number of samples = number of observations
      sampled_indices <- sample(nrow(data), nrow(data), replace = TRUE)
      
      # Subset matrices and vector using sampled indices
      sampled_X_Fixed <- X_Fixed[sampled_indices, ]
      sampled_X_alpha <- if (!is.null(X_alpha)) X_alpha[sampled_indices, ] else NULL
      sampled_X_beta <- if (!is.null(X_beta)) X_beta[sampled_indices, ] else NULL
      sampled_y <- y[sampled_indices]
      
      model.boot <- maxLik::maxLik(p_poisweibull_rp,
                                   start = start,
                                   y = sampled_y,
                                   X_Fixed = sampled_X_Fixed,
                                   X_alpha = sampled_X_alpha,
                                   X_beta = sampled_X_beta,
                                   ndraws = ndraws,
                                   est_method = method,
                                   method = method,
                                   control = list(iterlim = max.iters, 
                                                  printLevel = print.level))
      boot_err[, i] <- fit$estimate - model.boot$estimate
    }
    
    stderr <- apply(boot_err, 1, sd)
    fit$bootstrapped_se <- stderr
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
  fit$beta_formula = beta_formula
  fit$ndraws = ndraws
  fit$bootstraps = if (!is.null(bootstraps)) bootstraps else NULL
  
  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}
