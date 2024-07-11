#' Poisson-Generalized-Exponential Regression
#' 
#' @name poisGE
#' @param formula an R formula.
#' @param data a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula}. This can be the data used for estimating the model or another dataframe,
#' @param ln.scale.formula an optional formula for using independent variables to estimate the natural log of the scale parameter.
#' @param ndraws the number of Halton draws to use for the integration over the lognormal distribution.
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param max.iters the maximum number of iterations to allow the optimization method to perform.
#' @param print.level determines the level of verbosity for printing details of the optimization as it is computed. A value of 0 does not print out any information, a value of 1 prints minimal information, and a value of 2 prints the most information.
#'
#' @details
#' This implements maximum simulated likelihood (MSL) to estimate a Poisson-Generalized Exponential regression model. The regression model has the flexibility to model the scale parameter as a function of independent variables, similar to the generalized negative binomial.
#'
#' Details of the distribution can be found with the function \code{\link{dpge}}
#' The t-statistics and p-values for the coefficients related to ln(sigma) are,
#' by default, testing if the coefficients are different from a value of 0. This
#' has little practical meaning given that they are coefficients for ln(sigma).
#' They are not testing if the coefficients have statistical significance in
#' terms of improvement over a Poisson model. The Likelihood-Ratio test results
#' provided in the output provide a test comparing if the Poisson-Lognormal
#' model provides a statistically significant improvement in model fit over the
#' Poisson model.
#'
#' @import maxLik  stats modelr
#' @importFrom MASS glm.nb
#' @include ppoisGE.R
#' @export
#' @examples
#' \donttest{
#'
#' ## Generalized Poisson-Lognormal
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
#'
#' poisge.mod <- poisGE(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                         ShouldWidth04 + AADTover10k,
#'                         ln.scale.formula = ~ speed50,
#'                         data=washington_roads, 
#'                         ndraws = 100, 
#'                         method = 'nm')
#' summary(poisge.mod)
#' }
poisGE <- function(formula, data,  ln.scale.formula = NULL, ndraws=1500, 
                   method = 'BHHH', max.iters=200, print.level=0) {
  
  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- colnames(X)
  
  # Use the Poisson as starting values
  p_model <- glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  
  start <- append(start, 0) # Append ln(shape)
  
  if (is.null(ln.scale.formula)){
    full_start <- append(start, 0) # Append ln(scale)
  }else{
    scale_X <- stats::model.matrix(ln.scale.formula, data)
    S <- as.matrix(scale_X)
    scale_names <- nlme::Names(ln.scale.formula, data)
    s_coefs <- rep(0, length(scale_names))
    full_start <- c(start, s_coefs)
  }
  
  modparams <- as.numeric(length(start)-1) # save the number of model coefficients, not including scale
  
  reg.run <- function(beta, y, X, est_method){
    
    params_split <- split(beta,ceiling(seq_along(beta) / modparams))
    
    coefs <- as.vector(unlist(params_split[1]))
    dist_params <- as.vector(unlist(params_split[2]))
    shape <- exp(dist_params[1])
    
    
    if (is.null(ln.scale.formula)){
      scale <- exp(dist_params[2])
    }else{
      scale_pars <- dist_params[2:length(dist_params)]
      scale <- exp(S %*% scale_pars)
    }
    
    predicted <- exp(X %*% coefs)
    
    probs <- dpge(x=y, mean=predicted, shape=shape, scale=scale, ndraws=ndraws)
    
    ll <- sum(log(probs))
    if (est_method == 'bhhh' | est_method == 'BHHH'){
      return(log(probs))
    } else{return(ll)}
  }
  
  fit <- maxLik::maxLik(reg.run,
                        start = full_start,
                        y = y,
                        X = X,
                        est_method = method,
                        method = method,
                        control = list(iterlim = max.iters, printLevel = print.level))
  
  beta_est <- fit$estimate
  
  betas <- split(beta_est,ceiling(seq_along(beta_est) / modparams))
  beta_pred <- as.vector(unlist(betas[1]))
  
  pred <- exp(X %*% beta_pred)
  
  beta_dist_params <- as.vector(unlist(betas[2]))
  shape <- exp(beta_dist_params[1])
  fit$shape <- shape
  if (is.null(ln.scale.formula)){
    scale <- exp(exp(beta_dist_params[2]))
    fit$scale <- scale
  }
  else{
    scale_coefs <- beta_dist_params[2:length(beta_dist_params)]
    fit$scale_coefs <- scale_coefs
    scales <- exp(S %*% scale_coefs)
    fit$scales <- scales
  }
  x_names <- append(x_names, 'ln(shape)')
  if (is.null(ln.scale.formula)){
    x_names <- append(x_names, 'ln(scale)')
  }else{
    x_names <- append(x_names, paste0('ln(scale):', scale_names))
  }
  
  names(fit$estimate) <- x_names
  fit$beta_pred <- beta_pred # save coefficients for predictions
  fit$formula <- formula
  fit$ln.scale.formula <- ln.scale.formula
  fit$modelType <- "Poisson-Generalized-Exponential"
  fit$predictions <- pred
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum # The log-likelihood of the model
  
  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}