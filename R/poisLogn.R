#' Poisson-Lognormal Regression
#' 
#' @name poisLogn
#' @param formula an R formula.
#' @param data a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula}. This can be the data used for estimating the model or another dataframe,
#' @param ln.sigma.formula an optional formula for using independent variables to estimate the natural log of the standard deviation parameter (makes the model a generalized Poisson-Lognormal).
#' @param ndraws the number of Halton draws to use for the integration over the lognormal distribution.
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param max.iters the maximum number of iterations to allow the optimization method to perform.
#' @param print.level determines the level of verbosity for printing details of the optimization as it is computed. A value of 0 does not print out any information, a value of 1 prints minimal information, and a value of 2 prints the most information.
#'
#' @details
#' This implements maximum simulated likelihood (MSL) to estimate a Poisson-Lognormal regression model. The regression model has the flexibility to model the dispersion parameter \eqn{\sigma} from the lognormal distribution as a function of independent variables, similar to the generalized negative binomial.
#'
#' #' The compound Probability Mass Function(PMF) for the Poisson-Lognormal distribution is:
#' \deqn{f(y|\lambda,\theta,\alpha)=\int_0^\infty \frac{\lambda^y x^y e^{-\lambda x}}{y!}\frac{exp\left(-\frac{ln^2(x)}{2\sigma^2} \right)}{x\sigma\sqrt{2\pi}}dx}
#'
#' Where \eqn{\sigma} is a parameter for the lognormal distribution with the restriction \eqn{\sigma>0}, and \eqn{y} is a non-negative integer.
#'
#' The expected value of the distribution is:
#' \deqn{E[y]=e^{X\beta+\sigma^2/2} = \mu e^{\sigma^2/2}}
#' 
#' When `ln.sigma.formula` is used, the parameter \eqn{\sigma} is modeled as:
#' \deqn{ln(\sigma)=\beta_0+\beta_1 x_1 + \cdots + \beta_n x_n}
#' 
#' Thus, the resulting value for the parameter \eqn{\sigma} is:
#' \deqn{\sigma=e^{\beta_0+\beta_1 x_1 + \cdots + \beta_n x_n}}
#' 
#' The t-statistics and p-values for the coefficients related to ln(sigma) are, by default, testing if the coefficients are different from a value of 0. This has little practical meaning given that they are coefficients for ln(sigma). They are not testing if the coefficients have statistical significance in terms of improvement over a Poisson model. The Likelihood-Ratio test results provided in the output provide a test comparing if the Poisson-Longomal model provides a statistically significant improvement in model fit over the Poisson model.
#'
#' @import nlme  maxLik  MASS  stats modelr
#' @include ppoislogn.R
#' @export
#' @examples
#' \donttest{
#'
#' ## Generalized Poisson-Lognormal
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
#'
#' poslogn.mod <- poisLogn(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                         ShouldWidth04 + AADTover10k,
#'                         ln.sigma.formula = ~ lnaadt + AADTover10k,
#'                         data=washington_roads, 
#'                         ndraws = 100, 
#'                         method = 'BHHH')
#' summary(poslogn.mod)
#' }
poisLogn <- function(formula, data,  ln.sigma.formula = NULL, ndraws=1500, 
                     method = 'BHHH', max.iters=200, print.level=0) {
  
  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- colnames(X)
  
  # Use the Poisson as starting values
  p_model <- MASS::glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  s <- sqrt(log(1/p_model$theta+1))
  
  if (is.null(ln.sigma.formula)){
    if (!is.null(s)) {
      full_start <- append(start, s)
    }
  }else{
    sigma_X <- stats::model.matrix(ln.sigma.formula, data)
    sigma_names <- nlme::Names(ln.sigma.formula, data)
    s_coefs <- rep(0, length(sigma_names))
    s_coefs[1] = s
    if (!is.null(s_coefs)) {
      full_start <- c(start, s_coefs)
    }
  }
  
  modparams <- as.numeric(length(start)) # save the number of model coefficients, not including sigma
  
  reg.run <- function(beta, y, X, est_method){
    
    params_split <- split(beta,ceiling(seq_along(beta) / modparams))
    
    coefs <- as.vector(unlist(params_split[1]))
    dist_params <- as.vector(unlist(params_split[2]))
    
    if (is.null(ln.sigma.formula)){
      log_sigma <- dist_params[1]
    }else{
      sigma_coefs <- dist_params
      log_sigma <- sigma_X %*% sigma_coefs
    }
    
    sigma <- exp(log_sigma)
    
    predicted <- exp(X %*% coefs)
    
    # print(paste('The min value of the mean is ', min(predicted), " and the max is ", max(predicted)))
    
    # print(paste('The min value of the standard deviation is ', min(sigma), " and the max is ", max(sigma)))
    
    probs <- dpLnorm(x=y, mean=predicted, sigma=sigma, ndraws=ndraws)
    
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
  if (is.null(ln.sigma.formula)){
    sigma <- exp(beta_est[length(beta_est)])
    fit$sigma <- sigma
    pred <- exp(X %*% beta_pred + sigma^2/2)
  }
  else{
    sigma_pars <- as.vector(unlist(betas[2]))
    fit$sigma_coefs <- sigma_pars
    sigmas <- exp(sigma_X %*% sigma_pars)
    fit$sigma <- sigmas
    sigma_sq <- sigmas^2
    pred <- exp(X %*% beta_pred + sigma_sq/2)
  }
  
  if (is.null(ln.sigma.formula)){
    x_names <- append(x_names, 'ln(sigma)')
  }else{
    x_names <- append(x_names, paste0('ln(sigma):', sigma_names))
  }
  
  names(fit$estimate) <- x_names
  fit$beta_pred <- beta_pred # save coefficients for predictions
  fit$formula <- formula
  fit$ln.sigma.formula <- ln.sigma.formula
  fit$modelType <- "Poisson-Lognormal"
  fit$predictions <- pred
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum # The log-likelihood of the model
  
  # Estimate Poisson model for tests and pseudo R^2
  pois_mod <- glm(formula, data, family = poisson(link = "log"))
  base_mod <- glm(y ~ 1, family = poisson(link = "log"))
  
  LLpoisson <- sum(dpois(pois_mod$y, pois_mod$fitted.values, log=TRUE))
  LLbase <- sum(dpois(base_mod$y, base_mod$fitted.values, log=TRUE))
  
  fit$LR <- -2*(LLpoisson - fit$LL) # LR Statistic
  fit$LRdof <- length(x_names) - length(pois_mod$coefficients) # LR Degrees of Freedom
  if (fit$LR>0) {
    fit$LR_pvalue <- pchisq(fit$LR, fit$LRdof, lower.tail=FALSE)  # LR p-Value
  }else{
    fit$LR_pvalue <- 1
  }
  
  # Compute McFadden's Pseudo R^2, based on a Poisson intercept-only model
  fit$PseudoR2 <- 1-fit$LL/LLbase
  
  
  # Print out key model metrics
  LRpval <- ifelse(fit$LR_pvalue<0.0001, "<0.0001", round(fit$LR_pvalue,4))
  print('The Likelihood Ratio (LR) Test for H0: NB is No Better than Poisson')
  print(paste('LR = ', round(fit$LR,4)))
  print(paste('LR degrees of freedom = ', fit$LRdof))
  print(paste('LR p-value = ', LRpval))
  print(paste("Macfadden's Pseudo R^2 = ", round(fit$PseudoR2,4)))
  
  # Note that this has the predictions, residuals, observed outcome, LR test, and pseudo-r^2 stored with the model
  
  return(fit)
}
