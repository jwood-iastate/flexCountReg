#' Function for estimating a Poisson-Inverse-Gaussian regression model
#'
#' @name poisInvGaus
#' @param formula an R formula.
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param data a dataframe that has all of the variables in the \code{formula}.
#' @param form optional parameter indicating which formulation to use. Options include "Type 1" which is the standard form (and is the default) or "Type 2" which follows the formulation by Dean et. al. (1987).
#' @param max.iters the maximum number of iterations to allow the optimization method to perform.
#'
#' @import nlme  maxLik  MASS  stats modelr
#' @include pinvgaus.R
#'
#' @details
#' The Poisson-Inverse-Gaussian regression model is based on the Poisson-Inverse-Gaussian Distribution. 
#' @seealso [dpinvgaus()] for additional details of the distribution.
#' 
#' The expected value of the distribution in the regression utilizes a log-link function. Thus, the mean is:
#' \deqn{\mu=e^{X\beta}}
#'
#' The variance function for the Type 1 distribution (which is the default) is:
#' \deqn{\sigma^2=\mu+\eta\mu}
#' 
#' While the variance for the Type 2 distribution is:
#' \deqn{\sigma^2=\mu+\eta\mu^2}
#' 
#' The parameter \eqn{\eta} is estimated as the natural logarithm transformed value, ln(eta), to ensure that \eqn{\eta>0}.
#'
#' @examples
#' \donttest{
#'
#' ## Poisson-Inverse-Gaussian Model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
#' 
#' poisinvgaus.mod <- poisInvGaus(Total_crashes ~ lnaadt + lnlength + speed50 + 
#'                                 ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 method="nm",
#'                                 max.iters = 1000)
#' summary(poisinvgaus.mod)}
#' @export
poisInvGaus <- function(formula, data, form ="Type 1", method = 'BHHH', max.iters = 1000) {

  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- colnames(X)

  # Use the Negative Binomial as starting values
  p_model <- MASS::glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  eta <- 1/p_model$theta # intital eta

  full_start <- append(start, log(eta))

  modparams <- as.numeric(length(full_start))

  reg.run <- function(beta, y, X, est_method){
    pars <- length(beta)-1

    coefs <- as.vector(unlist(beta[1:pars]))
    eta <- exp(unlist(beta[length(beta)]))

    predicted <- exp(X %*% coefs) 

    probs <- dpinvgaus(y, predicted, eta, form)

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
                control = list(iterlim = max.iters))

  beta_est <- fit$estimate
  npars <- length(beta_est)-1

  beta_pred <- as.vector(unlist(beta_est[1:npars]))
  fit$beta_pred <- beta_pred # save coefficients for predictions

  mu <- exp(X %*% beta_pred)
  fit$eta <- exp(beta_est[length(beta_est)])

  x_names <- append(x_names, 'ln(eta)')
  names(fit$estimate) <- x_names

  fit$predictions <- mu

  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum # The log-likelihood of the model
  fit$modelType <- "poisInvGaus"
  fit$form <- form

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
  print('The Likelihood Ratio (LR) Test for H0: Poisson-Inverse-Gaussian is No Better than the Poisson')
  print(paste('LR = ', round(fit$LR,4)))
  print(paste('LR degrees of freedom = ', fit$LRdof))
  print(paste('LR p-value = ', LRpval))
  print(paste("Macfadden's Pseudo R^2 = ", round(fit$PseudoR2,4)))

  # Note that this has the predictions, residuals, observed outcome, LR test, and pseudo-r^2 stored with the model

  return(fit)
}
