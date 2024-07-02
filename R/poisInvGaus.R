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

  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}
