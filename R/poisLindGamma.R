#' Function for estimating a Poisson-Lindley-Gamma (i.e., Negative Binomial-Lindley) regression model
#'
#' @name poisLindGamma
#' @param formula an R formula.
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}}.
#' @param data a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula}.
#' @param ndraws the number of Halton draws to use for the integration over the gamma distribution.
#' @param max.iters the maximum number of iterations to allow the optimization method to perform.
#'
#' @import nlme  maxLik  MASS  stats modelr
#' @include plindGamma.R
#'
#' @details
#' The Poisson-Lindley-Gamma regression is based on a compound Poisson-Lindley-Gamma distribution. Details of the distribution can be seen at \code{\link[flexCountReg]{dplindGamma}}.
#'
#' The mean for the regression model is:
#' \deqn{\mu=e^{X\beta}}
#'
#' The variance function is defined as:
#' \deqn{\sigma^2=\mu+\left(\alpha+1-\frac{2}{(\theta+2)^2}\right)\mu^2}
#'
#' It should be noted that the p-value for the parameters `ln(theta)` and `ln(alpha)` in the model summary are testing if the parameter `theta` and `alpha` are equal to a value of 1.
#'
#'
#' @examples
#' \donttest{
#'
#' ## Poisson-Lindley-Gamma Model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
#' poislindgamma.mod <- poisLindGamma(Animal ~ lnaadt + lnlength + speed50 +
#'                                 ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 method="nm",
#'                                 ndraws=100,
#'                                 max.iters = 1000)
#' summary(poislindgamma.mod)}
#' @export
poisLindGamma <- function(formula, data, method = 'BHHH', ndraws=1500, max.iters = 1000) {

  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- colnames(X)

  # Use the Poisson as starting values
  p_model <- MASS::glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  start <- append(start, 0) # inital log(theta)
  full_start <- append(start, 0) # initial value for log(alpha)

  modparams <- as.numeric(length(full_start))

  reg.run <- function(beta, y, X, est_method){
    pars <- length(beta)-2
    distpars <- as.vector(unlist(beta[(pars+1):length(beta)]))

    coefs <- as.vector(unlist(beta[1:pars]))
    theta <- exp(distpars[1])
    alpha <- exp(distpars[2])

    predicted <- exp(X %*% coefs)

    probs <- dplindGamma(y, mean=predicted, theta=theta, alpha=alpha, ndraws=ndraws)

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
  npars <- length(beta_est)-2
  beta_pred <- as.vector(unlist(beta_est[1:npars]))
  distpars <- as.vector(unlist(beta_est[(npars+1):length(beta_est)]))

  fit$beta_pred <- beta_pred # save coefficients for predictions

  mu <- exp(X %*% beta_pred)
  fit$theta <- exp(distpars[1])
  fit$alpha <- exp(distpars[2])

  x_names <- append(x_names, 'ln(theta)')
  x_names <- append(x_names, 'ln(alpha)')
  names(fit$estimate) <- x_names

  fit$predictions <- mu

  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum # The log-likelihood of the model
  fit$modelType <- "poisLindGamma"

  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}
