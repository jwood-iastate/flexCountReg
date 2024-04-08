#' Function for estimating a Poisson-Lindley regression model
#'
#' @name poisLind
#' @param formula an R formula.
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param data a dataframe that has all of the variables in the \code{formula}.
#' @param max.iters the maximum number of iterations to allow the optimization method to perform,
#'
#' @import nlme  maxLik  MASS  stats modelr
#' @include plind.R
#'
#' @details
#' The Poisson-Lindley regression is based on a compound Poisson-Lindley distribution. It handles count outcomes with high levels of zero observations (or other high densities at low outcome values) that standard count regression methods, including the negative binomial, may struggle to adequately capture or model.
#'
#' The compound Probability Mass Function(PMF) for the Poisson-Lindley (PL) distribution is:
#' \deqn{f(y|\theta,\lambda)=\frac{\theta^2\lambda^y(\theta+\lambda+ y+1)}{(\theta+1)(\theta+\lambda)^{y+2}}}
#'
#' Where \eqn{\theta} and \eqn{\lambda} are distribution parameters with the restrictions that \eqn{\theta>0} and \eqn{\lambda>0}, and \eqn{y} is a non-negative integer.
#'
#' The expected value of the distribution is:
#' \deqn{\mu=\frac{\lambda(\theta+2)}{\theta(\theta+1)}}
#'
#' If a log-link function is used, the mean is:
#' \deqn{\mu=e^{X\beta}=\frac{\lambda(\theta+2)}{\theta(\theta+1)}}
#'
#' Thus, the parameter \eqn{\lambda} in the PL distribution when applied to regression analysis is:
#' \deqn{\lambda=\frac{\mu\theta(\theta+1)}{\theta+2}}
#'
#' The variance function is defined as:#'
#' \deqn{\sigma^2=\mu+\left(1-\frac{2}{(\theta+2)^2}\right)\mu^2}
#'
#' It should be noted that the p-value for the parameter `ln(theta)` in the model summary is testing if the parameter `theta` is equal to a value of 1. This has no practical meaning. The Likelihood-Ratio (LR) test compares the Poisson-Lindley regression with a Poisson regression with the same independent variables. Thus, the PR test result indicates the statistical significance for the improvement in how well the model fits the data over a Poisson regression. This indicates the statistical significance of the `theta` parameter.
#'
#'
#' @examples
#' \donttest{
#'
#' ## Poisson-Lindley Model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
#' poislind.mod <- poisLind(Animal ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 method="nm",
#'                                 max.iters = 1000)
#' summary(poislind.mod)}
#' @export
poisLind <- function(formula, data, method = 'BHHH', max.iters = 1000) {

  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- colnames(X)

  # Use the Negative Binomial as starting values
  p_model <- MASS::glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  t <- 116761/exp(-9.04*p_model$theta)
  theta <- ifelse(t<100,t,1) # Approximate Initial Theta

  full_start <- append(start, log(theta))

  modparams <- as.numeric(length(full_start))

  reg.run <- function(beta, y, X, est_method){
    pars <- length(beta)-1

    coefs <- as.vector(unlist(beta[1:pars]))
    theta <- exp(unlist(beta[length(beta)]))

    predicted <- exp(X %*% coefs) # not including Lindley adjustment

    probs <- dplind(y, mean=predicted, theta=theta)

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
  fit$theta <- exp(beta_est[length(beta_est)])

  x_names <- append(x_names, 'ln(theta)')
  names(fit$estimate) <- x_names

  fit$predictions <- mu

  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum # The log-likelihood of the model
  fit$modelType <- "poisLind"

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
  print('The Likelihood Ratio (LR) Test for H0: Poisson-Lindley is No Better than the Poisson')
  print(paste('LR = ', round(fit$LR,4)))
  print(paste('LR degrees of freedom = ', fit$LRdof))
  print(paste('LR p-value = ', LRpval))
  print(paste("Macfadden's Pseudo R^2 = ", round(fit$PseudoR2,4)))

  # Note that this has the predictions, residuals, observed outcome, LR test, and pseudo-r^2 stored with the model

  return(fit)
}
