#' Function for estimating a Sichel regression model
#'
#' @name sichel
#' @param formula an R formula.
#' @param method a method to use for optimization in the maximum likelihood
#'   estimation. For options, see \code{\link[maxLik]{maxLik}}.
#' @param data a dataframe that has all of the variables in the \code{formula}.
#' @param max.iters the maximum number of iterations to allow the optimization
#'   method to perform.
#'
#' @import nlme  maxLik  MASS  stats modelr
#' @include psichel.R
#'
#' @details
#' The Sichel regression model is based on the Sichel Distribution .
#' 
#' The expected value of the distribution in the regression utilizes a log-link
#' function. Thus, the mean is: \deqn{\mu=e^{X\beta}}
#'
#' The variance is:
#' \deqn{
#'    \sigma^2= 
#'      \mu + 
#'      \left(\frac{2\sigma(\gamma+1)}{c} + \frac{1}{c^2} - 1\right) \mu^2}
#' 
#' The parameter \eqn{\sigma} is estimated as the natural logarithm transformed
#' value, ln(sigma), to ensure that \eqn{\sigma>0}.
#'
#' @seealso \code{\link{dsichel}()} for more information on the distribution.
#' 
#' @examples
#' \donttest{
#'
#' # Sichel Model
#' data("washington_roads")
#' 
#' sichel.mod <- sichel(Total_crashes ~ lnaadt + lnlength,
#'                                 data=washington_roads,
#'                                 method="NM",
#'                                 max.iters = 1000)
#' summary(sichel.mod)}
#' @export
sichel <- function(formula, data, method = 'BHHH', max.iters = 1000) {

  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- colnames(X)

  # Use the Negative Binomial as starting values
  p_model <- MASS::glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  lnsigma <- 0 # initial ln(sigma)
  gamma <- 1 # initial gamma

  full_start <- append(start, lnsigma)
  full_start <- append(full_start, gamma)

  modparams <- as.numeric(length(full_start))

  reg.run <- function(beta, y, X, est_method){
    pars <- length(beta)-2
    distpars <- tail(beta, 2)

    coefs <- as.vector(unlist(beta[1:pars]))
    sigma <- exp(unlist(distpars[1]))
    gamma <- unlist(distpars[2])

    predicted <- exp(X %*% coefs) 

    probs <- dsichel(y, mu=predicted, sigma=sigma, gamma=gamma)

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
  distpars <- tail(beta_est,2)
  fit$beta_pred <- beta_pred # save coefficients for predictions

  mu <- exp(X %*% beta_pred)
  fit$sigma <- exp(distpars[1])
  fit$gamma <- distpars[2]

  x_names <- append(x_names, 'ln(sigma)')
  x_names <- append(x_names, 'gamma')
  names(fit$estimate) <- x_names

  fit$predictions <- mu

  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum # The log-likelihood of the model
  fit$modelType <- "Sichel"

  # Estimate Poisson model for tests and pseudo R^2
  pois_mod <- glm(formula, data, family = poisson(link = "log"))
  base_mod <- glm(y ~ 1, family = poisson(link = "log"))

  LLpoisson <- sum(dpois(pois_mod$y, pois_mod$fitted.values, log=TRUE))
  LLbase <- sum(dpois(base_mod$y, base_mod$fitted.values, log=TRUE))

  fit$LR <- -2*(LLpoisson - fit$LL) # LR Statistic
  fit$LRdof <- 
    length(x_names) - length(pois_mod$coefficients) # LR Degrees of Freedom
  if (fit$LR > 0) {
    fit$LR_pvalue <- pchisq(fit$LR, fit$LRdof, lower.tail=FALSE)  # LR p-Value
  }else{
    fit$LR_pvalue <- 1
  }

  # Compute McFadden's Pseudo R^2, based on a Poisson intercept-only model
  fit$PseudoR2 <- 1-fit$LL/LLbase


  # Print out key model metrics
  LRpval <- ifelse(fit$LR_pvalue<0.0001, "<0.0001", round(fit$LR_pvalue,4))
  msg <- paste('The Likelihood Ratio (LR) Test for H0:',
               'Sichel is No Better than the Poisson')
  print(msg)
  print(paste('LR = ', round(fit$LR,4)))
  print(paste('LR degrees of freedom = ', fit$LRdof))
  print(paste('LR p-value = ', LRpval))
  print(paste("Macfadden's Pseudo R^2 = ", round(fit$PseudoR2,4)))

  return(fit)
}
