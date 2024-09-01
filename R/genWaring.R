#' Function for estimating a Generalized Waring regression model
#'
#' @name genWaring
#' @param formula an R formula.
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}} documentation.
#' @param data a dataframe that has all of the variables in  `formula`.
#' @param max.iters the maximum number of iterations to allow the optimization method to perform.
#'
#' @details
#' # Generalized Waring Probability Mass Function, Mean, and Variance
#' The following are the versions of the PMF, mean, and variance used in this function. This is adjusted from the typical formulation by replacing parameter \code{k} with \eqn{\mu}
#' \deqn{PMF=\frac{\Gamma(\alpha+y)\Gamma(k+y)\Gamma(\rho+k)\Gamma(\alpha+\rho)}{y!\Gamma(\alpha)\Gamma(k)\Gamma(\rho)\Gamma(\alpha+k+\rho+y)}}
#' \deqn{\mu=e^{X\beta}=\frac{\alpha k}{\rho-1}}
#' \deqn{\sigma^2=\frac{\alpha k(\alpha+k+\rho-1)}{(\rho-1)^2(\rho-2)}}
#' 
#' The distribution parameters are often considered to capture the randomness (parameter \deqn{\alpha}), proneness (parameter \deqn{}k), and liability (parameter \deqn{\rho}) of the data.
#'
#' #' If we use:
#' \deqn{\alpha=\frac{\mu k}{\rho-1}}
#' 
#' The PMF becomes:
#' 
#' \deqn{PMF=\frac{\Gamma\left(\frac{\mu k}{\rho-1}+y\right)\Gamma(k+y)\Gamma(\rho+k)\Gamma\left(\frac{\mu k}{\rho-1}+\rho\right)}{y!\Gamma\left(\frac{\mu k}{\rho-1}\right)\Gamma(k)\Gamma(\rho)\Gamma\left(\frac{\mu k}{\rho-1}+k+\rho+y\right)}}
#' 
#' #' This results in a regression model where:
#' \deqn{\mu=e^{X\beta}}
#' \deqn{\sigma^2=\frac{\mu k^2\left(\frac{\mu k}{\rho-1}+k+\rho-1\right)}{(\rho-1)^3(\rho-2)}=\left(\frac{k^3+\rho k^2- k^2}{(\rho-1)^3(\rho-2)}\right)\mu+\left(\frac{k^3}{(\rho-1)^4(\rho-2)}\right)\mu^2}
#'
#' Note that when \deqn{p=1} or \deqn{p=2}, the distribution is undefined.
#' 
#' @import maxLik stats
#' @importFrom MASS glm.nb
#' @include Generalized-Waring.R
#' @examples
#' \donttest{
#'
#' # Generalized Waring Model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
#' genwaring.mod <- genWaring(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                                 ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 method='BFGS')
#' summary(genwaring.mod)}
#' @export
genWaring <- function(formula, data, method = 'BHHH', max.iters = 1000) {
  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(model.matrix(formula, data))
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- colnames(X)
  x_names <- append(x_names, 'ln(k)')
  x_names <- append(x_names, 'ln(rho)')
  
  # Use the Poisson as starting values
  p_model <- glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  start <- append(start, -0.5) # add initial starting value for ln(k)
  full_start <- append(start, -0.1) # add initial starting value for ln(rho)
  names(full_start) <- x_names

  reg.run <- function(beta, y, X) {
    pars <- length(beta) - 2
    
    coefs <- as.vector(unlist(beta[1:pars]))
    dispars <- tail(beta, 2)
    k <- exp(dispars[1])
    rho <- exp(dispars[2])
    
    predicted <- exp(X %*% coefs)
    
    probs <- dgwar(y, predicted, k, rho)
    
    ll <- sum(log(probs))
    if (method == 'bhhh' | method == 'BHHH') {
      return(log(probs))
    } else {
      return(ll)
    }
  }
  
  fit <- maxLik::maxLik(reg.run,
                        start = full_start,
                        y = y,
                        X = X,
                        method = method,
                        control = list(iterlim = max.iters))
  
  beta_est <- fit$estimate
  npars <- length(beta_est) - 2
  
  beta_pred <- as.vector(unlist(beta_est[1:npars]))
  fit$beta_pred <- beta_pred
  
  mu <- exp(X %*% beta_pred)
  distpars <- tail(beta_est, 2)
  fit$k <- exp(distpars[1])
  fit$rho <- exp(distpars[2])
  
  fit$predictions <- mu
  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum
  fit$modelType <- "genWaring"
  
  # Create and return a flexCountReg object
  obj <- .createFlexCountReg(model = fit, data = data, call = match.call())
  return(obj) 
}
