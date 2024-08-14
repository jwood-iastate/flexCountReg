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
#' \deqn{PMF=\frac{\Gamma(\alpha+\gamma)\Gamma(\alpha+\rho+x)\Gamma(\rho+\gamma)}{\Gamma(\alpha+\rho)\Gamma(\gamma+\rho+x)\Gamma(\alpha+\gamma+\rho+x)}\left(\frac{\rho}{\alpha+\rho}\right)}
#' \deqn{\mu=e^{X\beta}=\frac{\alpha\rho}{\gamma(\alpha+\rho)}}
#' \deqn{\sigma^2=\frac{\alpha\rho(\alpha+\rho+\gamma)}{\gamma^2(\alpha+\rho)^2(\alpha+\rho+1)}}
#'
#' Where \eqn{\alpha} and \eqn{\rho} are distribution parameters with the constraints that \eqn{\alpha\geq 0} and \eqn{\rho\geq 0}, \eqn{X} is a matrix of predictors, and \eqn{\beta} is a vector of coefficients.
#'
#' If we use:
#' \deqn{\gamma=\frac{\mu(\alpha+\rho)}{\alpha\rho}}
#' 
#' The PMF becomes:
#' 
#' \deqn{PMF=\frac{\Gamma(\alpha + \rho) \Gamma\left(\alpha + \frac{2\mu}{\rho} + x\right) \Gamma\left(\frac{2\mu}{\rho} + \rho\right)}{\Gamma\left(\alpha + \frac{2\mu}{\rho}\right) \Gamma\left(\frac{2\mu}{\rho} + \rho + x\right) \Gamma\left(\alpha + \rho + \frac{2\mu}{\rho} + x\right)} \left( \frac{\frac{2\mu}{\rho}}{\alpha + \frac{2\mu}{\rho}} \right)}
#' 
#' #' This results in a regression model where:
#' \deqn{\mu=e^{X\beta}}
#' \deqn{\sigma^2=\mu\left(1-\frac{1}{\alpha+\rho+1}\right)+\mu^2\frac{(\alpha+\rho)^2}{\alpha\rho(\alpha+\rho+1)}}
#'
#' @import maxLik stats
#' @importFrom MASS glm.nb
#' @examples
#' \donttest{
#'
#' # Generalized Waring Model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
#' genwaring.mod <- genWaring(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                                 ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads)
#' summary(genwaring.mod)}
#' @export
genWaring <- function(formula, data, method = 'BHHH', max.iters = 1000) {
  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(model.matrix(formula, data))
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- colnames(X)
  
  # Use the Poisson as starting values
  p_model <- glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  start <- append(start, 1) # add initial starting value for ln(alpha)
  full_start <- append(start, 1) # add initial starting value for ln(rho)
  
  gw_prob <- function(y, mu, alpha, rho){
    g <- 2*mu/rho
    p <- gamma(alpha + rho)*gamma(alpha+g+y)*gamma(g+rho)*(g/(alpha+g))/(gamma(alpha+g)*gamma(g+rho+y)*gamma(alpha+rho+g+y))
    return(p)
  }
  
  reg.run <- function(beta, y, X) {
    pars <- length(beta) - 2
    
    coefs <- as.vector(unlist(beta[1:pars]))
    dispars <- tail(beta, 2)
    alpha <- exp(dispars[1])
    rho <- exp(dispars[2])
    
    predicted <- exp(X %*% coefs)
    
    probs <- gw_prob(y, predicted, alpha, rho)
    
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
  fit$alpha <- exp(distpars[1])
  fit$rho <- exp(distpars[2])
  
  x_names <- append(x_names, 'ln(alpha)')
  x_names <- append(x_names, 'ln(rho)')
  for (i in 1:length(x_names)) {
    names(fit$estimate)[i] <- x_names[i]
  }
  
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
