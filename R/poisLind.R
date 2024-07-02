#' Function for estimating a Poisson-Lindley regression model
#'
#' @name poisLind
#' @param formula an R formula.
#' @param method a method to use for optimization in the maximum likelihood 
#' estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param data a dataframe that has all of the variables in the \code{formula}.
#' @param max.iters the maximum number of iterations to allow the optimization 
#' method to perform,
#'
#' @import nlme  maxLik  MASS  stats modelr
#' @include plind.R
#'
#' @details
#' The Poisson-Lindley regression is based on a compound Poisson-Lindley 
#' distribution. It handles count outcomes with high levels of zero 
#' observations (or other high densities at low outcome values) that standard 
#' count regression methods, including the negative binomial, may struggle to 
#' adequately capture or model.
#'
#' The compound Probability Mass Function(PMF) for the Poisson-Lindley (PL) 
#' distribution is:
#' \deqn{f(y|\theta,\lambda)=\frac{\theta^2\lambda^y(\theta+\lambda+ 
#' y+1)}{(\theta+1)(\theta+\lambda)^{y+2}}}
#'
#' Where \eqn{\theta} and \eqn{\lambda} are distribution parameters with the 
#' restrictions that \eqn{\theta>0} and \eqn{\lambda>0}, and \eqn{y} is a 
#' non-negative integer.
#'
#' The expected value of the distribution is:
#' \deqn{\mu=\frac{\lambda(\theta+2)}{\theta(\theta+1)}}
#'
#' If a log-link function is used, the mean is:
#' \deqn{\mu=e^{X\beta}=\frac{\lambda(\theta+2)}{\theta(\theta+1)}}
#'
#' Thus, the parameter \eqn{\lambda} in the PL distribution when applied to 
#' regression analysis is:
#' \deqn{\lambda=\frac{\mu\theta(\theta+1)}{\theta+2}}
#' 
#' Using the replacement and simplifying results in:
#' \deqn{f(y \mid \theta, \mu) = \\frac{\theta^2 (\mu \theta (\theta+1))^y 
#' (\theta^2 (1+\mu) + \theta (2+\mu) + (\theta+2) (y+1))}{(\theta+1) 
#' (\theta+2)^{y+1} (\theta^2 (1+\mu) + \theta (2+\mu))^{y+2}}}
#' And
#' \deqn{LL=2 \log(\theta) + y (\log(\mu) + \log(\theta) + \log(\theta+1)) + 
#' \log(\theta^2 (1+\mu) + \theta (2+\mu) + (\theta+2) (y+1)) - \log(\theta+1) 
#' - (y+1) \log(\theta+2) - (y+2) \log(\theta^2 (1+\mu) + \theta (2+\mu))}
#'
#' The variance function is defined as:
#' \deqn{\sigma^2=\mu+\left(1-\frac{2}{(\theta+2)^2}\right)\mu^2}
#' 
#' It should be noted that the p-value for the parameter `ln(theta)` in the 
#' model summary is testing if the parameter `theta` is equal to a value of 1. 
#' This has no practical meaning. The Likelihood-Ratio (LR) test compares the 
#' Poisson-Lindley regression with a Poisson regression with the same 
#' independent variables. Thus, the PR test result indicates the statistical 
#' significance for the improvement in how well the model fits the data over a 
#' Poisson regression. This indicates the statistical significance of the 
#' `theta` parameter.
#'
#'
#' @examples
#' \donttest{
#'
#' ## Poisson-Lindley Model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
#' poislind.mod <- poisLind(Animal ~ lnaadt + lnlength + speed50 +
#'                                 ShouldWidth04 + AADTover10k,
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

  reg.run <- function(beta, y, X){
    pars <- length(beta)-1

    coefs <- as.vector(unlist(beta[1:pars]))
    theta <- exp(unlist(beta[length(beta)]))

    predicted <- exp(X %*% coefs) # not including Lindley adjustment

    probs <- dplind(y, mean=predicted, theta=theta)

    ll <- sum(log(probs))
    if (method == 'bhhh' | method == 'BHHH'){
      return(log(probs))
    } else{return(ll)}
  }
  
  # Define gradient and Hessian functions
  gradFun <- function(beta, y, X) {
    n.coefs <- length(beta)-1
    coefs <- beta[1:n.coefs]
    ln_theta <- unlist(beta[length(beta)])
    theta <- exp(ln_theta)
    
    mu <- exp(X %*% coefs)
    
    d_db_factor <- y+(theta^2*mu+theta*mu)/(theta^2*(1+mu)+theta*(2+mu)+
                      (theta+2)*(y+1)) - (y+2)*(theta^2*mu+theta*mu)/(theta^2*(1+mu)+theta*(2+mu))
    
    d_dlnT_1 <- 2/theta+y*(1/theta+1/(theta+1))
    d_dlnT_2 <- (2*theta*(1+mu)+2+mu+y+1)/(theta^2*(mu+1)+theta*(2+mu)+(theta+2)*(y+1))
    d_dlnT_3 <- -1/(theta+1)-(y+1)*(1/(theta+2))
    d_dlnT_4 <- -(y+2)*(2*theta*(1+mu)+2+mu)/(theta^2*(mu+1)+theta*(2+mu))
    
    d_dlnT <- d_dlnT_1+d_dlnT_2+d_dlnT_3+d_dlnT_4
    
    gradX <- X
    
    for (i in 1:n.coefs){
      gradX[,i] <- X[,i]*d_db_factor
    }

    return(cbind(gradX, d_dlnT))
  }

  fit <- maxLik::maxLik(reg.run,
                start = full_start,
                y = y,
                X = X,
                grad = if (method == 'BHHH') gradFun else function(...) colSums(gradFun(...)),
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

  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}
