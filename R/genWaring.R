#' Function for estimating a Generalized Waring regression model
#'
#' @name genWaring
#' @param formula an R formula.
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see `maxLik` documentation.
#' @param data a dataframe that has all of the variables in  `formula` and `rpar_formula`,
#' @param max.iters the maximum number of iterations to allow the optimization method to perform.
#'
#' @details
#' ## Generalized Waring Probability Mass Function, Mean, and Variance
#' The following are the versions of the PMF, mean, and variance used in this function. This is adjusted from the typical formulation by replacing parameter \code{k} with \eqn{\mu}
#' \deqn{PMF=\frac{\Gamma(\alpha+\rho)}{\Gamma(\alpha)\Gamma(\rho)}\left(\frac{\alpha}{\mu(\rho-1)}\right)^{\alpha-1}\left(1+\frac{\alpha}{\mu(\rho-1)}\right)^{-(\alpha+\rho)}}
#' \deqn{\mu=e^{X\beta}}
#' \deqn{\sigma^2=\mu+\left(\frac{\frac{\mu(\rho-1)}{\alpha}+1}{\rho-2}\right)\mu+\left(\frac{\frac{\mu(\rho-1)}{\alpha}+\rho+1}{\frac{\mu(\rho-1)}{\alpha}(\rho-2)}\right)\mu^2}
#' Where \eqn{\alpha} and \eqn{\rho} are distribution parameters with the constraints that \eqn{\alpha\geq 0} and \eqn{\rho\geq 0}, \eqn{X} is a matrix of predictors, and \eqn{\beta} is a vector of coefficients.
#'
#'@import nlme  maxLik  MASS stats
#' @examples
#' \donttest{
#'
#' ## Generalized Waring Model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
#' genwaring.mod <- genWaring(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                                 ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads)
#' summary(genwaring.mod)}
#' @export
genWaring <- function(formula, data, method = 'BHHH', max.iters = 1000) {
  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- colnames(X)

  # Use the Poisson as starting values
  p_model <- MASS::glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  start <- append(start, 0) # add initial starting value for ln(alpha)
  full_start <- append(start, 0) # add initial starting value for ln(rho)

  modparams <- as.numeric(length(full_start))

  gw_prob <- function(y, mu, alpha, Rho){
    v <- alpha/(mu*(Rho-1))
    p <- gamma(alpha+Rho)/(gamma(alpha)*gamma(Rho)) * v^(alpha-1) * (1+v)^(-1*(alpha+Rho))
    return(p)
  }

  reg.run <- function(beta, y, X, est_method){
    pars <- length(beta)-2

    coefs <- as.vector(unlist(beta[1:pars]))
    dispars <- tail(beta,2)
    alpha <- exp(dispars[1])
    Rho <- exp(dispars[2])

    predicted <- exp(X %*% coefs) # not including Lindley adjustment

    probs <- gw_prob(y, predicted, alpha, Rho)

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
  fit$beta_pred <- beta_pred # save coefficients for predictions

  mu <- exp(X %*% beta_pred)
  distpars <- tail(beta_est,2)
  fit$alpha <- exp(distpars[1])
  fit$Rho <- exp(distpars[2])

  x_names <- append(x_names, 'ln(alpha)')
  x_names <- append(x_names, 'ln(rho)')
  for (i in 1:length(x_names)){
    names(fit$estimate)[i] <- x_names
  }
  names(fit$estimate) <- x_names

  fit$predictions <- mu

  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum # The log-likelihood of the model
  fit$modelType <- "genWaring"

  return(fit)
}
