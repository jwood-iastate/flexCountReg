#' Function for estimating a Poisson Lognormal regression model via maximum simulated likelihood
#'
#' @name poisLogn
#' @param formula an R formula.
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param data a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula}.
#' @param ndraws the number of Halton draws to use for estimating the random parameters,
#' @param max.iters the maximum number of iterations to allow the optimization method to perform,
#'
#' @import nlme  maxLik  MASS  stats  randtoolbox
#' @importFrom utils head  tail
#'
#' @examples
#' \donttest{
#'
#' ## Poisson-Lognormal Model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
#' poislog.mod <- poisLogn(Total_crashes ~ lnaadt + lnlength + speed50 + ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 ndraws=50)
#' summary(poislog.mod)}
#' @export
poisLogn <- function(formula, data, ndraws = 1500, method = 'BHHH', max.iters = 1000) {

  mod_df <- stats::model.frame(formula, data)
  X <- stats::model.matrix(formula, data)
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- nlme::Names(formula, data)

  p_model <- MASS::glm.nb(formula, data)
  start <- as.vector(p_model$coefficients)
  s <- sqrt(log(1/p_model$theta+1))
  start <- append(start, s) # Add starting values for sigma

  n <- as.integer(ndraws)

  poisson_prob <- function(observed, predicted) {
    log_probs <- dpois(observed, predicted, log = TRUE)

    return(log_probs)
  }

  p_poisson_lognormal <- function(p, y, X, n, est_method){

    coefs <- as.array(head(p,-1))
    sigma <- abs(tail(p,1))
    h = randtoolbox::halton(n)

    probs <- as.vector(rep(0, length(y)))

    i_lognorm <- stats::qnorm(h, mean=0, sd=sigma)

    for (i in i_lognorm){
      coefs_i <- coefs
      coefs_i[1] <- coefs_i[1] + i
      mu <- exp(X %*% coefs_i)
      p_prob <- poisson_prob(y, mu)
      probs <- probs + exp(p_prob)/n
    }

    ll <- sum(log(probs))
    if (est_method == 'bhhh' | est_method == 'BHHH'){
      return(log(probs))
    } else{return(ll)}
  }

  fit <- maxLik::maxLik(p_poisson_lognormal,
                start = start,
                y = y,
                X = X,
                n = n,
                est_method = method,
                method = method,
                control = list(iterlim = max.iters))

  beta_est <- fit$estimate
  beta_pred <- head(beta_est, -1)
  x_names <- append(x_names, 'sigma')
  names(fit$estimate) <- x_names
  fit$estimate['sigma'] <- abs(fit$estimate['sigma'])

  sigma_sq <- fit$estimate['sigma']^2

  fit$formula <- formula
  fit$beta_pred <- beta_pred
  fit$sigma <- fit$estimate['sigma']
  fit$modelType <- "posLogn"
  fit$predictions <- exp(X %*% beta_pred + sigma_sq/2)
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  # Note that this has the predictions, residuals, and observed outcome stored with the model

  return(fit)
}
