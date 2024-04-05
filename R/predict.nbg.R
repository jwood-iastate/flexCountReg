#' Function for generating predictions based on the random parameters negative binomial with multiple optional methods
#'
#' @name predict.nbg
#' @param model a model object estimated using the \code{rpnb} function,
#' @param data a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula}. This can be the data used for estimating the model or another dataframe,
#' @param method the method to be used in generating the predictions (options include \code{Simulated}, \code{Approximate}, or \code{Individual}).
#' @note the method option \code{Individual} requires that the outcome be observed for all observations.
#'
#' @import stats
#'
#' @examples
#' \donttest{
#'
#' ## Random Parameters Negative Binomial model
#' data("washington_roads")
#' nb <- rpnb(Total_crashes ~ - 1 + lnlength + speed50,
#'            rpar_formula = ~ lnaadt,
#'            data = washington_roads)
#' summary(nb)
#'
#' rpnb.predictions <- predict.rpnb(nb, washington_roads, method=Individual)
#' print(rpnb.predictions)}
#' @export
predict.nbg <- function(model, data){
  # This function takes in an nbg.reg model object and a dataframe
  # The function returns a dataframe with predictions, observed outcome, and residuals for the data that was input

  mod_df <- stats::model.frame(model$formula, data)
  X <- stats::model.matrix(model$formula, data)
  y <- as.numeric(stats::model.response(mod_df))

  beta_pred <- model$beta_pred

  predictions <- exp(X %*% beta_pred)
  return(predictions)
}
