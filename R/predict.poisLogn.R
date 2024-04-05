#' Function for generating predictions based on the Poisson Lognormal model
#'
#' @name predict.poisLogn
#' @param model a Poisson-Lognormal model object,
#' @param data a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula}. This can be the data used for estimating the model or another dataframe,
#' @param method the method to be used in generating the predictions (options include \code{Simulated}, \code{Approximate}, or \code{Individual}).
#' @note the method option \code{Individual} requires that the outcome be observed for all observations.
#'
#' @import stats
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
#' summary(poislog.mod)
#'
#' predicted <- predict.poisLogn(poislog.mod, washington_roads)
#' print(predicted)}
#'
predict.poisLogn <- function(model, data){
  # This function takes in an nbg.reg model object and a dataframe
  # The function returns a dataframe with predictions, observed outcome, and residuals for the data that was input

  mod_df <- stats::model.frame(model$formula, data)
  X <- stats::model.matrix(model$formula, data)
  y <- as.numeric(stats::model.response(mod_df))

  sigma <- model$sigma
  beta_pred <- model$beta_pred

  predictions <- exp(X %*% beta_pred+sigma^2/2)
  return(predictions)
}
