#' flexCountReg Class
#'
#' A class to represent various objects created by the flexCountReg package.
#'
#' @slot model The fitted model object (can be any type).
#' @slot data The data used for fitting the model (data frame).
#' @slot formula The R formula used in the regression model (main formula, not including random parameters, etc.)
#' @slot call The matched call (language).
#' @slot additional Any additional information (list).
#' @importFrom methods setClass
#' 
#' @export
setClass(
  Class = "flexCountReg",
  slots = list(
    model = "ANY",
    data = "data.frame",
    formula="formula",
    call = "language",
    additional = "list"
  )
)

