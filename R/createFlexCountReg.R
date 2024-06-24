# Internal constructor function for creating flexCountReg objects.
.createFlexCountReg <- function(model = NULL, data = NULL, call = NULL, ...) {
  formula <- formula(as.formula(model$formula))

  structure(list(
    call = call,
    model = model,
    data = data,
    formula=formula
  ), class = "flexCountReg")
}
