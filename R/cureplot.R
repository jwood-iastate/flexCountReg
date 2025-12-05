#' Cumulative Residuals (CURE) Plot for Count Models
#'
#' This function generates a Cumulative Residuals (CURE) plot for count models,
#' including those with random parameters, estimated using the flexCountReg
#' package.
#'
#' @name cureplot
#' @param model A model object estimated using this R package.
#' @param data Optional dataframe. If not provided, the data used to fit the
#'    model will be used.
#' @param indvar Optional independent variable name (character string). This is
#'   the continuous independent variable to plot the cumulative residuals
#'   against. If not provided, the plot will be against the predicted values.
#' @param method Optional parameter to pass to the predict function. This is
#'   only used for random parameters models (e.g., "Simulated" or "Individual").
#'   For further details, see \code{\link{predict.flexCountReg}}.
#' @param n_resamples Number of resamples for potential resampling in the CURE
#'   plot confidence bands. Default is 0 (no bands).
#' @param ... Additional arguments passed to \code{\link[cureplots]{cure_plot}}.
#' 
#' @importFrom cureplots calculate_cure_dataframe cure_plot
#' @importFrom stats model.frame model.response formula predict
#' @export
#' @examples
#' \donttest{
#' ## Example using a Negative Binomial model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
#' 
#' nb_model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + 
#'                             ShouldWidth04 + AADTover10k,
#'                             data = washington_roads, family = 'nb2', 
#'                             method = 'NM', max.iters = 500)
#'                             
#' # 1. Plot against fitted values (default) with confidence bands
#' cureplot(nb_model, n_resamples = 20)
#' 
#' # 2. Plot against a specific covariate (e.g., lnlength)
#' cureplot(nb_model, indvar = "lnlength", n_resamples = 20)
#' }
#' 
cureplot <- function(model, data = NULL, indvar = NULL, 
                     method = "Simulated", n_resamples = 0, ...) {
  
  # 1. Handle Data
  if (is.null(data)) {
    data <- model$data
  } else {
    data <- as.data.frame(data)
  }
  
  # 2. Extract Response Variable (y)
  # We use the formula from the model object to ensure we get the correct
  # response even if transformations were used in the formula, though strict
  # matching is safer.
  mf <- stats::model.frame(model$formula, data = data)
  y <- stats::model.response(mf)
  
  # 3. Generate Predictions
  # This relies on predict.flexCountReg handling the class/method logic
  predictions <- stats::predict(model, data = data, method = method)
  
  if(length(predictions) != length(y)){
    msg <- paste0("Length of predictions does not match length of observed", 
                  "data. Check for missing values.")
    warning(msg)
  }
  
  # 4. Calculate Raw Residuals
  resids <- y - predictions
  
  # 5. Determine Independent Variable (X-axis)
  if (!is.null(indvar)) {
    if (!indvar %in% names(data)) {
      warning(paste("Variable", indvar, "not found in the provided data."))
    }
    indvar_values <- data[[indvar]]
    x_name <- indvar
  } else {
    indvar_values <- predictions
    x_name <- "Predicted Values"
  }
  
  # 6. Create CURE Dataframe
  # Using cureplots package structure
  cure.df <- cureplots::calculate_cure_dataframe(covariate = indvar_values, 
                                                 residuals = resids)
  
  # Rename the first column to the variable name so the plot label is correct
  names(cure.df)[1] <- x_name
  
  # 7. Generate Plot
  cureplots::cure_plot(cure.df, n_resamples = n_resamples, ...)
}