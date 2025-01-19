#' Cumulative Residuals (CURE) Plot for Count Models
#'
#' This function generates a Cumulative Residuals (CURE) plot for count models,
#' including those with random parameters, estimated using the flexCountReg
#' package.
#'
#' @name cureplot
#' @param model A model object estimated using this R package.
#' @param data Optional dataframe. If not provided, the data used to fit the
#'   model will be used.
#' @param indvar Optional independent variable name. This is the continuous
#'   independent variable to plot the cumulative residuals against. If not
#'   provided, the plot will be against the predicted values.
#' @param method Optional parameter to pass to the predict function. This is
#'   only used for random parameters models. For further details, see
#'   \code{\link{predict.flexCountReg}}.
#' @param n_resamples Number of resamples for potential resampling in the CURE
#'   plot, based on the options in \code{\link[cureplots]{cure_plot}}.
#' @importFrom cureplots calculate_cure_dataframe cure_plot
#' @importFrom stats model.frame model.response formula
#' @export
#' @examples
#' 
#' ## Example using a Negative Binomial model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
#' nb_model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + 
#'                             ShouldWidth04 + AADTover10k,
#'                             data = washington_roads, family = 'nb2', 
#'                             method = 'NM', max.iters = 3000)
#' cureplot(nb_model)
#' 
cureplot <- function(
    model, data = NULL, indvar = NULL, method = "Simulated", n_resamples = 0) {

  # Use the modeling data if a dataframe is not provided
  if (is.null(data)) {
    data <- model$data
  }
  
  model.obj <- model$model  # The model object within the flexCountReg object
  formula <- formula(model)  # Get the model formula to extract the outcome
  
  # Extract the outcome
  y <- model.response(model.frame(formula, data))  
  
  predictions <- predict(model, data, method)
  resids <- y - predictions
  
  if (!is.null(indvar)) {
    indvar_values <- data[[indvar]]
    # assign(indvar, indvar_values, envir = .GlobalEnv)
  } else {
    indvar_values <- predictions
  }
  
  # Ensure that indvar_values and resids. have the same length
  if (length(indvar_values) != length(resids)) {
    stop("The length of the independent variable and residuals do not match.")
  }
  
  # Create CURE Dataframe
  cure.df <- cureplots::calculate_cure_dataframe(indvar_values, resids)
  names(cure.df)[1] <- ifelse(!is.null(indvar),indvar,"Predictions") # ensure naming of variable is correct
  
  # Generate CURE plot
  cureplots::cure_plot(cure.df, n_resamples = n_resamples)
}
