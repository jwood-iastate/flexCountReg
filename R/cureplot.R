#' Cumulative Residuals (CURE) Plot for Count Models
#'
#' This function generates a Cumulative Residuals (CURE) plot for count models,
#' including those with random parameters, estimated using the flexCountReg
#' package.
#'
#' @name cureplot
#' @param model A model object estimated using this R package.
#' @param data Optional dataframe. If not provided, the rows used to fit the
#'    model are used. Rows with missing values in the outcome, the
#'    predictors or \code{indvar} are left out.
#' @param indvar Optional independent variable name (character string). This is
#'   the continuous independent variable to plot the cumulative residuals
#'   against. If not provided, the plot will be against the predicted values.
#' @param method Optional prediction method passed to \code{predict()} for
#'   random-parameter models (e.g. "Simulated", "Exact" or "Individual") and
#'   negative multinomial models ("Standard" or "EB"). The default
#'   (\code{NULL}) uses the estimator's default. For further details, see
#'   \code{\link{predict.flexCountReg}}.
#' @param n_resamples Number of curves of randomly reordered residuals drawn
#'   in grey for comparison. Default is 0 (none). The red limits
#'   (\eqn{\pm 1.96} times the standard deviation of the cumulative residuals
#'   computed by \pkg{cureplots}) are drawn whatever the value of
#'   `n_resamples`; they are a descriptive reference, not a calibrated
#'   simultaneous band.
#' @param ... Additional arguments passed to \code{\link[cureplots]{cure_plot}}.
#' 
#' @importFrom cureplots calculate_cure_dataframe cure_plot
#' @importFrom stats model.frame model.response formula predict
#' 
#' @returns A CURE plot generated with \pkg{cureplots}.
#' 
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
                     method = NULL, n_resamples = 0, ...) {
  
  # 1. Data: the rows used in estimation by default
  data <- if (is.null(data)) .fcr_training_rows(model) else
    as.data.frame(data)
  if (!is.null(indvar) && !indvar %in% names(data)) {
    stop("Variable \"", indvar, "\" not found in the data.", call. = FALSE)
  }
  
  # 2. Outcome and predictions, one value per row of `data` (NA for rows
  # with missing values), so they stay aligned. `method` is passed only for
  # random-parameter and negative multinomial models.
  y <- .fcr_response(model, data)
  if (is.null(y)) {
    stop("The data do not contain the model's outcome.", call. = FALSE)
  }
  predictions <- if (is.null(method)) {
    stats::predict(model, newdata = data)
  } else {
    if (!inherits(model, c("countreg_rp", "negmultinom"))) {
      stop("`method` applies only to countreg_rp() and negmultinom() ",
           "models.", call. = FALSE)
    }
    stats::predict(model, newdata = data, method = method)
  }
  predictions <- as.numeric(predictions)
  if (length(predictions) != length(y)) {
    stop("Predictions and outcomes have different lengths.", call. = FALSE)
  }
  
  # 3. Plotting variable and one common mask of complete rows
  if (!is.null(indvar)) {
    indvar_values <- data[[indvar]]
    x_name <- indvar
  } else {
    indvar_values <- predictions
    x_name <- "Predicted Values"
  }
  ok <- !is.na(y) & !is.na(predictions) & !is.na(indvar_values)
  resids <- y[ok] - predictions[ok]
  indvar_values <- indvar_values[ok]
  
  # 4. CURE data frame and plot. cureplots names the covariate column after
  # the argument expression and reports it in a message, and cure_plot()
  # reports that a data frame was supplied; both messages describe this
  # wrapper's internal objects, so they are muffled (other messages pass).
  cure.df <- .fcr_muffle_cure_messages(cureplots::calculate_cure_dataframe(
    indvar_values, resids))
  names(cure.df)[1] <- x_name
  .fcr_muffle_cure_messages(
    cureplots::cure_plot(cure.df, n_resamples = n_resamples, ...))
}

# Evaluate a cureplots call without its informational messages about the
# covariate column and the supplied data frame.
.fcr_muffle_cure_messages <- function(expr) {
  withCallingHandlers(expr, message = function(m) {
    if (grepl("^(Covariate: |CURE data frame was provided)",
              conditionMessage(m))) {
      invokeRestart("muffleMessage")
    }
  })
}
