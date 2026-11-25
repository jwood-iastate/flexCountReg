#' Calculate Akaike Information Criterion (AIC)
#'
#' This function calculates the Akaike Information Criterion (AIC) for a given
#' model.
#'
#' @param LL Numeric value representing the log-likelihood of the model.
#' @param nparam Numeric value representing the number of parameters in the
#'   model.
#' @returns Numeric value representing the AIC.
#' @details For a fitted flexCountReg model, use \code{AIC(model)} instead,
#'   which takes the log-likelihood and number of parameters from the model.
#'
#' The AIC is calculated using the formula:
#' \deqn{AIC = -2 \cdot LL + 2 \cdot nparam}
#' Where \eqn{LL} is the log-likelihood of the model and \eqn{nparam} is the
#' number of parameters.
#' @examples
#' LL <- -120.5
#' nparam <- 5
#' myAIC(LL, nparam)
#' 
#' @export
myAIC <- function(LL, nparam){
  return(-2 * LL + 2 * nparam)
}

#' Calculate Bayesian Information Criterion (BIC)
#'
#' This function calculates the Bayesian Information Criterion (BIC) for a given
#' model.
#'
#' @param LL Numeric value representing the log-likelihood of the model.
#' @param nparam Numeric value representing the number of parameters in the
#'   model.
#' @param n Numeric value representing the number of observations.
#' @return Numeric value representing the BIC.
#' @details For a fitted flexCountReg model, use \code{BIC(model)} instead,
#'   which also uses the number of panels for panel models.
#'
#' The BIC is calculated using the formula:
#' \deqn{BIC = -2 \cdot LL + nparam \cdot \log(n)}
#' Where \eqn{LL} is the log-likelihood of the model, \eqn{nparam} is the number
#' of parameters, and \eqn{n} is the number of observations.
#' @examples
#' LL <- -120.5
#' nparam <- 5
#' n <- 100
#' myBIC(LL, nparam, n)
#' 
#' @export
myBIC <- function(LL, nparam, n){
  return(-2 * LL + nparam * log(n))
}

#' Calculate Root Mean Squared Error (RMSE)
#'
#' This function calculates the Root Mean Squared Error (RMSE) between observed
#' and predicted values.
#'
#' @param y Numeric vector representing the observed values, or a data frame
#'   (see the section on yardstick).
#' @param mu Numeric vector representing the predicted values (with a data
#'   frame, the column of observed values).
#' @param ... with a data frame, further arguments of the yardstick metric
#'   (the column of predictions first).
#' @return Numeric value representing the RMSE (with a data frame, the
#'   yardstick tibble).
#' @section yardstick: With a data frame, \code{rmse(data, truth, estimate)}
#'   calls the \pkg{yardstick} metric of the same name, and the function can
#'   be used in \code{yardstick::metric_set()}. Loading flexCountReg after
#'   tidymodels (which masks yardstick's function) therefore does not change
#'   tidymodels code.
#' @details The RMSE is calculated using the formula:
#' \deqn{RMSE = \sqrt{\frac{1}{n} \sum_{i=1}^{n} (y_i - \mu_i)^2}}
#' Where \eqn{y} is the vector of observed values and \eqn{\mu} is the vector of
#' predicted values.
#' 
#' @examples
#' y <- c(1, 2, 3)
#' mu <- c(1.1, 1.9, 3.2)
#' rmse(y, mu)
#' 
#' @export
rmse <- structure(function(y, mu, ...) {
  if (missing(y) || is.data.frame(y)) {
    rlang::check_installed("yardstick", reason = "for data frame metrics.")
    return(.fcr_yardstick(yardstick::rmse, if (!missing(y)) y,
                          if (!missing(mu)) rlang::enquo(mu),
                          rlang::enquos(...)))
  }
  sqrt(mean((y - mu)^2))
}, direction = "minimize", range = c(0, Inf),
class = c("numeric_metric", "metric", "function"))

#' Calculate Mean Absolute Error (MAE)
#'
#' This function calculates the Mean Absolute Error (MAE) between observed and
#' predicted values.
#'
#' @param y Numeric vector representing the observed values, or a data frame
#'   (see the section on yardstick).
#' @param mu Numeric vector representing the predicted values (with a data
#'   frame, the column of observed values).
#' @param ... with a data frame, further arguments of the yardstick metric
#'   (the column of predictions first).
#' @return Numeric value representing the MAE (with a data frame, the
#'   yardstick tibble).
#' @section yardstick: With a data frame, \code{mae(data, truth, estimate)}
#'   calls the \pkg{yardstick} metric of the same name, and the function can
#'   be used in \code{yardstick::metric_set()}.
#' @details The MAE is calculated using the formula:
#' \deqn{MAE = \frac{1}{n} \sum_{i=1}^{n} |y_i - \mu_i|}
#' Where \eqn{y} is the vector of observed values and \eqn{\mu} is the vector of
#' predicted values.
#' @examples
#' y <- c(1, 2, 3)
#' mu <- c(1.1, 1.9, 3.2)
#' mae(y, mu)
#' 
#' @export
mae <- structure(function(y, mu, ...) {
  if (missing(y) || is.data.frame(y)) {
    rlang::check_installed("yardstick", reason = "for data frame metrics.")
    return(.fcr_yardstick(yardstick::mae, if (!missing(y)) y,
                          if (!missing(mu)) rlang::enquo(mu),
                          rlang::enquos(...)))
  }
  mean(abs(y - mu))
}, direction = "minimize", range = c(0, Inf),
class = c("numeric_metric", "metric", "function"))

# Call a yardstick metric on the arguments of rmse() or mae() given a data
# frame (`data`) first, or only named arguments (as metric_set() does).
.fcr_yardstick <- function(metric, data, truth, dots) {
  args <- c(if (!is.null(data)) list(data = data),
            if (!is.null(truth)) list(truth), as.list(dots))
  rlang::eval_tidy(rlang::call2(metric, !!!args))
}
