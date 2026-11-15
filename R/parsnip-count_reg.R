#' Count regression model specification for tidymodels
#'
#' \code{count_reg()} defines a count regression model for the
#' \pkg{parsnip} package, fitted with the \code{"flexCountReg"} engine
#' (\code{\link{countreg}}). The count distribution is a main argument, so it
#' can be tuned with \pkg{tune}.
#'
#' @param mode the model mode; only \code{"regression"} is available.
#' @param family the count distribution, one of the \code{family} values of
#'   \code{\link{countreg}} (for example \code{"NB2"}, \code{"PIG2"} or
#'   \code{"PLN"}). Default \code{"NB2"}.
#' @param engine the fitting engine; only \code{"flexCountReg"}.
#'
#' @details
#' Other \code{countreg()} arguments are engine arguments, set with
#' \code{parsnip::set_engine("flexCountReg", ...)}: for example
#' \code{method}, \code{dis_param_formula_1}, \code{underreport_formula},
#' \code{offset} or \code{stderr}. Variables used only in those auxiliary
#' formulas, and offset columns, must be kept in the data by any recipe (in a
#' workflow, add them with \code{workflows::add_variables()} or give them a
#' non-predictor role). A formula \code{offset()} term can be used with
#' \code{workflows::add_model(spec, formula = y ~ x + offset(log_exposure))}.
#' Case weights from tidymodels are passed to \code{countreg()} as frequency
#' weights.
#'
#' Predictions (\code{type = "numeric"}) are the expected counts from
#' \code{\link{predict.flexCountReg}} and have one row per row of
#' \code{new_data}, with \code{NA} for rows with missing values.
#'
#' Requires the \pkg{parsnip} package. The engine is registered when
#' flexCountReg or parsnip is loaded.
#'
#' @returns A \code{count_reg} model specification.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("parsnip", quietly = TRUE)) {
#'   spec <- count_reg(family = "NB2") |>
#'     parsnip::set_engine("flexCountReg", method = "BFGS")
#'   fit <- parsnip::fit(spec, Total_crashes ~ lnaadt + lnlength,
#'                       data = washington_roads)
#'   predict(fit, new_data = washington_roads[1:5, ])
#' }
#' }
#' @export
count_reg <- function(mode = "regression", family = NULL,
                      engine = "flexCountReg") {
  rlang::check_installed("parsnip", reason = "to use `count_reg()`.")
  .fcr_register_count_reg()
  args <- list(family = rlang::enquo(family))
  parsnip::new_model_spec("count_reg", args = args, eng_args = NULL,
                          mode = mode, user_specified_mode = !missing(mode),
                          method = NULL, engine = engine,
                          user_specified_engine = !missing(engine))
}

#' Count distribution tuning parameter
#'
#' A \pkg{dials} parameter for the \code{family} argument of
#' \code{\link{count_reg}}, for tuning over count distributions.
#'
#' @param values the candidate families.
#' @returns A \pkg{dials} qualitative parameter.
#' @examples
#' if (requireNamespace("dials", quietly = TRUE)) count_family()
#' @export
count_family <- function(values = c("POISSON", "NB1", "NB2", "NBP", "GP1",
                                    "GP2", "PIG1", "PIG2", "PL", "PLN",
                                    "COM")) {
  rlang::check_installed("dials", reason = "to use `count_family()`.")
  dials::new_qual_param(type = "character", values = values,
                        label = c(family = "Count distribution"))
}

#' Fit wrapper used by the parsnip engine
#'
#' Internal function called by \pkg{parsnip}; use \code{\link{count_reg}}
#' or \code{\link{countreg}} instead.
#'
#' @param formula,data,family,weights,... passed to \code{\link{countreg}}.
#' @returns A \code{countreg()} model.
#' @keywords internal
#' @export
count_reg_fit <- function(formula, data, family = "NB2", weights = NULL, ...) {
  if (!is.null(weights)) weights <- as.numeric(weights)
  countreg(formula, data = data, family = family, weights = weights, ...)
}

# Register the count_reg model and its flexCountReg engine with parsnip.
.fcr_register_count_reg <- function() {
  if (!requireNamespace("parsnip", quietly = TRUE)) return(invisible(FALSE))
  current <- parsnip::get_model_env()
  if ("count_reg" %in% current$models) return(invisible(TRUE))

  parsnip::set_new_model("count_reg")
  parsnip::set_model_mode("count_reg", "regression")
  parsnip::set_model_engine("count_reg", "regression", "flexCountReg")
  parsnip::set_dependency("count_reg", "flexCountReg", "flexCountReg",
                          mode = "regression")
  parsnip::set_model_arg(
    model = "count_reg", eng = "flexCountReg", parsnip = "family",
    original = "family",
    func = list(pkg = "flexCountReg", fun = "count_family"),
    has_submodel = FALSE)
  parsnip::set_fit(
    model = "count_reg", eng = "flexCountReg", mode = "regression",
    value = list(interface = "formula",
                 protect = c("formula", "data", "weights"),
                 func = c(pkg = "flexCountReg", fun = "count_reg_fit"),
                 defaults = list()))
  parsnip::set_encoding(
    model = "count_reg", eng = "flexCountReg", mode = "regression",
    options = list(predictor_indicators = "none", compute_intercept = FALSE,
                   remove_intercept = FALSE, allow_sparse_x = FALSE))
  parsnip::set_pred(
    model = "count_reg", eng = "flexCountReg", mode = "regression",
    type = "numeric",
    value = list(pre = NULL, post = NULL, func = c(fun = "predict"),
                 args = list(object = rlang::expr(object$fit),
                             new_data = rlang::expr(new_data),
                             type = "response")))
  invisible(TRUE)
}
