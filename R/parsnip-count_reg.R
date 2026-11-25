#' Count regression model specification for tidymodels
#'
#' \code{count_reg()} defines a count regression model for the
#' \pkg{parsnip} package. Every estimator of flexCountReg is an engine, so
#' all of them can be used in workflows, resampling, tuning and workflow
#' sets. The count distribution (\code{family}) is a main argument.
#'
#' @param mode the model mode; only \code{"regression"} is available.
#' @param family the count distribution, one of the \code{family} values of
#'   \code{\link{countreg}} (for example \code{"NB2"}, \code{"PIG2"} or
#'   \code{"PLN"}) or a family from \code{\link{new_count_family}()}.
#'   \code{NULL} uses the default of the engine's estimator
#'   (\code{"NB2"} for \code{countreg()}). The engines \code{"renb"},
#'   \code{"repois"} and \code{"poislind_re"} have a fixed distribution and
#'   do not accept it; \code{"negmultinom"} and \code{"countreg_fe"} accept
#'   the families of their estimators; \code{"countreg_mv"} accepts one
#'   family per outcome.
#' @param engine the fitting engine; see the section \emph{Engines}.
#' @param object a \code{count_reg()} specification.
#' @param parameters a one-row tibble or named list of main arguments to
#'   update.
#' @param fresh logical; replace the arguments (\code{TRUE}) or modify them?
#' @param ... engine arguments to update.
#'
#' @section Engines:
#' \tabular{lll}{
#'   \strong{Engine} \tab \strong{Estimator} \tab \strong{Model} \cr
#'   \code{"flexCountReg"} \tab \code{\link{countreg}} \tab any family, with
#'     zero-inflated, hurdle and underreporting components \cr
#'   \code{"countreg_rp"} \tab \code{\link{countreg_rp}} \tab random
#'     parameters \cr
#'   \code{"countreg_mix"} \tab \code{\link{countreg_mix}} \tab finite
#'     mixtures \cr
#'   \code{"countreg_lc"} \tab \code{\link{countreg_lc}} \tab latent classes
#'     \cr
#'   \code{"countreg_fe"} \tab \code{\link{countreg_fe}} \tab fixed effects
#'     \cr
#'   \code{"countreg_hier"} \tab \code{\link{countreg_hier}} \tab
#'     hierarchical models \cr
#'   \code{"countreg_alloc"} \tab \code{\link{countreg_alloc}} \tab total
#'     plus allocation (several outcome columns) \cr
#'   \code{"countreg_mv"} \tab \code{\link{countreg_mv}} \tab multivariate
#'     (several outcomes) \cr
#'   \code{"renb"} \tab \code{\link{renb}} \tab random-effects negative
#'     binomial \cr
#'   \code{"repois"} \tab \code{\link{repois}} \tab random-effects Poisson
#'     \cr
#'   \code{"poislind_re"} \tab \code{\link{poislind_re}} \tab random-effects
#'     Poisson-Lindley \cr
#'   \code{"negmultinom"} \tab \code{\link{negmultinom}} \tab negative
#'     multinomial
#' }
#' All other arguments of the estimator are engine arguments, set with
#' \code{parsnip::set_engine()}, for example
#' \code{set_engine("countreg_rp", rpar_formula = ~ -1 + x, ndraws = 200)}.
#' Arguments with a finite set of choices have tuning parameters (see
#' \code{\link{count_family}}).
#'
#' The estimators' own \code{engine} argument (the numerical method) is
#' called \code{fit_engine} here, for example
#' \code{set_engine("countreg_hier", fit_engine = "laplace")}.
#'
#' The engine argument \code{predict_args} (a list) holds arguments passed
#' to the model's \code{predict()} method for numeric predictions and
#' intervals, for example \code{predict_args = list(method = "Exact")} for
#' \code{"countreg_rp"}, \code{list(newgroup = "mean_fe")} for
#' \code{"countreg_fe"} or \code{list(marginal = TRUE)} for
#' \code{"countreg_hier"}.
#'
#' @section Formulas, recipes and extra columns:
#' The engines use the formula interface. In a workflow, the model receives
#' only the columns the preprocessor declares as outcomes and predictors
#' (not those with other roles), so every column the model needs must be a
#' predictor: offsets, grouping variables, and variables of auxiliary
#' formulas (\code{dis_param_formula_1}, \code{rpar_formula}, ...). Two ways
#' to say which of them enter the mean model:
#' \itemize{
#'   \item Give the model formula, as in
#'   \code{workflows::add_model(spec, formula = y ~ x + offset(log_exp))},
#'   with \code{workflows::add_variables()} or a recipe as the preprocessor.
#'   This works for every model, including smooths (\code{s(x)}), formula
#'   offsets, random effects (\code{(1 | site)}), \code{tv()} terms and
#'   several outcomes (\code{cbind(y1, y2) ~ x}).
#'   \item Without a model formula, all predictors enter the mean model
#'   (\code{y ~ .}), except that columns named in \code{offset},
#'   \code{group_var}, \code{panel_id}, \code{cluster}, \code{fixef},
#'   \code{cre_group} or \code{time}, and the variables of
#'   \code{rpar_formula} for \code{"countreg_rp"} (which are random instead
#'   of fixed), are left out of the \code{.}. The outcome is then named
#'   \code{..y} in the fitted model.
#' }
#' With a recipe, do not transform the grouping columns (for example with
#' \code{step_dummy()} or \code{step_normalize()}) or an offset.
#'
#' @section Case weights:
#' Case weights are passed to the estimator's \code{weights}. By default
#' they are frequency weights. parsnip passes the engine only their values,
#' so for importance weights with \code{"flexCountReg"} or
#' \code{"countreg_rp"} set the engine argument
#' \code{weight_type = "sampling"}; a warning is given for weights that are
#' not whole numbers without it.
#'
#' @section Predictions:
#' \code{type = "numeric"} gives the expected counts of the estimator's
#' \code{predict()} method, with one row per row of \code{new_data} and
#' \code{NA} for rows with missing values. For models with several outcomes
#' there is one column per outcome or category
#' (\code{.pred_<name>}); \code{predict_args = list(type = "total")} gives
#' the total of \code{"countreg_alloc"} as \code{.pred}. All
#' engines except \code{"countreg_fe"} also have \code{type = "conf_int"} and
#' \code{type = "pred_int"} (see \code{\link{predict.flexCountReg}}), with a
#' pair of columns per outcome for models with several
#' (\code{.pred_lower_<name>}, \code{.pred_upper_<name>}).
#' Everything else is available with \code{type = "raw"}, whose \code{opts}
#' are passed to \code{predict()}: for example
#' \code{predict(fit, new_data, type = "raw",
#' opts = list(type = "prob", counts = 0:5))} for count probabilities.
#'
#' Predictions that use the observed outcome (Empirical Bayes predictions of
#' random-parameter models, posterior class probabilities) need it in
#' \code{new_data}; a workflow removes the outcome, so use the model itself:
#' \code{predict(hardhat::extract_fit_engine(wf), newdata = ...)}. All other
#' functions of the package (\code{\link{margEffTable}},
#' \code{\link{cureplot}}, \code{\link{rootogram}}, \code{vcov()}, ...) also
#' take the extracted model.
#'
#' @section Several outcomes:
#' The engines \code{"countreg_alloc"} and \code{"countreg_mv"} can be
#' fitted and predicted in workflows, but the resampling functions of
#' \pkg{tune} compute metrics for one outcome, so use
#' \code{\link{fit_resamples_counts}()} to validate them.
#'
#' @section Loading:
#' Requires the \pkg{parsnip} package. The engines are registered when
#' flexCountReg is loaded; load it before using a saved model or workflow.
#' Parallel workers of \pkg{tune} load it themselves.
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
#'   predict(fit, new_data = washington_roads[1:5, ], type = "pred_int")
#'
#'   # Another estimator
#'   re <- count_reg() |>
#'     parsnip::set_engine("repois", group_var = "ID", method = "BFGS")
#'   fit_re <- parsnip::fit(re, Total_crashes ~ lnaadt + lnlength,
#'                          data = washington_roads)
#'   tidy(fit_re)
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

#' @rdname count_reg
#' @method update count_reg
#' @importFrom stats update
#' @export
update.count_reg <- function(object, parameters = NULL, family = NULL,
                             fresh = FALSE, ...) {
  args <- list(family = rlang::enquo(family))
  parsnip::update_spec(object = object, parameters = parameters,
                       args_enquo_list = args, fresh = fresh,
                       cls = "count_reg", ...)
}

#' Tuning parameters for count regression models
#'
#' \pkg{dials} parameters for \code{\link{count_reg}}: the count
#' distribution (the main argument \code{family}) and engine arguments.
#' Arguments of \code{count_reg()} marked with \code{tune()} use them
#' automatically.
#'
#' @param values the candidate values.
#' @param range the smallest and largest number of components or classes.
#' @param trans a transformation (see \code{dials::new_quant_param()}).
#' @details
#' \tabular{lll}{
#'   \strong{Parameter} \tab \strong{Argument} \tab \strong{Engine} \cr
#'   \code{count_family()} \tab \code{family} \tab all with a family \cr
#'   \code{count_components()} \tab \code{components}, \code{classes} \tab
#'     \code{"countreg_mix"}, \code{"countreg_lc"} \cr
#'   \code{count_zero_type()} \tab \code{zero_type} \tab
#'     \code{"flexCountReg"} \cr
#'   \code{count_correlated()} \tab \code{correlated} \tab
#'     \code{"countreg_rp"}, \code{"countreg_mix"}, \code{"countreg_lc"} \cr
#'   \code{count_re_dist()} \tab \code{re_dist} \tab
#'     \code{"countreg_hier"} (one term) \cr
#' }
#' The families supported by \code{"countreg_fe"} (\code{"POISSON"} and
#' \code{"NB2"} by default here) and \code{"negmultinom"} (\code{"NM1"},
#' \code{"NM2"} and \code{"NMP"}) are their own default values. Any other engine argument
#' can be tuned with an explicit grid, or with \code{update()} on the
#' parameter set.
#' @returns A \pkg{dials} parameter.
#' @examples
#' if (requireNamespace("dials", quietly = TRUE)) {
#'   count_family()
#'   count_components()
#' }
#' @name count_family
NULL

#' @rdname count_family
#' @export
count_family <- function(values = c("POISSON", "NB1", "NB2", "NBP", "GP1",
                                    "GP2", "PIG1", "PIG2", "PL", "PLN",
                                    "COM")) {
  .fcr_qual_param("character", values, c(family = "Count distribution"))
}

#' @rdname count_family
#' @export
count_components <- function(range = c(2L, 4L), trans = NULL) {
  rlang::check_installed("dials", reason = "to use `count_components()`.")
  dials::new_quant_param(type = "integer", range = range,
                         inclusive = c(TRUE, TRUE), trans = trans,
                         label = c(components = "# Mixture Components"))
}

#' @rdname count_family
#' @export
count_zero_type <- function(values = c("inflated", "hurdle")) {
  .fcr_qual_param("character", values, c(zero_type = "Zero Model"))
}

#' @rdname count_family
#' @export
count_correlated <- function(values = c(FALSE, TRUE)) {
  .fcr_qual_param("logical", values,
                  c(correlated = "Correlated Random Parameters"))
}

#' @rdname count_family
#' @export
count_re_dist <- function(values = c("normal", "uniform", "triangular",
                                     "gamma")) {
  .fcr_qual_param("character", values,
                  c(re_dist = "Random Effect Distribution"))
}

.fcr_qual_param <- function(type, values, label) {
  rlang::check_installed("dials", reason = "for tuning parameters.")
  dials::new_qual_param(type = type, values = values, label = label)
}

# tune parameters of the engine arguments (the registered main arguments
# come from parsnip's method)
.fcr_tunable_count_reg <- function(x, ...) {
  res <- NextMethod()
  if (is.null(x$engine)) return(res)
  eng <- .fcr_parsnip_engines()[[x$engine]]
  dial <- function(fun, ...) {
    list(list(pkg = "flexCountReg", fun = fun, ...))
  }
  for (nm in names(eng$tunable)) {
    i <- which(res$name == nm)
    if (length(i)) res$call_info[i] <- dial(eng$tunable[[nm]])
  }
  if (!is.null(eng$families)) {
    res$call_info[res$name == "family"] <- dial("count_family",
                                                values = eng$families)
  }
  res
}

# The parsnip engines: the estimator of each, its tunable engine arguments,
# its families (when not those of count_family()) and whether its predict()
# has intervals.
.fcr_parsnip_engines <- function() {
  mix <- c(correlated = "count_correlated")
  list(
    flexCountReg = list(fun = "countreg",
                        tunable = c(zero_type = "count_zero_type")),
    countreg_rp = list(fun = "countreg_rp",
                       tunable = c(correlated = "count_correlated")),
    countreg_mix = list(fun = "countreg_mix",
                        tunable = c(components = "count_components", mix)),
    countreg_lc = list(fun = "countreg_lc",
                       tunable = c(classes = "count_components", mix)),
    countreg_fe = list(fun = "countreg_fe", intervals = FALSE,
                       families = c("POISSON", "NB2")),
    countreg_hier = list(fun = "countreg_hier",
                         tunable = c(re_dist = "count_re_dist")),
    countreg_alloc = list(fun = "countreg_alloc"),
    countreg_mv = list(fun = "countreg_mv"),
    renb = list(fun = "renb"),
    repois = list(fun = "repois"),
    poislind_re = list(fun = "poislind_re"),
    negmultinom = list(fun = "negmultinom", families = c("NM1", "NM2", "NMP")))
}

#' Fit and predict wrappers used by the parsnip engines
#'
#' Internal functions called by \pkg{parsnip}; use \code{\link{count_reg}}
#' or the estimators instead.
#'
#' @param formula,data,family,weights,... passed to the estimator.
#'   \code{family = NULL} uses the estimator's default.
#' @param model the name of the estimator.
#' @param predict_args a list of arguments for the model's \code{predict()}
#'   method.
#' @param fit_engine the estimator's \code{engine} argument.
#' @param object a fitted model.
#' @param new_data data to predict.
#' @param format \code{"numeric"}, \code{"interval"} or \code{"raw"}.
#' @returns \code{count_reg_fit()}: the estimator's model.
#'   \code{count_reg_predict()}: its predictions.
#' @keywords internal
#' @export
count_reg_fit <- function(formula, data, family = NULL, weights = NULL, ...,
                          model = "countreg", predict_args = NULL,
                          fit_engine = NULL) {
  fun <- get(model, envir = topenv(), mode = "function")
  fun_args <- names(formals(fun))
  dots <- list(...)
  # set_engine() has an argument called `engine` itself
  if (!is.null(fit_engine)) dots$engine <- fit_engine
  bad <- setdiff(names(dots), fun_args)
  if (length(bad)) {
    stop("Engine argument(s) not used by ", model, "(): ", toString(bad),
         ".", call. = FALSE)
  }
  if (!is.null(predict_args) && !is.list(predict_args)) {
    stop("`predict_args` must be a list.", call. = FALSE)
  }
  if (!is.null(family) && !"family" %in% fun_args) {
    stop(model, "() has a fixed distribution: do not set `family` with ",
         "this engine.", call. = FALSE)
  }
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (is.null(dots$weight_type) &&
        any(abs(weights - round(weights)) > 1e-8, na.rm = TRUE)) {
      warning("Case weights that are not whole numbers are used as if they ",
              "were frequencies, so standard errors, likelihoods and BIC ",
              "are those of the expanded data",
              if ("weight_type" %in% fun_args)
                "; for importance or sampling weights set the engine argument weight_type = \"sampling\"",
              ".", call. = FALSE)
    }
  }
  formula <- .fcr_parsnip_formula(formula, data, dots, model)
  # The call holds the values of the arguments, so the stored call can be
  # evaluated again (bootstrap refits, update()) wherever it is
  fit_call <- as.call(c(
    list(call("::", quote(flexCountReg), as.name(model)),
         formula = formula, data = quote(data)),
    if (!is.null(family)) list(family = family),
    if (!is.null(weights)) list(weights = quote(weights)), dots))
  out <- eval(fit_call)
  out$parsnip <- list(predict_args = predict_args, weights = weights)
  out
}

# Replace `.` on the right side of a formula by the columns of `data`,
# without the outcome and the columns that have a role of their own.
.fcr_parsnip_formula <- function(formula, data, args, model) {
  if (!inherits(formula, "formula") || length(formula) != 3L ||
      !"." %in% all.vars(formula[[3L]])) {
    return(formula)
  }
  role <- function(x) {
    if (inherits(x, "formula")) all.vars(x) else
      if (is.character(x)) x else character(0)
  }
  drop <- unlist(lapply(args[intersect(names(args), c(
    "offset", "group_var", "panel_id", "cluster", "fixef", "cre_group",
    "time"))], role))
  if (model == "countreg_rp") drop <- c(drop, role(args$rpar_formula))
  rhs <- formula[[3L]]
  # Variables written in the formula itself with `.` (such as the group of
  # (1 | site)) are repeated only if they are terms of their own, as in lm()
  bars <- unlist(lapply(.fcr_find_bars(rhs), function(b) all.vars(b[[3L]])))
  keep <- setdiff(names(data), c(all.vars(formula[[2L]]), drop, bars))
  dot <- if (!length(keep)) 1 else
    Reduce(function(a, b) call("+", a, b), lapply(keep, as.name))
  if (length(keep) > 1L && !identical(rhs, quote(.))) dot <- call("(", dot)
  formula[[3L]] <- do.call(substitute, list(rhs, list(. = dot)))
  formula
}

# The `a | b` calls in an expression.
.fcr_find_bars <- function(e) {
  if (!is.call(e)) return(NULL)
  if (identical(e[[1L]], as.name("|")) || identical(e[[1L]], as.name("||"))) {
    return(list(e))
  }
  unlist(lapply(as.list(e)[-1L], .fcr_find_bars), recursive = FALSE)
}

#' @rdname count_reg_fit
#' @export
count_reg_predict <- function(object, new_data, ...,
                              format = c("numeric", "interval", "raw")) {
  format <- match.arg(format)
  args <- list(...)
  extra <- object$parsnip$predict_args
  if (format != "raw") {
    args <- c(args, extra[setdiff(names(extra), names(args))])
  }
  if (format != "raw" && !nrow(new_data)) {
    if (format == "interval") {
      nm <- if (identical(args$type, "total")) NULL else
        object$outcomes %||% object$categories
      nm <- if (is.null(nm)) c(".pred_lower", ".pred_upper") else
        paste0(c(".pred_lower_", ".pred_upper_"), rep(nm, each = 2L))
      return(tibble::as_tibble(stats::setNames(
        rep(list(numeric()), length(nm)), nm)))
    }
    nm <- if (!identical(args$type %||% "counts", "counts")) NULL else
      object$outcomes %||% object$categories
    if (is.null(nm)) return(numeric(0))
    return(tibble::as_tibble(stats::setNames(
      rep(list(numeric()), length(nm)), paste0(".pred_", nm))))
  }
  p <- do.call(stats::predict,
               c(list(object, newdata = as.data.frame(new_data)), args))
  if (format == "raw") return(p)
  if (format == "interval") {
    if (is.data.frame(p)) {
      return(tibble::tibble(.pred_lower = p$lwr, .pred_upper = p$upr))
    }
    # Several outcomes: a pair of columns for each
    p <- unlist(lapply(names(p), function(nm) stats::setNames(
      list(p[[nm]]$lwr, p[[nm]]$upr),
      paste0(c(".pred_lower_", ".pred_upper_"), nm))), recursive = FALSE)
    return(tibble::as_tibble(p))
  }
  if (is.matrix(p)) {
    nm <- colnames(p) %||% seq_len(ncol(p))
    p <- tibble::as_tibble(as.data.frame(p), .name_repair = "minimal")
    names(p) <- paste0(".pred_", nm)
    return(p)
  }
  as.vector(p)
}

# parsnip's augment() for models with several outcomes: parsnip's own method
# expects one column of predictions.
.fcr_augment_multi <- function(x, new_data, ...) {
  nd <- tibble::as_tibble(new_data)
  p <- stats::predict(x, new_data = new_data)
  res <- lapply(sub("^\\.pred_", "", names(p)), function(v) {
    if (v %in% names(nd)) nd[[v]] - p[[paste0(".pred_", v)]]
  })
  names(res) <- sub("^\\.pred_", ".resid_", names(p))
  tibble::as_tibble(c(as.list(p), Filter(Negate(is.null), res), as.list(nd)),
                    .name_repair = "unique_quiet")
}

# Refit a model fitted by count_reg_fit() (stored call with `data` and
# `weights` as symbols) with its own data
.fcr_update_parsnip <- function(object, call, env) {
  if (identical(call$data, quote(data))) {
    d <- object[["data"]] %||% object$data_used
    if (is.null(d)) {
      stop("The model was fit with keep_data = FALSE; give `data` to ",
           "update().", call. = FALSE)
    }
    # Group-mean columns are added again by the refit
    d <- as.data.frame(d)[setdiff(names(d), object$model$cre$cols)]
    # (countreg_alloc() and countreg_mv() store only the estimation rows)
    w <- object$parsnip$weights
    if (is.null(object[["data"]]) && length(w) != nrow(d)) {
      w <- w[object$rows_used]
    }
    env <- list2env(list(data = d, weights = w), parent = env)
  }
  out <- eval(call, env)
  out$parsnip <- object$parsnip
  # With other data the stored weights no longer apply
  out$parsnip["weights"] <- list(if (identical(call$weights, quote(weights)))
    get("weights", env))
  out
}

#' Update and refit a flexCountReg model
#'
#' As \code{\link[stats]{update}}, also for models fitted with tidymodels
#' (see \code{\link{count_reg}}), which are refitted with their stored
#' training data.
#'
#' @param object a fitted model.
#' @param formula. changes to the formula (see
#'   \code{\link[stats]{update.formula}}).
#' @param ... arguments to change or add.
#' @param evaluate if \code{FALSE}, return the call instead of the model.
#' @returns The refitted model, or the call.
#' @method update flexCountReg
#' @export
update.flexCountReg <- function(object, formula., ..., evaluate = TRUE) {
  # As update.default(), with its own match.call(): forwarding `...` to it
  # would store the arguments as ..1, ..2 in the call
  call <- object$call
  if (is.null(call)) {
    stop("The model has no stored call, so it cannot be updated.",
         call. = FALSE)
  }
  extras <- match.call(expand.dots = FALSE)$...
  from_parsnip <- !is.null(object$parsnip)
  if (!missing(formula.)) {
    old <- if (from_parsnip) call$formula else stats::formula(object)
    call$formula <- stats::update(old, formula.)
  }
  # An argument set to NULL is removed from the call
  for (nm in names(extras)) call[[nm]] <- extras[[nm]]
  if (from_parsnip && "data" %in% names(extras) &&
      !"weights" %in% names(extras)) {
    call$weights <- NULL
  }
  if (!evaluate) return(call)
  if (from_parsnip) {
    return(.fcr_update_parsnip(object, call, parent.frame()))
  }
  eval(call, parent.frame())
}

#' @rdname update.flexCountReg
#' @method update countreg_alloc
#' @export
update.countreg_alloc <- update.flexCountReg

#' @rdname update.flexCountReg
#' @method update countreg_mv
#' @export
update.countreg_mv <- update.flexCountReg

# Register the count_reg model and its engines with parsnip.
.fcr_register_count_reg <- function() {
  if (!requireNamespace("parsnip", quietly = TRUE)) return(invisible(FALSE))
  current <- parsnip::get_model_env()
  if (!"count_reg" %in% current$models) {
    parsnip::set_new_model("count_reg")
    parsnip::set_model_mode("count_reg", "regression")
  }
  # Idempotent: register each engine unless it is already registered (another
  # package may have registered a count_reg model without these engines)
  done <- tryCatch(parsnip::get_from_env("count_reg")$engine,
                   error = function(e) NULL)
  engines <- .fcr_parsnip_engines()
  for (eng in setdiff(names(engines), done)) {
    .fcr_register_engine(eng, engines[[eng]])
  }
  registerS3method("tunable", "count_reg", .fcr_tunable_count_reg,
                   envir = asNamespace("generics"))
  for (cl in c("_countreg_alloc", "_countreg_mv")) {
    registerS3method("augment", cl, .fcr_augment_multi,
                     envir = asNamespace("generics"))
  }
  invisible(TRUE)
}

.fcr_register_engine <- function(eng, info) {
  parsnip::set_model_engine("count_reg", "regression", eng)
  parsnip::set_dependency("count_reg", eng, "flexCountReg",
                          mode = "regression")
  # Registered for every engine: an engine with a fixed distribution stops
  # with a message instead of ignoring `family`
  parsnip::set_model_arg(
    model = "count_reg", eng = eng, parsnip = "family",
    original = "family",
    func = list(pkg = "flexCountReg", fun = "count_family"),
    has_submodel = FALSE)
  parsnip::set_fit(
    model = "count_reg", eng = eng, mode = "regression",
    value = list(interface = "formula",
                 protect = c("formula", "data", "weights"),
                 func = c(pkg = "flexCountReg", fun = "count_reg_fit"),
                 defaults = if (info$fun == "countreg") list() else
                   list(model = info$fun)))
  parsnip::set_encoding(
    model = "count_reg", eng = eng, mode = "regression",
    options = list(predictor_indicators = "none", compute_intercept = FALSE,
                   remove_intercept = FALSE, allow_sparse_x = FALSE))
  pred <- function(type, args) {
    parsnip::set_pred(
      model = "count_reg", eng = eng, mode = "regression", type = type,
      value = list(
        pre = NULL, post = NULL,
        func = c(pkg = "flexCountReg", fun = "count_reg_predict"),
        args = c(list(object = rlang::expr(object$fit),
                      new_data = rlang::expr(new_data)), args)))
  }
  pred("numeric", list(format = "numeric"))
  pred("raw", list(format = "raw"))
  if (!isFALSE(info$intervals)) {
    pred("conf_int", list(format = "interval", interval = "confidence",
                          level = quote(level)))
    pred("pred_int", list(format = "interval", interval = "prediction",
                          level = quote(level)))
  }
}
