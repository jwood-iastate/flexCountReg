# Compile an analysis specification into R language objects. The same
# expressions are evaluated by the background worker and deparsed into the
# exported script, so the code shown and downloaded is the code that ran.
# Column names enter only as symbols or strings (never pasted into text that
# is parsed), so any name, including spaces, quotes or backticks, is safe.

.sym <- function(x) as.name(x)

# a + b + c from a list of language objects (1 when empty)
.plus <- function(terms) {
  if (!length(terms)) return(1)
  Reduce(function(a, b) call("+", a, b), terms)
}

# Formula as a call, e.g. y ~ x + offset(log(e)); one-sided when lhs is NULL
.formula_call <- function(lhs, rhs_terms) {
  rhs <- .plus(rhs_terms)
  if (is.null(lhs)) call("~", rhs) else call("~", lhs, rhs)
}

# Expressions that load and prepare the data (assigning `analysis_data`)
compile_prep <- function(spec) {
  d <- spec$data
  r <- spec$recipe
  out <- list()
  if (identical(d$source, "bundled")) {
    out$load <- bquote(analysis_data <- local({
      env <- new.env()
      utils::data(list = .(d$name), package = "flexCountReg", envir = env)
      env[[.(d$name)]]
    }))
  } else {
    out$load <- quote(analysis_data <- utils::read.csv(
      data_file, check.names = FALSE, stringsAsFactors = FALSE))
  }
  for (f in r$factors) {
    out[[paste0("factor_", f)]] <- bquote(
      analysis_data[[.(f)]] <- factor(analysis_data[[.(f)]]))
  }
  out
}

# The model formula (as a call)
compile_formula <- function(spec) {
  r <- spec$recipe
  terms <- lapply(r$predictors, .sym)
  if (!is.null(r$exposure)) {
    terms <- c(terms, list(call("offset", call("log", .sym(r$exposure)))))
  }
  if (!is.null(r$log_offset)) {
    terms <- c(terms, list(call("offset", .sym(r$log_offset))))
  }
  .formula_call(.sym(r$response), terms)
}

# The fitting call, assigned to `fit`
compile_fit <- function(spec) {
  r <- spec$recipe
  m <- spec$model
  args <- list(formula = compile_formula(spec), data = quote(analysis_data))
  fam <- unname(m$family)
  args$family <- if (length(fam) == 1) fam else as.call(c(.sym("c"),
                                                            as.list(fam)))
  if (length(m$dis1)) {
    args$dis_param_formula_1 <- .formula_call(NULL, lapply(m$dis1, .sym))
  }
  if (!is.null(r$weights)) args$weights <- r$weights
  if (m$estimator == "countreg") {
    fun <- quote(flexCountReg::countreg)
  } else {
    lc <- m$estimator == "countreg_lc"
    fun <- if (lc) quote(flexCountReg::countreg_lc) else
      quote(flexCountReg::countreg_mix)
    if (length(fam) == 1) args[[if (lc) "classes" else "components"]] <-
      as.integer(m$components)
    if (length(m$mixing)) {
      args[[if (lc) "membership_formula" else "mixing_formula"]] <-
        .formula_call(NULL, lapply(m$mixing, .sym))
    }
    if (length(m$shared_mean) || isTRUE(m$shared_dis1)) {
      sh <- list()
      if (length(m$shared_mean)) {
        sh$mean <- if (length(m$shared_mean) == 1) m$shared_mean else
          as.call(c(.sym("c"), as.list(m$shared_mean)))
      }
      if (isTRUE(m$shared_dis1)) sh$dis1 <- TRUE
      args$shared <- as.call(c(.sym("list"), sh))
    }
    if (!is.null(r$group)) args$group_var <- r$group
    args$n_starts <- as.integer(m$n_starts)
    args$seed <- as.integer(m$seed)
  }
  args$method <- m$method
  args$max.iters <- as.integer(m$max_iters)
  if (!identical(m$stderr, "normal")) args$stderr <- m$stderr
  call("<-", .sym("fit"), as.call(c(fun, args)))
}

# A data.frame() call that rebuilds a small table (scenarios)
compile_data_frame <- function(df) {
  cols <- lapply(df, function(v) {
    if (is.factor(v)) {
      call("factor", as.call(c(.sym("c"), as.list(as.character(v)))),
           levels = as.call(c(.sym("c"), as.list(levels(v)))))
    } else {
      as.call(c(.sym("c"), as.list(v)))
    }
  })
  as.call(c(.sym("data.frame"), cols, check.names = FALSE))
}

# Post-estimation expressions requested in the app (diagnostics and
# predictions), each producing a named object
compile_post <- function(spec, options = list()) {
  out <- list(summary = quote(fit_summary <- summary(fit)))
  if (isTRUE(options$rootogram)) {
    out$rootogram <- quote(rootogram_data <- flexCountReg::rootogram(
      fit, overflow = TRUE))
  }
  if (isTRUE(options$quantile_residuals)) {
    out$residuals <- quote(quantile_residuals <- {
      set.seed(1)
      stats::residuals(fit, type = "quantile")
    })
  }
  if (!is.null(options$cure_var)) {
    out$cure <- bquote(cure_plot <- flexCountReg::cureplot(
      fit, indvar = .(options$cure_var)))
  }
  if (!is.null(options$scenarios) && nrow(options$scenarios)) {
    out$scenarios <- call("<-", .sym("scenarios"),
                          compile_data_frame(options$scenarios))
    out$predictions <- quote(scenario_predictions <- cbind(
      scenarios, expected_count = stats::predict(fit, newdata = scenarios)))
  }
  out
}

# Compile an analysis specification
#
compile_analysis <- function(spec, options = list()) {
  list(prep = compile_prep(spec), fit = list(fit = compile_fit(spec)),
       post = compile_post(spec, options))
}

.deparse <- function(e) paste(deparse(e, width.cutoff = 72L), collapse = "\n")

# Reproducible R script of an analysis
#
analysis_script <- function(spec, options = list(), fit_label = NULL,
                            fit_id = NULL, versions = fcr_versions()) {
  a <- compile_analysis(spec, options)
  hdr <- c(
    "# Reproducible analysis generated by the flexCountReg app.",
    "# It runs without the app; only flexCountReg (and its dependencies)",
    "# are needed.",
    if (!is.null(fit_label)) paste0("# Fit: ", fit_label,
                                    if (!is.null(fit_id))
                                      paste0(" (request ", fit_id, ")")),
    paste0("# Versions: ", paste(names(versions), versions, sep = " ",
                                 collapse = ", ")),
    paste0("# Data: ", if (spec$data$source == "bundled")
      paste0("bundled dataset ", spec$data$name, " from flexCountReg") else
        paste0("uploaded file ", spec$data$name)),
    "")
  body <- c("library(flexCountReg)", "")
  if (spec$data$source != "bundled") {
    body <- c(body, "# Path to your copy of the data file:",
              paste0("data_file <- ", .deparse(spec$data$name)), "")
  }
  sec <- function(title, exprs) {
    if (!length(exprs)) return(character(0))
    c(paste("#", title), vapply(exprs, .deparse, character(1)), "")
  }
  body <- c(body, sec("Data", a$prep), sec("Model", a$fit))
  post <- a$post
  body <- c(body, sec("Results", post["summary"]), "fit_summary", "")
  rest <- post[setdiff(names(post), "summary")]
  body <- c(body, sec("Diagnostics and predictions", rest))
  if ("predictions" %in% names(rest)) body <- c(body, "scenario_predictions")
  paste(c(hdr, body), collapse = "\n")
}

# Evaluate compiled expressions in a fresh environment (used by the worker
# and by tests). `data_file` is the path of an uploaded data file.
run_compiled <- function(exprs, data_file = NULL, envir = NULL) {
  env <- envir %||% new.env(parent = globalenv())
  if (!is.null(data_file)) assign("data_file", data_file, envir = env)
  for (e in exprs) eval(e, env)
  env
}
