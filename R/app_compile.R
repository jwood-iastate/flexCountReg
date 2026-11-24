# Compile what the app runs into R language objects. Every action (loading
# data, a fit, a diagnostic, a prediction, a tolerance interval, a
# comparison) becomes a list of expressions that the app evaluates and
# appends to the session script, so the script is exactly the code that ran.
# Column names enter only as symbols or strings (never as pasted text that
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

.cvec <- function(x) {
  if (length(x) == 1) x else as.call(c(.sym("c"), as.list(x)))
}

# ---- Data -------------------------------------------------------------------------

# Path of an uploaded file in the exported bundle, relative to the script:
# data/<object>.csv. It comes from the object name (data_1, newdata_1, ...),
# which is unique in a session, and not from the name of the uploaded file,
# so two uploads with the same file name never collapse into one file.
.app_upload_file <- function(object) paste0("data/", object, ".csv")

# Expressions that load a dataset into `object`. For an uploaded file the
# script names the file by its path in the bundle
# (`<object>_file <- "data/<object>.csv"`, see .app_upload_file()); the app
# instead points that variable to its own copy of the upload. Returns
# list(display, run): the expressions shown in the script and those the app
# evaluates.
compile_data_load <- function(source, name, object,
                              file = .app_upload_file(object)) {
  obj <- .sym(object)
  if (identical(source, "bundled")) {
    e <- bquote(.(obj) <- local({
      env <- new.env()
      utils::data(list = .(name), package = "flexCountReg", envir = env)
      env[[.(name)]]
    }))
    return(list(display = list(e), run = list(e)))
  }
  fvar <- .sym(paste0(object, "_file"))
  read <- bquote(.(obj) <- utils::read.csv(.(fvar), check.names = FALSE,
                                           stringsAsFactors = FALSE))
  list(display = list(bquote(.(fvar) <- .(file)), read), run = list(read))
}

# ---- Models -----------------------------------------------------------------------

# The model formula (as a call). Random parameters of countreg_rp() are
# left out of the fixed formula (they are in rpar_formula). Interactions
# (pairs of predictors) follow the main effects as a:b; with countreg_rp an
# interaction that involves a random parameter stays a fixed term of this
# formula (the variable's own effect is the random parameter). Smooth
# predictors enter as s(x). A hierarchical model adds its random-effect
# terms in lme4's notation: (1 + slopes | group), with || for uncorrelated
# effects, and (1 | group2) (crossed) or (1 | group:group2) (nested).
compile_formula <- function(spec) {
  r <- spec$recipe
  m <- spec$model
  preds <- r$predictors
  rp <- m$estimator == "countreg_rp"
  if (rp) preds <- setdiff(preds, m$rpar)
  terms <- lapply(preds, function(p) {
    if (p %in% r$smooth) call("s", .sym(p)) else .sym(p)
  })
  for (pr in r$interactions) {
    pr <- unlist(pr)
    terms <- c(terms, list(call(":", .sym(pr[1]), .sym(pr[2]))))
  }
  if (!is.null(r$exposure)) {
    terms <- c(terms, list(call("offset", call("log", .sym(r$exposure)))))
  }
  if (!is.null(r$log_offset)) {
    terms <- c(terms, list(call("offset", .sym(r$log_offset))))
  }
  if (m$estimator == "countreg_hier" && !is.null(r$group)) {
    g <- .sym(r$group)
    lhs <- .plus(c(list(1), lapply(m$re_slopes, .sym)))
    bar <- if (length(m$re_slopes) && !isTRUE(m$re_correlated)) "||" else "|"
    terms <- c(terms, list(call("(", call(bar, lhs, g))))
    if (!is.null(r$group2)) {
      g2 <- .sym(r$group2)
      if (identical(m$group2_type, "nested")) g2 <- call(":", g, g2)
      terms <- c(terms, list(call("(", call("|", 1, g2))))
    }
  }
  rhs <- .plus(terms)
  # A random intercept replaces the fixed one
  if (rp && isTRUE(m$rpar_intercept)) rhs <- call("-", rhs, 1)
  call("~", .sym(r$response), rhs)
}

# Expressions that fit the model into `fit_name` (and its summary into
# `<fit_name>_summary`), from the data object of the specification. When
# the recipe converts columns to factors, a copy `<fit_name>_data` is made
# first.
compile_fit <- function(spec, fit_name) {
  r <- spec$recipe
  m <- spec$model
  out <- list()
  dat <- .sym(spec$data$object)
  if (length(r$factors)) {
    fd <- .sym(paste0(fit_name, "_data"))
    out$copy <- call("<-", fd, dat)
    for (f in r$factors) {
      out[[paste0("factor_", f)]] <- bquote(
        .(fd)[[.(f)]] <- factor(.(fd)[[.(f)]]))
    }
    dat <- fd
  }
  args <- list(formula = compile_formula(spec))
  est <- m$estimator
  fun <- call("::", .sym("flexCountReg"), .sym(est))
  if (est %in% c("renb", "poisLind.re", "repois", "negmultinom")) {
    args$group_var <- r$group
    args$data <- dat
    if (est == "negmultinom") args$family <- m$nm_family
  } else if (est == "countreg_fe") {
    args$fixef <- .formula_call(NULL, lapply(c(r$group, r$group2), .sym))
    args$data <- dat
    args$family <- .cvec(unname(m$family))
  } else {
    if (est == "countreg_rp") {
      rterms <- lapply(m$rpar, .sym)
      args$rpar_formula <- if (isTRUE(m$rpar_intercept)) {
        .formula_call(NULL, if (length(rterms)) rterms else list())
      } else {
        call("~", .plus(c(list(call("-", 1)), rterms)))
      }
    }
    args$data <- dat
    args$family <- .cvec(unname(m$family))
    if (length(m$dis1)) {
      args$dis_param_formula_1 <- .formula_call(NULL, lapply(m$dis1, .sym))
    }
  }
  if (!is.null(r$weights)) args$weights <- r$weights
  if (identical(m$weight_type, "sampling")) args$weight_type <- "sampling"
  if (length(m$cre)) {
    args$cre <- .formula_call(NULL, lapply(m$cre, .sym))
  }
  if (est == "countreg_hier" && !identical(m$re_dist, "normal")) {
    args$re_dist <- m$re_dist
  }
  if (est == "countreg_rp") {
    if (!identical(m$rpar_dist, "n")) args$rpardists <- m$rpar_dist
    if (isTRUE(m$correlated)) args$correlated <- TRUE
    args$ndraws <- as.integer(m$ndraws)
    if (!is.null(r$group)) args$group_var <- r$group
  }
  if (.app_is_mixture(est)) {
    lc <- est == "countreg_lc"
    if (length(m$family) == 1) {
      args[[if (lc) "classes" else "components"]] <- as.integer(m$components)
    }
    if (length(m$mixing)) {
      args[[if (lc) "membership_formula" else "mixing_formula"]] <-
        .formula_call(NULL, lapply(m$mixing, .sym))
    }
    if (length(m$shared_mean) || isTRUE(m$shared_dis1)) {
      sh <- list()
      if (length(m$shared_mean)) sh$mean <- .cvec(m$shared_mean)
      if (isTRUE(m$shared_dis1)) sh$dis1 <- TRUE
      args$shared <- as.call(c(.sym("list"), sh))
    }
    if (!is.null(r$group)) args$group_var <- r$group
    args$n_starts <- as.integer(m$n_starts)
    args$seed <- as.integer(m$seed)
  }
  if (est == "countreg_fe") {
    # countreg_fe(): `method` is the estimation algorithm (automatic here);
    # standard errors are cluster-robust or Hessian-based
    if (isTRUE(m$fe_jackknife)) args$jackknife <- TRUE
    args$stderr <- if (identical(m$stderr, "cluster")) "cluster" else
      "hessian"
    if (identical(m$stderr, "cluster") && !is.null(r$cluster)) {
      args$cluster <- .formula_call(NULL, list(.sym(r$cluster)))
    }
    args$max.iters <- as.integer(m$max_iters)
  } else {
    args$method <- if (est == "countreg_hier") {
      if (identical(m$method, "NM")) "Nelder-Mead" else "BFGS"
    } else m$method
    args$max.iters <- as.integer(m$max_iters)
    if (!identical(m$stderr, "normal")) args$stderr <- m$stderr
    if (identical(m$stderr, "cluster")) {
      args$cluster <- .formula_call(NULL, list(.sym(r$cluster)))
    }
  }
  fit_sym <- .sym(fit_name)
  out$fit <- call("<-", fit_sym, as.call(c(fun, args)))
  out$summary <- call("<-", .sym(paste0(fit_name, "_summary")),
                      call("summary", fit_sym))
  out
}

# ---- Outputs ------------------------------------------------------------------------

# A data.frame() call that rebuilds a small table (scenarios)
compile_data_frame <- function(df) {
  cols <- lapply(df, function(v) {
    if (is.factor(v)) {
      call("factor", .cvec(as.character(v)), levels = .cvec(levels(v)))
    } else {
      as.call(c(.sym("c"), as.list(v)))
    }
  })
  as.call(c(.sym("data.frame"), cols, check.names = FALSE))
}

# The rows predicted: an edited scenario table (<out>_scenarios), or new
# data prepared for the fit (<out>_data, from req$prep)
.scenario_sym <- function(req, out) {
  .sym(paste0(out, if (!is.null(req$prep)) "_data" else "_scenarios"))
}
.scenario_exprs <- function(req, sc) {
  if (!is.null(req$prep)) return(req$prep)
  list(call("<-", sc, compile_data_frame(req$scenarios)))
}

# Expressions of one output request on a fitted model. `req` is a list with
# `type` and its settings; `out` is the name of the result object.
#   rootogram  - max_count (NULL: largest observed count)
#   quantile   - randomized quantile residuals (seed 1) and their plots
#   cure       - var (NULL: fitted values)
#   smooths    - the estimated smooth effects (plot_smooths())
#   predict    - scenarios (data frame), interval, level, method, mixture,
#                seed
#   tolerance  - scenarios, content, confidence, side, scope, tol_method,
#                seed
#   compare    - fits (object names), labels
#   margeff    - measure, indicator, se, bootstraps, vars (NULL: the
#                defaults of margEffTable()), seed (bootstrap only)
# `seed` (frozen by app_prepare_output()) is passed to predict() and
# tolerance_interval(), so a parametric bootstrap behind an interval gives
# the same limits whenever the code is run. An interval prediction is kept
# whole in `<out>_interval` (its "interval" attribute records the
# uncertainty method and the bootstrap refits accepted) before it is bound
# to the scenarios in `<out>`.
compile_output <- function(req, fit_name, out) {
  f <- .sym(fit_name)
  o <- .sym(out)
  switch(req$type,
    rootogram = {
      a <- list(f)
      if (!is.null(req$max_count)) a$max_count <- as.integer(req$max_count)
      a$overflow <- TRUE
      a$main <- req$title %||% fit_name
      list(call("<-", o, as.call(c(quote(flexCountReg::rootogram), a))))
    },
    quantile = list(
      bquote(.(o) <- {
        set.seed(1)
        stats::residuals(.(f), type = "quantile")
      }),
      bquote(graphics::par(mfrow = c(1, 2))),
      bquote(graphics::plot(stats::fitted(.(f)), .(o), xlab = "Fitted mean",
                            ylab = "Quantile residual")),
      quote(graphics::abline(h = 0, lty = 2)),
      bquote(stats::qqnorm(.(o))),
      bquote(stats::qqline(.(o), lty = 2)),
      quote(graphics::par(mfrow = c(1, 1)))),
    smooths = list(
      bquote(graphics::par(mfrow = grDevices::n2mfrow(
        length(.(f)$model$smooth_info)))),
      bquote(.(o) <- flexCountReg::plot_smooths(.(f), seed = 1)),
      quote(graphics::par(mfrow = c(1, 1)))),
    cure = {
      a <- list(f)
      if (!is.null(req$var)) a$indvar <- req$var
      list(call("<-", o, as.call(c(quote(flexCountReg::cureplot), a))),
           call("print", o))
    },
    predict = {
      sc <- .scenario_sym(req, out)
      pa <- list(f, newdata = sc)
      if (!is.null(req$method)) pa$method <- req$method
      pre <- list()
      if (!identical(req$interval %||% "none", "none")) {
        pa$interval <- req$interval
        pa$level <- req$level
        if (!is.null(req$seed)) pa$seed <- as.integer(req$seed)
        iv <- .sym(paste0(out, "_interval"))
        pre <- list(call("<-", iv, as.call(c(quote(stats::predict), pa))))
        pred <- call("cbind", sc, iv)
      } else if (isTRUE(req$mixture)) {
        pred <- as.call(list(
          .sym("data.frame"), sc,
          expected_count = as.call(c(quote(stats::predict), pa)),
          mean = as.call(c(quote(stats::predict), pa,
                           list(type = "component"))),
          prob = as.call(list(quote(flexCountReg::component_probs), f,
                              newdata = sc)),
          check.names = FALSE))
      } else {
        pred <- call("cbind", sc,
                     expected_count = as.call(c(quote(stats::predict), pa)))
      }
      c(.scenario_exprs(req, sc), pre, list(call("<-", o, pred)))
    },
    tolerance = {
      sc <- .scenario_sym(req, out)
      ta <- list(f, newdata = sc, content = req$content,
                 confidence = req$confidence, side = req$side,
                 scope = req$scope, method = req$tol_method)
      if (!is.null(req$seed)) ta$seed <- as.integer(req$seed)
      c(.scenario_exprs(req, sc),
        list(call("<-", o, as.call(c(quote(flexCountReg::tolerance_interval),
                                     ta)))))
    },
    margeff = {
      # digits = Inf keeps the numbers unrounded; the page rounds them for
      # display only. A bootstrap draws resamples, so its frozen seed is
      # set first.
      a <- list(f)
      if (length(req$vars)) a$vars <- .cvec(req$vars)
      a$measure <- req$measure %||% "auto"
      a$indicator <- req$indicator %||% "pseudo_elasticity"
      a$se <- req$se %||% "delta"
      if (identical(a$se, "bootstrap")) {
        a$bootstraps <- as.integer(req$bootstraps %||% 200L)
      }
      a$digits <- Inf
      c(if (!is.null(req$seed)) list(call("set.seed", as.integer(req$seed))),
        list(call("<-", o, as.call(c(quote(flexCountReg::margEffTable),
                                     a)))))
    },
    compare = {
      fl <- as.call(c(.sym("list"), stats::setNames(
        lapply(req$fits, .sym), req$labels)))
      list(bquote(.(o) <- local({
        fits <- .(fl)
        data.frame(
          fit = names(fits),
          model = vapply(fits, function(m) class(m)[1], ""),
          logLik = vapply(fits, function(m) as.numeric(stats::logLik(m)), 0),
          df = vapply(fits, function(m) attr(stats::logLik(m), "df"), 0),
          AIC = vapply(fits, stats::AIC, 0),
          BIC = vapply(fits, stats::BIC, 0),
          units = vapply(fits, function(m) attr(stats::logLik(m), "nobs"), 0),
          row.names = NULL)
      })))
    },
    stop("Unknown output type ", req$type, call. = FALSE))
}

# ---- Session script -------------------------------------------------------------------

.deparse <- function(e) paste(deparse(e, width.cutoff = 72L), collapse = "\n")

# The session script: every logged entry (list(title, display)) in order
session_script <- function(log, versions = fcr_versions()) {
  hdr <- c(
    "# R code of an analysis run in the flexCountReg app.",
    "# It reproduces everything the app ran, in order, without the app;",
    "# only flexCountReg (and its dependencies) are needed.",
    paste0("# Versions: ", paste(names(versions), versions, sep = " ",
                                 collapse = ", ")),
    "", "library(flexCountReg)", "")
  body <- unlist(lapply(log, function(e) {
    # A title is a one-line comment, whatever it contains
    title <- gsub("[[:space:]]+", " ", e$title)
    c(paste0("# ---- ", title, " ----"),
      vapply(e$display, .deparse, character(1)), "")
  }))
  paste(c(hdr, body), collapse = "\n")
}

# The code a draft fit would add (shown before it is run)
draft_script <- function(spec, fit_name = "fit_next") {
  paste(vapply(compile_fit(spec, fit_name), .deparse, character(1)),
        collapse = "\n")
}

# Evaluate expressions in an environment, returning it
run_compiled <- function(exprs, envir = new.env(parent = globalenv())) {
  for (e in exprs) eval(e, envir)
  envir
}
