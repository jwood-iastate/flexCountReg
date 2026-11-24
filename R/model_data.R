# Shared data preparation for estimation and prediction.
#
# `prepare_model_data()` builds every model matrix used by an estimator from a
# single, complete-case subset of `data`, so the outcome, the main design
# matrix, the offset, the weights, the group id and every auxiliary design
# matrix always refer to the same rows. It also returns a "blueprint" that
# `forge_new_data()` uses to rebuild the same columns from new data at
# prediction time.

#' @importFrom stats as.formula terms delete.response model.frame model.matrix
#' @importFrom stats model.response model.offset complete.cases na.pass na.fail
#' @importFrom stats .getXlevels
NULL

# All variable names used in a set of (possibly NULL) formulas.
.fcr_vars <- function(...) {
  fs <- Filter(Negate(is.null), list(...))
  unique(unlist(lapply(fs, .fcr_formula_vars)))
}

# TRUE when a formula contains an offset() term.
.fcr_has_formula_offset <- function(formula) {
  !is.null(attr(stats::terms(stats::as.formula(formula)), "offset"))
}

# Right-hand-side design matrix for an auxiliary formula, or for the trained
# terms of one (which carry the data-dependent bases of poly(), ns(),
# scale() etc. in their "predvars"). Rows with missing values are kept as NA
# rows (na.pass), so the caller controls row filtering. The trained terms
# are returned in attribute "terms".
.fcr_rhs_matrix <- function(fml, mf_data, xlev = NULL,
                            drop_intercept = FALSE, contrasts = NULL) {
  if (is.null(fml)) return(NULL)
  tt <- if (inherits(fml, "terms")) fml else
    stats::delete.response(stats::terms(stats::as.formula(fml)))
  mf <- .fcr_model_frame_xlev(tt, mf_data, xlev)
  X <- stats::model.matrix(tt, mf, contrasts.arg = contrasts)
  ctr <- attr(X, "contrasts")
  if (drop_intercept) X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  attr(X, "xlevels") <- stats::.getXlevels(tt, mf)
  attr(X, "contrasts") <- ctr
  attr(X, "terms") <- attr(mf, "terms")
  X
}

# model.frame() with training factor levels. Re-applying the levels drops any
# contrasts attribute of a factor in the new data (with a warning); the
# training contrasts are passed to model.matrix() afterwards, so that
# warning is expected and muffled.
.fcr_model_frame_xlev <- function(tt, data, xlev) {
  withCallingHandlers(
    stats::model.frame(tt, data, na.action = stats::na.pass, xlev = xlev),
    warning = function(w) {
      if (grepl("contrasts dropped from factor", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    })
}

# Expand a `.` in a formula against the columns of `data` (as lm() does).
.fcr_expand_dot <- function(fml, data) {
  if (is.null(fml)) return(NULL)
  fml <- stats::as.formula(fml)
  if (!"." %in% all.vars(fml)) return(fml)
  out <- stats::formula(stats::terms(fml, data = data))
  environment(out) <- environment(fml)
  out
}

# Prepare all data used by an estimator.
#
# formula       main model formula; may contain offset() terms.
# data          data frame.
# aux_formulas  named list of one-sided formulas (NULL entries allowed).
# offset        character vector of offset column names (log scale), or NULL.
#               Only one of `offset` and a formula offset() may be used.
# weights       NULL, a column name, or a numeric vector of length nrow(data).
# group_var     NULL or column name(s) identifying clusters.
prepare_model_data <- function(formula, data, aux_formulas = list(),
                               offset = NULL, weights = NULL,
                               group_var = NULL) {
  data <- as.data.frame(data)
  # Expand `.` before any checks; the expanded forms go into the blueprint
  formula <- .fcr_expand_dot(formula, data)
  aux_formulas <- lapply(Filter(Negate(is.null), aux_formulas),
                         .fcr_expand_dot, data = data)

  # Smooth terms (s(), te(), ...): the parametric part goes through
  # model.frame() as before; the smooths are set up below
  main_specs <- NULL
  formula_pf <- formula
  if (.fcr_has_smooth(formula)) {
    sp <- .fcr_smooth_split(formula)
    formula_pf <- sp$pf
    main_specs <- sp$specs
  }
  formula_offset <- .fcr_has_formula_offset(formula_pf)
  if (!is.null(offset) && formula_offset) {
    stop("Specify the offset either with offset() in the formula or with the ",
         "`offset` argument, not both.", call. = FALSE)
  }
  if (!is.null(offset)) {
    if (!is.character(offset)) {
      stop("`offset` must be a character vector of column names.",
           call. = FALSE)
    }
    miss <- setdiff(offset, names(data))
    if (length(miss)) {
      stop("Offset column(s) not found in `data`: ", toString(miss),
           call. = FALSE)
    }
  }

  w_col <- NULL
  if (!is.null(weights)) {
    if (is.character(weights)) {
      if (length(weights) != 1L) {
        stop("`weights` must be a single column name or a numeric vector.",
             call. = FALSE)
      }
      w_col <- weights
    } else if (is.numeric(weights)) {
      if (length(weights) != nrow(data)) {
        stop("A numeric `weights` vector must have one value per row of ",
             "`data`.", call. = FALSE)
      }
    } else {
      stop("`weights` must be a column name or a numeric vector.",
           call. = FALSE)
    }
  }

  vars <- unique(c(.fcr_formula_vars(formula),
                   do.call(.fcr_vars, unname(aux_formulas)),
                   offset, w_col, group_var))
  miss <- setdiff(vars, names(data))
  if (length(miss)) {
    stop("Variable(s) not found in `data`: ", toString(miss), call. = FALSE)
  }

  # One row filter for everything
  keep <- stats::complete.cases(data[, vars, drop = FALSE])
  if (is.numeric(weights)) keep <- keep & !is.na(weights)
  if (!any(keep)) stop("No complete cases in `data`.", call. = FALSE)
  d <- data[keep, , drop = FALSE]

  # Missing raw values were removed above; values made invalid by a
  # transformation (e.g. log of a negative) are reported by the checks below
  mf <- stats::model.frame(formula_pf, d, na.action = stats::na.pass)
  y <- stats::model.response(mf)
  if (!is.null(y) && (!is.numeric(y) || is.matrix(y))) {
    stop("The response must be a numeric count vector, not ",
         class(y)[1], ". Convert factors with as.numeric(as.character()).",
         call. = FALSE)
  }
  if (!is.null(y)) .fcr_check_counts(y)
  tt <- stats::terms(mf)
  X <- stats::model.matrix(tt, mf)
  .fcr_check_design(X, "formula")
  smooth <- list(main = list(), aux = list())
  if (length(main_specs)) {
    built <- .fcr_smooth_build(main_specs, d)
    smooth$main <- .fcr_smooth_shift(built$smooths, ncol(X))
    contr <- attr(X, "contrasts")
    X <- cbind(X, built$B)
    attr(X, "contrasts") <- contr
  }

  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, nrow(d))
  if (!is.null(offset)) {
    off <- off + rowSums(as.matrix(d[, offset, drop = FALSE]))
  }
  off <- as.numeric(off)
  if (any(!is.finite(off))) {
    stop("Offsets must be finite (use the log of the exposure, e.g. ",
         "log(length)).", call. = FALSE)
  }

  het <- c("het_mean", "het_var")
  aux <- Map(function(f, nm) {
    if (!.fcr_has_smooth(f)) {
      m <- .fcr_rhs_matrix(f, d, drop_intercept = nm %in% het)
      .fcr_check_design(m, nm)
      return(m)
    }
    sp <- .fcr_smooth_split(f)
    m <- .fcr_rhs_matrix(sp$pf, d, drop_intercept = nm %in% het)
    .fcr_check_design(m, nm)
    built <- .fcr_smooth_build(sp$specs, d)
    smooth$aux[[nm]] <<- .fcr_smooth_shift(built$smooths, ncol(m))
    keep_attr <- attributes(m)[c("xlevels", "contrasts", "terms")]
    m <- cbind(m, built$B)
    attributes(m)[names(keep_attr)] <- keep_attr
    m
  }, aux_formulas, names(aux_formulas))

  w <- if (is.null(weights)) {
    rep(1, nrow(d))
  } else if (is.character(weights)) {
    as.numeric(d[[weights]])
  } else {
    as.numeric(weights[keep])
  }
  if (any(!is.finite(w))) stop("`weights` must be finite.", call. = FALSE)
  if (any(w < 0)) stop("`weights` must be non-negative.", call. = FALSE)
  if (sum(w) <= 0) {
    stop("`weights` must have a positive total.", call. = FALSE)
  }

  g <- NULL
  if (!is.null(group_var)) {
    g <- droplevels(interaction(d[, group_var, drop = FALSE], drop = TRUE,
                                sep = ":"))
  }

  offset_label <- if (formula_offset) {
    vapply(attr(tt, "offset"), function(i) deparse(attr(tt, "variables")[[i + 1]]),
           character(1))
  } else {
    offset
  }

  list(y = y, X = X, offset = off, weights = w,
       formula = formula, aux_formulas = aux_formulas,
       group = g, aux = aux, rows_used = which(keep),
       n_dropped = sum(!keep), smooth = smooth,
       blueprint = list(terms = tt,
                        xlevels = stats::.getXlevels(tt, mf),
                        contrasts = attr(X, "contrasts"),
                        offset_cols = offset,
                        offset_label = offset_label,
                        aux_formulas = aux_formulas,
                        aux_xlevels = lapply(aux, attr, "xlevels"),
                        aux_contrasts = lapply(aux, attr, "contrasts"),
                        aux_terms = lapply(aux, attr, "terms"),
                        x_names = colnames(X),
                        smooth = smooth))
}

# Rebuild model matrices for new data from a blueprint. Returns matrices with
# exactly nrow(new_data) rows; rows with missing values become NA rows.
forge_new_data <- function(blueprint, new_data) {
  new_data <- as.data.frame(new_data)
  tt <- stats::delete.response(blueprint$terms)
  need <- unique(c(all.vars(tt), blueprint$offset_cols,
                   unlist(lapply(blueprint$aux_formulas, .fcr_formula_vars))))
  miss <- setdiff(need, names(new_data))
  if (length(miss)) {
    stop("`newdata` is missing required column(s): ", toString(miss),
         call. = FALSE)
  }

  mf <- .fcr_model_frame_xlev(tt, new_data, blueprint$xlevels)
  X <- stats::model.matrix(tt, mf, contrasts.arg = blueprint$contrasts)
  if (length(blueprint$smooth$main)) {
    X <- cbind(X, .fcr_smooth_predict(blueprint$smooth$main, new_data))
  }
  if (!identical(colnames(X), blueprint$x_names)) {
    stop("`newdata` produces different model columns than the training ",
         "data.", call. = FALSE)
  }

  off <- stats::model.offset(mf)
  if (is.null(off)) off <- rep(0, nrow(new_data))
  if (!is.null(blueprint$offset_cols)) {
    off <- off + rowSums(as.matrix(new_data[, blueprint$offset_cols,
                                            drop = FALSE]))
  }

  het <- c("het_mean", "het_var")
  # Trained terms (data-dependent bases fixed at their training values);
  # blueprints from earlier versions only have the formulas
  aux <- Map(function(f, nm) {
    m <- .fcr_rhs_matrix(blueprint$aux_terms[[nm]] %||% f, new_data,
                         blueprint$aux_xlevels[[nm]], nm %in% het,
                         contrasts = blueprint$aux_contrasts[[nm]])
    attr(m, "terms") <- NULL
    sm <- blueprint$smooth$aux[[nm]]
    if (length(sm)) m <- cbind(m, .fcr_smooth_predict(sm, new_data))
    m
  }, blueprint$aux_formulas, names(blueprint$aux_formulas))

  list(X = X, offset = as.numeric(off), aux = aux)
}

# The response must be finite, non-negative whole numbers.
.fcr_check_counts <- function(y) {
  y <- as.numeric(y)
  if (any(!is.finite(y))) {
    stop("The response must be finite counts (it has infinite or NaN ",
         "values).", call. = FALSE)
  }
  if (any(y < 0)) {
    stop("The response must be non-negative counts (", sum(y < 0),
         " negative value(s)).", call. = FALSE)
  }
  if (any(abs(y - round(y)) > sqrt(.Machine$double.eps))) {
    stop("The response must be whole numbers (counts); ",
         sum(abs(y - round(y)) > sqrt(.Machine$double.eps)),
         " value(s) are not.", call. = FALSE)
  }
  invisible(NULL)
}

# Design matrices after all formula transformations: finite values (a
# transformation such as log() of a non-positive value gives NaN or -Inf)
# and full column rank (no aliased columns).
.fcr_check_design <- function(X, what) {
  if (is.null(X) || !ncol(X)) return(invisible(NULL))
  bad <- colnames(X)[colSums(!is.finite(X)) > 0]
  if (length(bad)) {
    stop("The design matrix of `", what, "` has non-finite values in ",
         "column(s) ", toString(bad), " (e.g. log() of a non-positive ",
         "value).", call. = FALSE)
  }
  qx <- qr(X)
  if (qx$rank < ncol(X)) {
    aliased <- colnames(X)[qx$pivot[seq(qx$rank + 1, ncol(X))]]
    stop("The design matrix of `", what, "` is rank deficient: ",
         toString(aliased), " ", if (length(aliased) == 1) "is" else "are",
         " a linear combination of other columns (aliased). Remove ",
         if (length(aliased) == 1) "it" else "them", ".", call. = FALSE)
  }
  invisible(NULL)
}
