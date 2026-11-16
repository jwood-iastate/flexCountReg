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
  unique(unlist(lapply(fs, function(f) all.vars(stats::as.formula(f)))))
}

# TRUE when a formula contains an offset() term.
.fcr_has_formula_offset <- function(formula) {
  !is.null(attr(stats::terms(stats::as.formula(formula)), "offset"))
}

# Right-hand-side design matrix for an auxiliary formula. Rows with missing
# values are kept as NA rows (na.pass), so the caller controls row filtering.
.fcr_rhs_matrix <- function(fml, mf_data, xlev = NULL,
                            drop_intercept = FALSE) {
  if (is.null(fml)) return(NULL)
  tt <- stats::delete.response(stats::terms(stats::as.formula(fml)))
  mf <- stats::model.frame(tt, mf_data, na.action = stats::na.pass,
                           xlev = xlev)
  X <- stats::model.matrix(tt, mf)
  if (drop_intercept) X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  attr(X, "xlevels") <- stats::.getXlevels(tt, mf)
  X
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
  formula <- stats::as.formula(formula)
  aux_formulas <- Filter(Negate(is.null), aux_formulas)

  formula_offset <- .fcr_has_formula_offset(formula)
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

  vars <- unique(c(all.vars(formula),
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

  mf <- stats::model.frame(formula, d, na.action = stats::na.fail)
  tt <- stats::terms(mf)
  X <- stats::model.matrix(tt, mf)

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
  aux <- Map(function(f, nm) .fcr_rhs_matrix(f, d,
                                             drop_intercept = nm %in% het),
             aux_formulas, names(aux_formulas))

  w <- if (is.null(weights)) {
    rep(1, nrow(d))
  } else if (is.character(weights)) {
    as.numeric(d[[weights]])
  } else {
    as.numeric(weights[keep])
  }
  if (any(w < 0)) stop("`weights` must be non-negative.", call. = FALSE)

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

  list(y = stats::model.response(mf), X = X, offset = off, weights = w,
       group = g, aux = aux, rows_used = which(keep),
       n_dropped = sum(!keep),
       blueprint = list(terms = tt,
                        xlevels = stats::.getXlevels(tt, mf),
                        contrasts = attr(X, "contrasts"),
                        offset_cols = offset,
                        offset_label = offset_label,
                        aux_formulas = aux_formulas,
                        aux_xlevels = lapply(aux, attr, "xlevels"),
                        x_names = colnames(X)))
}

# Rebuild model matrices for new data from a blueprint. Returns matrices with
# exactly nrow(new_data) rows; rows with missing values become NA rows.
forge_new_data <- function(blueprint, new_data) {
  new_data <- as.data.frame(new_data)
  tt <- stats::delete.response(blueprint$terms)
  need <- unique(c(all.vars(tt), blueprint$offset_cols,
                   unlist(lapply(blueprint$aux_formulas, all.vars))))
  miss <- setdiff(need, names(new_data))
  if (length(miss)) {
    stop("`newdata` is missing required column(s): ", toString(miss),
         call. = FALSE)
  }

  mf <- stats::model.frame(tt, new_data, na.action = stats::na.pass,
                           xlev = blueprint$xlevels)
  X <- stats::model.matrix(tt, mf, contrasts.arg = blueprint$contrasts)
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
  aux <- Map(function(f, nm) {
    .fcr_rhs_matrix(f, new_data, blueprint$aux_xlevels[[nm]], nm %in% het)
  }, blueprint$aux_formulas, names(blueprint$aux_formulas))

  list(X = X, offset = as.numeric(off), aux = aux)
}
