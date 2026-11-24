# Time-varying effects (temporal instability): the tv() formula term, its
# bases, and the summaries temporal_effects() and temporal_tests().
#
# tv(x, time) returns a matrix of basis columns, one column per basis
# function, so it works in every formula that goes through model.frame()
# and model.matrix(). Its makepredictcall() method stores the training time
# values in the terms' "predvars", so prediction for new data uses the same
# columns. The column names are made readable afterwards by
# .fcr_tv_finish(): "tv(lnaadt, time = Year)Year2011" becomes
# "lnaadt:Year2011". tv(x, time, basis = "smooth") is rewritten to the
# varying-coefficient smooth s(time, by = x) before the model is set up.

#' Time-varying effects in model formulas
#'
#' `tv()` lets the effect of a variable change over time (temporal
#' instability). Used in a model formula, `tv(x, time = year)` replaces the
#' single coefficient of `x` by one effect for each value of `year`,
#' parameterized by a basis chosen for the question being asked. It can be
#' used in the mean formula and in the formulas of distribution parameters,
#' zero models and the other auxiliary formulas of every estimator, and
#' `tv(1, time = year)` gives a time-varying intercept.
#'
#' @param x the numeric variable whose effect varies, or `1` for the
#'   intercept.
#' @param time the time variable (numeric, such as a year, or a factor).
#'   Its distinct values in the estimation data define the periods, in
#'   increasing order (factor level order for a factor).
#' @param basis how the period effects \eqn{\beta_t} are parameterized:
#'   \describe{
#'   \item{`"step"` (default)}{the effect in the first period and the change
#'     at each later period: \eqn{\beta_t = \beta_1 + \sum_{k=2}^t
#'     \Delta_k}. A zero change \eqn{\Delta_k} means equal effects in
#'     periods \eqn{k-1} and \eqn{k}. Columns `x:time<first>` and
#'     `x:time_step<k>`.}
#'   \item{`"level"`}{a separate effect for each period (columns
#'     `x:time<t>`). A zero coefficient means no effect in that period.}
#'   \item{`"hinge"`}{the cumulative "spline indicators" of Shabab et al.
#'     (2026): \eqn{\beta_t = \sum_k \max(t - k + 1, 0)\,\delta_k}, with
#'     \eqn{t} the period's position (1, 2, ...). A zero \eqn{\delta_k}
#'     means that the trend of the earlier periods continues at period
#'     \eqn{k}, not that the effect stays the same: dropping \eqn{\delta_2}
#'     imposes \eqn{\beta_2 = 2\beta_1}. Columns `x:time_hinge<k>`.}
#'   \item{`"smooth"`}{a penalized spline in time, \eqn{\beta(t)}, the
#'     varying-coefficient smooth `s(time, by = x)` of mgcv, so the data
#'     choose between a constant, a linear trend and a flexible path.
#'     `time` must be numeric. Arguments in `...` (such as `k` and `bs`)
#'     are passed to [mgcv::s()]; `k` defaults to the number of periods, at
#'     most 10.}
#'   }
#'   The step, level and hinge bases give the same fit (the same
#'   log-likelihood) and the same period effects; they differ in what the
#'   individual coefficients and their tests mean.
#' @param levels the periods. Set from the estimation data; give it only to
#'   fix the periods in advance.
#' @param new_time what to do at prediction for a time value that was not in
#'   the estimation data: `"error"` (default), or `"nearest"` to use the
#'   effect of the nearest period (numeric time only), for example the last
#'   period's effect for a later year.
#' @param ... arguments passed to [mgcv::s()] with `basis = "smooth"`.
#'
#' @details
#' With an intercept in the model, `tv(1, time = year)` leaves out its first
#' column: the intercept is the first period's level, and the remaining
#' columns are the changes (`"step"`), the differences from the first
#' period (`"level"`) or the hinge terms from the second period
#' (`"hinge"`). A variable should appear either alone or in `tv()`, not
#' both.
#'
#' Use [temporal_effects()] for the effect in each period with its standard
#' error, and [temporal_tests()] to test whether and when the effects
#' change.
#'
#' @return A numeric matrix with one column per basis function (it is meant
#'   to be used inside a model formula).
#'
#' @references
#' Shabab, K. R., Bhowmik, T., Zaki, M. H. and Eluru, N. (2026). A joint
#' count and generalized ordered fractional split approach for addressing
#' temporal instability in road safety analysis. *Analytic Methods in
#' Accident Research*, 51, 100440.
#'
#' @seealso [temporal_effects()], [temporal_tests()]
#' @examples
#' data("washington_roads")
#' m <- countreg(Total_crashes ~ tv(lnaadt, time = Year) + speed50 +
#'                 offset(lnlength), data = washington_roads,
#'               family = "NB2", method = "BFGS")
#' temporal_effects(m)
#' temporal_tests(m)
#' @export
tv <- function(x, time, basis = c("step", "level", "hinge", "smooth"),
               levels = NULL, new_time = c("error", "nearest"), ...) {
  basis <- match.arg(basis)
  new_time <- match.arg(new_time)
  if (basis == "smooth") {
    stop("tv(basis = \"smooth\") can only be used in a model formula of ",
         "this package, where it becomes s(time, by = x).", call. = FALSE)
  }
  tname <- paste(deparse(substitute(time)), collapse = "")
  if (!is.numeric(x) && !is.logical(x)) {
    stop("tv(): `x` must be numeric (make separate columns for the levels ",
         "of a factor).", call. = FALSE)
  }
  x <- as.numeric(x)
  if (length(x) == 1L) x <- rep(x, length(time))
  if (length(x) != length(time)) {
    stop("tv(): `x` and `time` have different lengths.", call. = FALSE)
  }
  if (is.null(levels)) {
    levels <- if (is.factor(time)) {
      base::levels(droplevels(time))
    } else {
      sort(unique(time[!is.na(time)]))
    }
  }
  if (length(levels) < 2L) {
    stop("tv(): `time` has fewer than two distinct values.", call. = FALSE)
  }
  pos <- .fcr_tv_position(time, levels, new_time)
  spec <- .fcr_tv_spec(levels, basis, tname)
  out <- spec$B[pos, , drop = FALSE] * x
  dimnames(out) <- list(NULL, spec$names)
  structure(out, class = c("tv", "matrix", "array"), levels = levels,
            basis = basis)
}

# Positions of time values among the periods; NA time gives NA.
.fcr_tv_position <- function(time, levels, new_time = "error") {
  num <- is.numeric(levels)
  pos <- if (num) {
    match(as.numeric(time), levels)
  } else {
    match(as.character(time), as.character(levels))
  }
  bad <- is.na(pos) & !is.na(time)
  if (any(bad)) {
    if (new_time == "nearest" && num) {
      tv <- as.numeric(time[bad])
      pos[bad] <- vapply(tv, function(t) which.min(abs(levels - t)),
                         integer(1))
    } else {
      stop("tv(): time value(s) not in the estimation data: ",
           toString(utils::head(unique(time[bad]), 5)),
           if (num) ". Use tv(..., new_time = \"nearest\") to use the nearest period's effect.",
           call. = FALSE)
    }
  }
  pos
}

# Basis matrix (periods x columns) and column names of a basis.
.fcr_tv_spec <- function(levels, basis, tname) {
  nt <- length(levels)
  lab <- as.character(levels)
  t <- seq_len(nt)
  switch(basis,
    step = list(B = outer(t, t, `>=`) * 1,
                names = c(paste0(tname, lab[1]),
                          paste0(tname, "_step", lab[-1]))),
    level = list(B = diag(nt), names = paste0(tname, lab)),
    hinge = list(B = pmax(outer(t, t, function(a, b) a - b + 1), 0),
                 names = paste0(tname, "_hinge", lab)))
}

#' @importFrom stats makepredictcall
#' @exportS3Method stats::makepredictcall
makepredictcall.tv <- function(var, call) {
  fn <- call[[1L]]
  nm <- if (is.call(fn)) as.character(fn[[3L]]) else as.character(fn)
  if (!identical(nm, "tv")) return(NextMethod())
  call$levels <- attr(var, "levels")
  call
}

# Is an expression a call to tv()?
.fcr_is_tv <- function(e) {
  if (!is.call(e)) return(FALSE)
  fn <- e[[1L]]
  if (is.name(fn)) return(identical(as.character(fn), "tv"))
  is.call(fn) && identical(as.character(fn[[1L]]), "::") &&
    identical(as.character(fn[[3L]]), "tv")
}

# tv() call with its arguments matched by name.
.fcr_tv_match <- function(e) {
  e[[1L]] <- quote(tv)
  match.call(tv, e, expand.dots = TRUE)
}

.fcr_tv_basis_of <- function(mc) {
  b <- mc$basis
  if (is.null(b)) "step" else match.arg(eval(b), c("step", "level", "hinge", "smooth"))
}

# Is `e` the literal constant 1 (the intercept)?
.fcr_tv_is_one <- function(e) is.numeric(e) && length(e) == 1L && e == 1

# Rename the columns of tv() terms in a model matrix and drop the first
# column of tv(1, ...) when the model has an intercept.
.fcr_tv_finish <- function(X, tt) {
  vars <- as.list(attr(tt, "variables"))[-1L]
  is_tv <- vapply(vars, .fcr_is_tv, logical(1))
  if (!any(is_tv) || !ncol(X)) return(X)
  asg <- attr(X, "assign")
  ctr <- attr(X, "contrasts")
  labs <- attr(tt, "term.labels")
  cn <- colnames(X)
  drop <- integer(0)
  for (v in vars[is_tv]) {
    lab <- paste(deparse(v, width.cutoff = 500L), collapse = " ")
    mc <- .fcr_tv_match(v)
    one <- .fcr_tv_is_one(mc$x)
    pre <- if (one) "(Intercept):" else
      paste0(paste(deparse(mc$x, width.cutoff = 500L), collapse = " "), ":")
    hit <- grepl(lab, cn, fixed = TRUE)
    cn[hit] <- sub(lab, pre, cn[hit], fixed = TRUE)
    if (one && "(Intercept)" %in% colnames(X)) {
      ti <- match(lab, labs)
      if (!is.na(ti)) drop <- c(drop, which(asg == ti)[1L])
    }
  }
  colnames(X) <- cn
  if (length(drop)) {
    X <- X[, -drop, drop = FALSE]
    attr(X, "assign") <- asg[-drop]
    attr(X, "contrasts") <- ctr
  }
  X
}

# Rewrite tv(x, time, basis = "smooth", ...) as s(time, by = x, ...)
# (s(time, ...) for x = 1). With `data`, `k` defaults to the number of
# distinct time values, at most 10.
.fcr_tv_smooth_rewrite <- function(fml, data = NULL) {
  if (is.null(fml)) return(NULL)
  fml <- stats::as.formula(fml)
  walk <- function(e) {
    if (!is.call(e)) return(e)
    if (.fcr_is_tv(e)) {
      mc <- .fcr_tv_match(e)
      if (!identical(.fcr_tv_basis_of(mc), "smooth")) return(e)
      args <- as.list(mc)[-1L]
      extra <- args[setdiff(names(args),
                            c("x", "time", "basis", "levels", "new_time"))]
      if (is.null(extra$k) && !is.null(data)) {
        tvals <- eval(mc$time, data, environment(fml))
        extra$k <- min(10L, length(unique(tvals[!is.na(tvals)])))
      }
      s_args <- c(list(mc$time),
                  if (!.fcr_tv_is_one(mc$x)) list(by = mc$x), extra)
      return(as.call(c(list(quote(s)), s_args)))
    }
    for (i in seq_along(e)[-1L]) {
      if (!is.null(e[[i]])) e[[i]] <- walk(e[[i]])
    }
    e
  }
  out <- walk(fml)
  environment(out) <- environment(fml)
  out
}

# Does a formula contain tv(..., basis = "smooth")?
.fcr_has_tv_smooth <- function(fml) {
  if (is.null(fml)) return(FALSE)
  walk <- function(e) {
    if (!is.call(e)) return(FALSE)
    if (.fcr_is_tv(e) &&
        identical(.fcr_tv_basis_of(.fcr_tv_match(e)), "smooth")) {
      return(TRUE)
    }
    any(vapply(as.list(e)[-1L], walk, logical(1)))
  }
  walk(stats::as.formula(fml))
}

# Time variables of the tv() terms of formulas.
.fcr_tv_time_vars <- function(fmls) {
  out <- character(0)
  walk <- function(e) {
    if (!is.call(e)) return(invisible())
    if (.fcr_is_tv(e)) {
      out <<- c(out, all.vars(.fcr_tv_match(e)$time))
      return(invisible())
    }
    for (a in as.list(e)[-1L]) walk(a)
  }
  for (f in Filter(Negate(is.null), fmls)) walk(stats::as.formula(f))
  unique(out)
}

# ---------------------------------------------------------------------------
# Period effects and tests
# ---------------------------------------------------------------------------

.fcr_any_blueprint <- function(object) {
  object$blueprint %||% object$model$blueprint
}

# The tv() terms of a fit: one entry per term and coefficient block, with
# the period labels, the matrix A mapping the block's coefficients to the
# period effects, and the coefficient names.
.fcr_tv_terms <- function(object) {
  bp <- .fcr_any_blueprint(object)
  if (is.null(bp)) return(list())
  cf <- stats::coef(object)
  cn <- names(cf)
  term_sets <- c(list(bp$terms), unname(bp$aux_terms %||% list()))
  out <- list()
  seen <- character(0)
  for (tt in Filter(Negate(is.null), term_sets)) {
    pv <- as.list(attr(tt, "predvars") %||% attr(tt, "variables"))[-1L]
    for (v in pv[vapply(pv, .fcr_is_tv, logical(1))]) {
      mc <- .fcr_tv_match(v)
      lev <- eval(mc$levels)
      if (is.null(lev)) next
      basis <- .fcr_tv_basis_of(mc)
      tname <- paste(deparse(mc$time), collapse = "")
      one <- .fcr_tv_is_one(mc$x)
      xl <- if (one) "(Intercept)" else
        paste(deparse(mc$x, width.cutoff = 500L), collapse = " ")
      sp <- .fcr_tv_spec(lev, basis, tname)
      # Random parameters of countreg_rp() have the suffix ":Mean" (their
      # standard deviations belong to the basis columns, not the periods)
      for (sfx in c("", ":Mean")) {
        full <- paste0(xl, ":", sp$names, sfx)
        # Coefficient blocks holding these columns: the mean ("") or a
        # prefix such as "ln(alpha)"
        prefixes <- unique(c(
          if (any(full %in% cn)) "",
          unlist(lapply(full, function(f) {
            h <- cn[endsWith(cn, paste0(":", f))]
            substr(h, 1L, nchar(h) - nchar(f) - 1L)
          }))))
        for (pf in prefixes) {
          nm <- if (nzchar(pf)) paste0(pf, ":", full) else full
          present <- nm %in% cn
          if (!any(present)) next
          key <- paste(pf, sfx, xl, tname, basis)
          if (key %in% seen) next
          seen <- c(seen, key)
          B <- sp$B
          coef_names <- nm
          if (!present[1L] && one) {
            # tv(1, ...) with an intercept: the intercept is the first
            # column
            coef_names[1L] <- if (nzchar(pf)) paste0(pf, ":(Intercept)") else
              "(Intercept)"
            B[, 1L] <- 1
          }
          if (!all(coef_names %in% cn)) next
          part <- if (nzchar(sfx)) "random parameter mean" else
            if (nzchar(pf)) pf else "mean"
          out[[length(out) + 1L]] <- list(
            part = part, term = xl, time = tname, basis = basis,
            periods = as.character(lev), A = B, coef = coef_names)
        }
      }
    }
  }
  c(out, .fcr_tv_smooth_terms(object, bp, cn))
}

# Varying-coefficient smooths s(time, by = x) (from tv(basis = "smooth"))
# and time smooths s(time) as entries for temporal_effects().
.fcr_tv_smooth_terms <- function(object, bp, cn) {
  sm <- bp$smooth
  if (is.null(sm)) return(list())
  blocks <- c(list(main = sm$main), sm$aux %||% list())
  out <- list()
  for (bn in names(blocks)) {
    for (e in blocks[[bn]]) {
      if (length(e$vars) != 1L || is.null(e$sm$Xu)) next
      times <- sort(unique(as.vector(e$sm$Xu[, 1L]) + as.vector(e$sm$shift %||% 0)[1L]))
      if (length(times) > 50L) next
      nd <- stats::setNames(data.frame(times), e$vars)
      if (!is.null(e$by)) nd[[e$by]] <- 1
      P <- mgcv::PredictMat(e$sm, nd)
      sm_names <- paste0(e$label, ".", seq_along(e$cols))
      hit <- cn[endsWith(cn, sm_names[1L])]
      if (!length(hit)) next
      pf <- substr(hit[1L], 1L, nchar(hit[1L]) - nchar(sm_names[1L]))
      coef_names <- paste0(pf, sm_names)
      A <- P
      if (is.null(e$by)) {
        ic <- paste0(pf, "(Intercept)")
        if (ic %in% cn) {
          A <- cbind(1, P)
          coef_names <- c(ic, coef_names)
        }
      }
      if (!all(coef_names %in% cn)) next
      out[[length(out) + 1L]] <- list(
        part = if (nzchar(pf)) sub(":$", "", pf) else "mean", term = e$by %||% "(Intercept)",
        time = e$vars, basis = "smooth", periods = format(times),
        A = A, coef = coef_names)
    }
  }
  out
}

.fcr_tv_select <- function(object, term) {
  tvs <- .fcr_tv_terms(object)
  if (!length(tvs)) {
    stop("The model has no time-varying effects (tv() terms).", call. = FALSE)
  }
  if (!is.null(term)) {
    tvs <- Filter(function(z) z$term %in% term, tvs)
    if (!length(tvs)) {
      stop("No time-varying effect for: ", toString(term), call. = FALSE)
    }
  }
  tvs
}

#' Period effects of time-varying terms
#'
#' Computes the effect of each [tv()] term in each period, with its standard
#' error from the covariance of the coefficients
#' (\eqn{\mathrm{Var}(\hat\beta_t) = A V A^\top}, where \eqn{A} maps the
#' basis coefficients to the period effects), and pointwise confidence
#' intervals. The effects are the same for every basis; only the
#' coefficients differ.
#'
#' @param object a fitted model with [tv()] terms.
#' @param term optional names of the varying variables (`"(Intercept)"` for
#'   `tv(1, ...)`); all by default.
#' @param level confidence level.
#' @param vcov optional covariance matrix of the coefficients, such as a
#'   cluster-robust one; `vcov(object)` by default.
#' @param ... unused.
#'
#' @return A data frame with columns `part` (`"mean"` for the mean, otherwise
#'   the coefficient block, such as `"ln(alpha)"`), `term`, `time`,
#'   `basis`, `period`, `estimate`, `std.error`, `conf.low` and `conf.high`.
#' @seealso [tv()], [temporal_tests()]
#' @export
temporal_effects <- function(object, term = NULL, level = 0.95,
                             vcov = NULL, ...) {
  tvs <- .fcr_tv_select(object, term)
  cf <- stats::coef(object)
  V <- vcov %||% stats::vcov(object)
  z <- stats::qnorm(1 - (1 - level) / 2)
  rows <- lapply(tvs, function(e) {
    b <- cf[e$coef]
    est <- drop(e$A %*% b)
    Vb <- V[e$coef, e$coef, drop = FALSE]
    se <- sqrt(pmax(rowSums((e$A %*% Vb) * e$A), 0))
    data.frame(part = e$part, term = e$term, time = e$time,
               basis = e$basis, period = e$periods, estimate = est,
               std.error = se, conf.low = est - z * se,
               conf.high = est + z * se, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Tests of time-varying effects
#'
#' Wald tests of whether the effect of each [tv()] term changes over time:
#' a joint test that the effects are equal in all periods, and a test of
#' each change between adjacent periods. The tests concern the period
#' effects, so they are the same for every basis.
#'
#' @inheritParams temporal_effects
#' @param adjust p-value adjustment for the adjacent-period tests (see
#'   [stats::p.adjust()]); none by default.
#'
#' @return A data frame with columns `part`, `term`, `test` (`"all equal"`
#'   or the pair of periods), `estimate` (the change, for adjacent
#'   periods), `std.error`, `statistic` (\eqn{\chi^2} for the joint test,
#'   \eqn{z} otherwise), `df` and `p.value`.
#' @seealso [tv()], [temporal_effects()]
#' @export
temporal_tests <- function(object, term = NULL, vcov = NULL,
                           adjust = "none", ...) {
  tvs <- .fcr_tv_select(object, term)
  cf <- stats::coef(object)
  V <- vcov %||% stats::vcov(object)
  rows <- lapply(tvs, function(e) {
    np <- nrow(e$A)
    D <- cbind(-diag(np - 1L), 0) + cbind(0, diag(np - 1L))
    C <- D %*% e$A
    b <- cf[e$coef]
    Vb <- V[e$coef, e$coef, drop = FALSE]
    d <- drop(C %*% b)
    Vd <- C %*% Vb %*% t(C)
    ev <- eigen((Vd + t(Vd)) / 2, symmetric = TRUE)
    keep <- ev$values > max(ev$values) * 1e-10
    w <- sum((crossprod(ev$vectors[, keep, drop = FALSE], d))^2 /
               ev$values[keep])
    dfw <- sum(keep)
    se <- sqrt(pmax(diag(Vd), 0))
    zs <- d / se
    p <- 2 * stats::pnorm(-abs(zs))
    p <- stats::p.adjust(p, method = adjust)
    rbind(
      data.frame(part = e$part, term = e$term, test = "all equal",
                 estimate = NA_real_, std.error = NA_real_, statistic = w,
                 df = dfw, p.value = stats::pchisq(w, dfw, lower.tail = FALSE),
                 stringsAsFactors = FALSE),
      data.frame(part = e$part, term = e$term,
                 test = paste(e$periods[-1L], "vs", e$periods[-np]),
                 estimate = d, std.error = se, statistic = zs, df = 1,
                 p.value = p, stringsAsFactors = FALSE))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}
