# Covariance matrices from the stored scores: Hessian, robust and
# cluster-robust, after the fit; and estfun() / bread() methods for the
# sandwich package.
#
# Every estimator stores the scores of its independent units (rows for
# countreg(), panels for the panel models), multiplied by the weights:
# maxLik's gradientObs, with one row per log-likelihood contribution.

# The maxLik fit holding the scores and the Hessian.
.fcr_ml_fit <- function(object) {
  m <- object$model
  if (!is.null(m$maxLik_object)) m$maxLik_object else m
}

# Unit-level scores (units x parameters), named like coef().
.fcr_scores <- function(object) {
  G <- .fcr_ml_fit(object)$gradientObs
  if (is.null(G)) {
    stop("This model does not store per-unit scores, so sandwich ",
         "covariances are not available for it.", call. = FALSE)
  }
  G <- as.matrix(G)
  colnames(G) <- names(stats::coef(object))
  G
}

# The unit (row of the scores) of each analysis row.
.fcr_unit_index <- function(object, n_units) {
  m <- object$model
  n <- m$nobs %||% stats::nobs(object)
  if (n_units == n) return(seq_len(n))
  pi <- m$panel_index
  if (!is.null(pi) && length(pi) == n && max(pi) == n_units) {
    return(as.integer(pi))
  }
  stop("The rows of this model's scores cannot be matched to its data.",
       call. = FALSE)
}

# Weight of each unit (NULL without weights).
.fcr_unit_weights <- function(object, unit) {
  w <- object$model$weights
  if (is.null(w) || length(w) != length(unit)) return(NULL)
  w[!duplicated(unit)][order(unit[!duplicated(unit)])]
}

# Cluster ids of the analysis rows: a data frame with one column per
# clustering dimension. `cluster` is a one-sided formula (~ a, or ~ a + b
# for two-way clustering), a character vector of column names, or a vector
# (or data frame of vectors) with one value per row of `data` or per
# analysis row.
.fcr_cluster_frame <- function(cluster, data, rows_used) {
  n <- length(rows_used)
  if (inherits(cluster, "formula")) {
    labs <- attr(stats::terms(cluster), "term.labels")
    if (length(labs) == 0L) {
      stop("`cluster` must name at least one variable, as in ~ site.",
           call. = FALSE)
    }
    if (is.null(data)) {
      stop("A `cluster` formula needs the data: the model was fit with ",
           "`keep_data = FALSE`, so give `cluster` as a vector with one ",
           "value per row used in estimation.", call. = FALSE)
    }
    miss <- setdiff(all.vars(cluster), names(data))
    if (length(miss)) {
      stop("Cluster variable(s) not found in the data: ", toString(miss),
           ".", call. = FALSE)
    }
    cf <- lapply(labs, function(l)
      eval(str2lang(l), data, environment(cluster))[rows_used])
    names(cf) <- labs
  } else if (is.character(cluster) && !is.null(data) &&
             all(cluster %in% names(data)) && length(cluster) < n) {
    cf <- lapply(cluster, function(v) data[[v]][rows_used])
    names(cf) <- cluster
  } else {
    cl <- if (is.data.frame(cluster)) as.list(cluster) else list(cluster)
    n_data <- if (!is.null(data)) nrow(data) else NA
    cf <- lapply(cl, function(v) {
      if (length(v) == n) return(v)
      if (!is.na(n_data) && length(v) == n_data) return(v[rows_used])
      stop("`cluster` must be a formula, column name(s), or a vector with ",
           "one value per row of the data (or per row used in estimation, ",
           n, ").", call. = FALSE)
    })
    names(cf) <- if (is.data.frame(cluster)) names(cluster) else "cluster"
  }
  bad <- vapply(cf, anyNA, logical(1))
  if (any(bad)) {
    stop("Cluster variable(s) ", toString(names(cf)[bad]), " have missing ",
         "values in rows used in estimation.", call. = FALSE)
  }
  as.data.frame(cf, stringsAsFactors = FALSE, optional = TRUE)
}

# Cluster-robust covariance (Cameron, Gelbach and Miller 2011) from unit
# scores U (weights included), the inverse of the negative Hessian B, the
# unit of each analysis row and the cluster ids of each analysis row.
# Clusters must contain whole units (a panel may not be split). One-way:
#   V = B M B,  M = G / (G - 1) sum_g S_g S_g',  S_g = sum_{units in g} U_u.
# Multiway: V = sum over non-empty sets S of dimensions of
#   (-1)^(|S| + 1) V_S, with V_S clustered on the intersection of S.
# Units with zero weight are not counted.
.fcr_cluster_vcov <- function(U, B, unit, cf, w_unit = NULL) {
  n_units <- nrow(U)
  keep <- if (is.null(w_unit)) rep(TRUE, n_units) else w_unit > 0
  first <- match(seq_len(n_units), unit)
  ids <- lapply(names(cf), function(v) {
    x <- as.character(cf[[v]])
    xu <- x[first]
    if (any(x != xu[unit])) {
      stop("Cluster variable `", v, "` varies within a panel. Clusters must ",
           "contain whole panels (the model's independent units).",
           call. = FALSE)
    }
    xu
  })
  D <- length(ids)
  V <- matrix(0, ncol(U), ncol(U))
  n_cl <- integer(0)
  for (s in seq_len(D)) {
    for (S in utils::combn(D, s, simplify = FALSE)) {
      g <- do.call(paste, c(ids[S], sep = "\r"))[keep]
      Sg <- rowsum(U[keep, , drop = FALSE], g, reorder = FALSE)
      G <- nrow(Sg)
      adj <- if (G > 1) G / (G - 1) else 1
      V <- V + (-1)^(s + 1) * adj * (B %*% crossprod(Sg) %*% B)
      if (s == 1) n_cl <- c(n_cl, G)
    }
  }
  names(n_cl) <- names(cf)
  dimnames(V) <- list(colnames(U), colnames(U))
  attr(V, "n_clusters") <- n_cl
  V
}

#' Covariance matrix of the estimates
#'
#' \code{vcov()} returns the covariance matrix chosen at fit time, or
#' recomputes a Hessian-based, robust or cluster-robust one from the stored
#' scores, without refitting.
#'
#' @param object a fitted \code{flexCountReg} model.
#' @param type \code{NULL} or \code{"fit"} (default): the covariance chosen
#'   at fit time (\code{stderr}). \code{"hessian"}: the inverse of the
#'   negative Hessian. \code{"robust"}: the sandwich with one independent
#'   unit per row, or per panel for panel models. \code{"cluster"}:
#'   cluster-robust with the clusters in \code{cluster}.
#' @param cluster for \code{type = "cluster"}: a one-sided formula of the
#'   clustering variable(s) in the model's data (\code{~ corridor}; two or
#'   more variables, \code{~ corridor + year}, give multiway clustering),
#'   their names, or a vector with one cluster id per row used in
#'   estimation.
#' @param ... not used.
#'
#' @details Cluster-robust covariances use the scores of the model's
#'   independent units (rows, or panels for \code{countreg_rp()} with
#'   \code{group_var}, \code{renb()}, \code{repois()}, \code{poislind_re()}),
#'   summed within each cluster, with the small-sample factor
#'   \eqn{G/(G-1)} for \eqn{G} clusters. Clusters may be coarser than the
#'   panels (for example sites within corridors) but may not split a panel.
#'   Multiway clustering follows Cameron, Gelbach and Miller (2011). With
#'   frequency weights, a row with weight \eqn{w} stands for \eqn{w}
#'   observations in the same cluster.
#'
#'   The number of clusters is attached as the attribute
#'   \code{"n_clusters"}. The scores and the bread are also available to the
#'   sandwich package through \code{\link{estfun.flexCountReg}}.
#'
#' @references Cameron, A.C., Gelbach, J.B. and Miller, D.L. (2011). Robust
#'   inference with multiway clustering. \emph{Journal of Business and
#'   Economic Statistics}, 29(2), 238-249.
#'
#' @returns A square matrix with the parameter names as dimnames.
#' @examples
#' \donttest{
#' data("washington_roads")
#' m <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
#'               family = "NB2", method = "BFGS")
#' sqrt(diag(vcov(m)))
#' sqrt(diag(vcov(m, type = "cluster", cluster = ~ ID)))
#' }
#' @export
vcov.flexCountReg <- function(object, type = NULL, cluster = NULL, ...) {
  if (is.null(type) || identical(type, "fit")) {
    if (!is.null(cluster)) {
      stop("`cluster` is used with type = \"cluster\".", call. = FALSE)
    }
    return(object$model[["vcov"]])
  }
  type <- match.arg(type, c("hessian", "robust", "cluster"))
  fit <- .fcr_ml_fit(object)
  nm <- names(stats::coef(object))
  if (type == "hessian") {
    V <- .fcr_hessian_vcov(fit)
    dimnames(V) <- list(nm, nm)
    return(V)
  }
  U <- .fcr_scores(object)
  unit <- .fcr_unit_index(object, nrow(U))
  w_unit <- .fcr_unit_weights(object, unit)
  if (type == "robust") {
    if (!is.null(cluster)) {
      stop("`cluster` is used with type = \"cluster\".", call. = FALSE)
    }
    V <- .fcr_robust_vcov(fit, w_unit, object$model$weight_type %||%
                            "frequency")
    dimnames(V) <- list(nm, nm)
    return(V)
  }
  if (is.null(cluster)) {
    stop("type = \"cluster\" needs `cluster`, e.g. cluster = ~ site.",
         call. = FALSE)
  }
  rows <- object$model$rows_used %||% seq_along(unit)
  cf <- .fcr_cluster_frame(cluster, object$data, rows)
  .fcr_cluster_vcov(U, .fcr_hessian_vcov(fit), unit, cf, w_unit)
}

#' Scores and bread for the sandwich package
#'
#' Methods that make \code{flexCountReg} fits work with
#' \code{sandwich::vcovCL()}, \code{sandwich::vcovHC()},
#' \code{sandwich::sandwich()} and, through them,
#' \code{lmtest::coeftest()} and \code{modelsummary(vcov = ~ cluster)}.
#'
#' @param x a fitted \code{flexCountReg} model.
#' @param ... not used.
#'
#' @details \code{estfun()} returns the scores of the model's independent
#'   units: one row per row of the data used in estimation, or per panel for
#'   panel models, multiplied by the weights. \code{bread()} is
#'   \eqn{n (-H)^{-1}} with \eqn{n} the number of units, the convention of
#'   the sandwich package, so \code{sandwich::sandwich(m)} is
#'   \eqn{(-H)^{-1} \sum_i s_i s_i' (-H)^{-1}}.
#'
#'   The sandwich package treats each row of \code{estfun()} as an
#'   independent unit. That is right without weights and for sampling
#'   weights; with frequency weights, use \code{vcov(m, type = "robust")}
#'   or \code{type = "cluster"}, which count a row with weight \eqn{w} as
#'   \eqn{w} observations. For panel models, a cluster variable for
#'   \code{vcovCL()} must have one value per panel (in the order of the
#'   panel ids); \code{vcov(m, type = "cluster", cluster = ~ v)} takes the
#'   variable from the data.
#'
#' @returns \code{estfun()}: a matrix with one row per unit and one column
#'   per parameter. \code{bread()}: a square matrix.
#' @examples
#' \donttest{
#' data("washington_roads")
#' m <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
#'               family = "NB2", method = "BFGS")
#' sandwich::vcovCL(m, cluster = washington_roads$ID, type = "HC0")
#' }
#' @importFrom sandwich estfun bread
#' @export
estfun.flexCountReg <- function(x, ...) {
  .fcr_scores(x)
}

#' @rdname estfun.flexCountReg
#' @export
bread.flexCountReg <- function(x, ...) {
  B <- .fcr_hessian_vcov(.fcr_ml_fit(x))
  nm <- names(stats::coef(x))
  dimnames(B) <- list(nm, nm)
  nrow(.fcr_scores(x)) * B
}
