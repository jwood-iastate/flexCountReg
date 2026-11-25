# Fixed-effects count regression: countreg_fe().
#
# Three estimators, by structure:
#   "conditional"  Poisson, one fixed-effect dimension: the conditional
#                  (multinomial) likelihood, which does not contain the
#                  group effects. Newton-Raphson with closed-form score and
#                  Hessian. Identical slopes to the dummy-variable Poisson.
#   "demean"       Poisson, any number of dimensions: IRLS in which each
#                  weighted least-squares step demeans the working response
#                  and the regressors by alternating projections on the
#                  fixed-effect dimensions (Correia, Guimaraes and Zylkin
#                  2020; the method of ppmlhdfe and fixest::fepois()).
#   "concentrated" any family, one dimension: the dummy-variable likelihood
#                  with each group's intercept solved by a one-dimensional
#                  Newton step for given slopes and distribution parameters;
#                  the profile score follows by the envelope theorem.
#
# Groups whose counts are all zero carry no information on the slopes (their
# effect is -Inf) and are dropped, as are single-observation groups
# (perfectly fitted by their own effect); for several dimensions both are
# dropped repeatedly until none remain.

#' Fixed-effects count regression
#'
#' Estimates count regression models with fixed effects for one or more
#' grouping variables (sites, years, ...). Poisson models use the
#' conditional likelihood (one fixed-effect dimension) or IRLS with
#' alternating projections (several dimensions); other families, including
#' families made with \code{\link{new_count_family}()}, use the
#' dummy-variable likelihood with the group effects concentrated out.
#'
#' @param formula model formula of the mean, \code{y ~ x1 + x2}, with an
#'   optional \code{offset()}. The intercept is absorbed by the fixed
#'   effects.
#' @param fixef a one-sided formula of the fixed-effect variables, e.g.
#'   \code{~ site} or \code{~ site + year}.
#' @param data a data frame.
#' @param family the count family: \code{"POISSON"} (default) or any other
#'   family of \code{\link{countreg}()} (one fixed-effect dimension; enter
#'   further dimensions with few levels, such as years, as factors in
#'   \code{formula}), or a family made with \code{new_count_family()}.
#' @param offset optional name of an offset column (or use \code{offset()}
#'   in the formula).
#' @param weights optional name of a column of frequency weights.
#' @param method \code{"auto"} (default: \code{"conditional"} for Poisson
#'   with one dimension, \code{"demean"} for Poisson with several,
#'   \code{"concentrated"} otherwise), or one of those.
#' @param stderr \code{"cluster"} (default; cluster-robust by
#'   \code{cluster}) or \code{"hessian"}.
#' @param cluster the clustering variable(s) for \code{stderr = "cluster"}:
#'   a one-sided formula, column name(s) or a vector (see
#'   \code{\link{vcov.flexCountReg}}). The default is the first fixed-effect
#'   variable. Clusters must contain whole groups of that variable.
#' @param jackknife logical: apply the split-panel jackknife bias correction
#'   of Dhaene and Jochmans (2015) to the slopes (for families other than
#'   Poisson in short panels). Each group's observations, ordered by
#'   \code{time}, are split into a first and second half (sharing the
#'   middle period when the number of periods is odd); the corrected
#'   estimate is \eqn{2\hat\theta - (\hat\theta_1 + \hat\theta_2)/2}, for
#'   the slopes and the distribution parameters.
#' @param time optional column that orders the observations within groups
#'   for \code{jackknife} (default: the row order).
#' @param max.iters maximum number of iterations.
#' @param tol convergence tolerance (relative change of the log-likelihood
#'   or of the deviance).
#' @param keep_data logical; store the data in the fit (needed for
#'   predictions on the training rows without \code{newdata}).
#' @param sp optional fixed smoothing parameters for smooth terms in
#'   \code{formula}; by default they are selected automatically.
#'
#' @details Poisson fixed-effects estimates of the slopes are consistent
#'   with fixed T and need only a correct conditional mean (Wooldridge
#'   1999). For other families the concentrated (dummy-variable) estimates
#'   have an incidental-parameter bias of order 1/T. For NB2 the bias of the
#'   slopes is small (Allison and Waterman 2002), but the dispersion is
#'   strongly underestimated in short panels. In simulations with 200 groups
#'   (NB2, size 2), the bias of ln(alpha) was -1.08, -0.55 and -0.24 with
#'   T = 3, 5 and 10 periods, and slope biases were 1-2\%.
#'   \code{jackknife = TRUE} reduced the dispersion bias to -0.03 (T = 5) and
#'   0.06 (T = 10) but did not improve the slopes, and it failed with T = 3
#'   (its half-panels are too short; a warning is given below 4 periods). The conditional negative binomial of Hausman, Hall and
#'   Griliches is a different model (its fixed effect enters the
#'   dispersion) and is not provided.
#'
#'   Standard errors are cluster-robust by the first fixed-effect variable
#'   by default, with the factor \eqn{G/(G-1)} for \eqn{G} clusters. They
#'   equal those of \code{fixest::fepois()} and \code{fixest::fenegbin()}
#'   when the clustering is given explicitly, \code{vcov = ~ site} (fixest
#'   0.13 and later default to iid standard errors), with
#'   \code{ssc(adj = FALSE)}; by default fixest also multiplies the
#'   variance by \eqn{(n-1)/(n-K)}. \code{vcov(fit, type = "robust")}
#'   treats each group as an independent unit with the conditional and
#'   concentrated likelihoods (it is clustered by group, without
#'   \eqn{G/(G-1)}), and each row with the \code{"demean"} method. The log-likelihood is that of the dummy-variable model,
#'   and its degrees of freedom count the fixed effects, so AIC and BIC
#'   compare with models that estimate them.
#'
#'   \strong{Smooth terms.} \code{formula} may contain smooth terms
#'   (\code{s(x)}, \code{te()}, ...) of variables that change within
#'   groups; a smooth of a variable that is constant within groups is
#'   absorbed by the fixed effects and gives an error. The fit is then
#'   penalized, with the smoothing parameters selected as in
#'   \code{\link{countreg}()}. The smooth bases are built on the rows used in
#'   estimation (after dropping all-zero and single-observation groups), so
#'   with fixed smoothing parameters the estimates equal those of
#'   \code{mgcv::gam()} with group dummies. The jackknife is not available
#'   with smooth terms.
#'
#'   \code{predict()} uses each known group's estimated effect. A group not
#'   in the estimation data has no fixed effect: its prediction is
#'   \code{NA}, unless \code{newgroup = "mean_fe"} is given to use the
#'   mean of the estimated effects, weighted by the groups' sizes in the
#'   estimation sample (their numbers of rows, counted with the frequency
#'   weights). This mean is fixed at estimation, so the prediction for a new
#'   group does not depend on the other rows of \code{newdata}.
#'
#' @references
#' Allison, P.D. and Waterman, R.P. (2002). Fixed-effects negative binomial
#' regression models. \emph{Sociological Methodology}, 32, 247-265.
#'
#' Correia, S., Guimaraes, P. and Zylkin, T. (2020). Fast Poisson estimation
#' with high-dimensional fixed effects. \emph{The Stata Journal}, 20(1),
#' 95-115.
#'
#' Dhaene, G. and Jochmans, K. (2015). Split-panel jackknife estimation of
#' fixed-effect models. \emph{Review of Economic Studies}, 82(3), 991-1030.
#'
#' Wooldridge, J.M. (1999). Distribution-free estimation of some nonlinear
#' panel data models. \emph{Journal of Econometrics}, 90(1), 77-97.
#'
#' @returns An object of class \code{c("countreg_fe", "flexCountReg")}.
#'   \code{model$fixef} holds the estimated effects of each dimension,
#'   \code{model$dropped} the groups and rows that were dropped, and
#'   \code{model$not_identified} the regressors without variation within
#'   groups.
#' @examples
#' \donttest{
#' data("washington_roads")
#' fe <- countreg_fe(Total_crashes ~ lnaadt + offset(lnlength),
#'                   fixef = ~ ID, data = washington_roads)
#' summary(fe)
#' }
#' @export
countreg_fe <- function(formula, fixef, data, family = "POISSON",
                        offset = NULL, weights = NULL,
                        method = c("auto", "conditional", "demean",
                                   "concentrated"),
                        stderr = c("cluster", "hessian"), cluster = NULL,
                        jackknife = FALSE, time = NULL,
                        max.iters = 1000, tol = 1e-10, keep_data = TRUE,
                        sp = NULL) {
  method <- match.arg(method)
  stderr <- match.arg(stderr)
  family_spec <- if (inherits(family, "fcr_count_family")) family
  family <- .fcr_resolve_family(family)
  if (!inherits(fixef, "formula") || length(fixef) != 2L) {
    stop("`fixef` must be a one-sided formula, e.g. ~ site or ~ site + year.",
         call. = FALSE)
  }
  fe_vars <- attr(stats::terms(fixef), "term.labels")
  if (!length(fe_vars) || !all(fe_vars %in% all.vars(fixef))) {
    stop("`fixef` must list variables, e.g. ~ site + year.", call. = FALSE)
  }
  data <- as.data.frame(data)
  miss <- setdiff(fe_vars, names(data))
  if (length(miss)) {
    stop("Fixed-effect variable(s) not found in `data`: ", toString(miss),
         ".", call. = FALSE)
  }
  if (method == "auto") {
    method <- if (family != "POISSON") "concentrated" else
      if (length(fe_vars) == 1L) "conditional" else "demean"
  }
  if (method %in% c("conditional", "demean") && family != "POISSON") {
    stop("method = \"", method, "\" is for the Poisson family; use ",
         "\"concentrated\".", call. = FALSE)
  }
  if (method %in% c("conditional", "concentrated") && length(fe_vars) > 1L) {
    stop("method = \"", method, "\" takes one fixed-effect dimension. ",
         if (family == "POISSON") "Use method = \"demean\", or enter " else
           "Enter ",
         "further dimensions with few levels (such as years) as factors ",
         "in `formula`.", call. = FALSE)
  }
  .fcr_check_controls(max.iters = max.iters)

  # --- Data -----------------------------------------------------------------
  prep <- prepare_model_data(formula, data, offset = offset,
                             weights = weights, group_var = fe_vars)
  formula <- prep$formula
  if (prep$n_dropped > 0) {
    message(prep$n_dropped, " row(s) with missing values were removed.")
  }
  X <- prep$X[, colnames(prep$X) != "(Intercept)", drop = FALSE]
  y <- prep$y
  off <- prep$offset
  w <- prep$weights
  rows <- prep$rows_used
  fe <- lapply(fe_vars, function(v) as.character(data[[v]][rows]))
  names(fe) <- fe_vars

  # Drop all-zero and single-observation groups (repeatedly, across
  # dimensions), which carry no information on the slopes
  keep <- w > 0
  dropped <- list()
  repeat {
    drop_now <- rep(FALSE, length(y))
    for (v in fe_vars) {
      g <- fe[[v]]
      tot <- tapply(w * y * keep, g, sum)
      n_g <- tapply(keep, g, sum)
      # Observations of a group, counting frequency weights (a row with
      # weight 3 is three observations, as in the expanded data)
      n_obs <- tapply(w * keep, g, sum)
      zero <- names(tot)[tot == 0 & n_g > 0]
      single <- if (family == "POISSON" || method == "concentrated")
        names(n_obs)[n_g > 0 & n_obs <= 1]
      bad <- union(zero, single)
      if (length(bad)) {
        hit <- keep & g %in% bad
        drop_now <- drop_now | hit
        dropped[[v]] <- union(dropped[[v]], bad[bad %in% g[hit]])
      }
    }
    if (!any(drop_now)) break
    keep <- keep & !drop_now
  }
  n_drop_rows <- sum(w > 0 & !keep)
  if (n_drop_rows > 0) {
    message(n_drop_rows, " row(s) in ", sum(lengths(dropped)), " group(s) ",
            "with all-zero counts or a single observation were dropped: ",
            "they carry no information on the slopes.")
  }
  if (!any(keep)) stop("No groups with positive counts remain.", call. = FALSE)
  idx <- which(keep)
  # Smooth bases are rebuilt on the rows used in estimation (the dropped
  # groups do not enter the fit), so they span the data actually fitted;
  # the blueprint keeps the rebuilt smooths for prediction
  if (length(prep$smooth$main)) {
    sm_specs <- .fcr_smooth_split(formula)$specs
    built <- .fcr_smooth_build(sm_specs,
                               droplevels(data[rows[idx], , drop = FALSE]))
    old_cols <- unlist(lapply(prep$smooth$main, `[[`, "cols"))
    p0 <- min(old_cols) - 1L
    newB <- matrix(NA_real_, nrow(X), ncol(built$B),
                   dimnames = list(NULL, colnames(built$B)))
    newB[idx, ] <- built$B
    Xfull <- prep$X
    Xfull <- cbind(Xfull[, -old_cols, drop = FALSE], newB)
    prep$smooth$main <- .fcr_smooth_shift(built$smooths, p0)
    prep$blueprint$smooth$main <- prep$smooth$main
    prep$X <- Xfull
    X <- Xfull[, colnames(Xfull) != "(Intercept)", drop = FALSE]
  }
  X <- X[idx, , drop = FALSE]
  y <- y[idx]
  off <- off[idx]
  w <- w[idx]
  fe <- lapply(fe, function(g) factor(g[idx]))
  rows_fit <- rows[idx]

  # Regressors without variation within groups are not identified
  within_sd <- function(x) {
    xt <- x
    for (v in fe_vars) {
      for (it in 1:50) {
        m <- stats::ave(xt, fe[[v]])
        xt <- xt - m
        if (max(abs(m)) < 1e-12 * max(1, abs(x))) break
      }
    }
    max(abs(xt))
  }
  flat <- vapply(seq_len(ncol(X)), function(j)
    within_sd(X[, j]) < 1e-9 * max(1, abs(X[, j])), logical(1))
  # Smooth terms: a smooth of a variable that is constant within groups is
  # absorbed by the fixed effects; its basis columns are never dropped one
  # by one
  sm_ents <- prep$smooth$main
  sm_cols <- lapply(sm_ents, function(e) colnames(prep$X)[e$cols])
  for (k in seq_along(sm_ents)) {
    fl <- flat[colnames(X) %in% sm_cols[[k]]]
    if (any(fl)) {
      stop("The smooth ", sm_ents[[k]]$label, " has no variation within ",
           "the fixed-effect groups, so it is absorbed by the fixed effects. ",
           "Use smooths of variables that change within groups.",
           call. = FALSE)
    }
  }
  not_identified <- colnames(X)[flat]
  if (any(flat)) {
    warning("Regressor(s) without variation within the fixed-effect groups ",
            "are not identified and were removed: ",
            toString(not_identified), ".", call. = FALSE)
    X <- X[, !flat, drop = FALSE]
  }
  if (ncol(X) == 0L) {
    stop("No regressor varies within the fixed-effect groups.",
         call. = FALSE)
  }
  smooth <- NULL
  if (length(sm_ents)) {
    if (isTRUE(jackknife)) {
      stop("The split-panel jackknife is not available with smooth terms.",
           call. = FALSE)
    }
    for (k in seq_along(sm_ents)) {
      sm_ents[[k]]$cols <- match(sm_cols[[k]], colnames(X))
    }
    smooth <- list(ents = list(main = sm_ents), sp = sp)
  }

  fit <- switch(method,
    conditional = .fcr_fe_conditional(y, X, off, w, fe[[1]], max.iters, tol,
                                      smooth = smooth),
    demean = .fcr_fe_demean(y, X, off, w, fe, max.iters, tol,
                            smooth = smooth),
    concentrated = .fcr_fe_concentrated(y, X, off, w, fe[[1]], family,
                                        max.iters, tol, smooth = smooth))

  # --- Split-panel jackknife -------------------------------------------------
  jk <- NULL
  if (isTRUE(jackknife)) {
    ord <- if (is.null(time)) seq_along(y) else {
      if (!time %in% names(data)) {
        stop("`time` \"", time, "\" is not a column of `data`.",
             call. = FALSE)
      }
      data[[time]][rows_fit]
    }
    g1 <- fe[[1]]
    rk <- stats::ave(seq_along(y), g1, FUN = function(i) rank(ord[i],
                                                              ties.method = "first"))
    Tg <- stats::ave(rep(1, length(y)), g1, FUN = length)
    if (stats::median(Tg) < 4) {
      warning("The split-panel jackknife needs panels of at least 4 ",
              "periods: with ", stats::median(Tg), " (median), the half-",
              "panels are too short and the correction is unreliable.",
              call. = FALSE)
    }
    # Halves of each group's panel; with an odd number of periods the two
    # halves share the middle period (Dhaene and Jochmans 2015, Sec. 5)
    in_half <- list(rk <= ceiling(Tg / 2), rk > floor(Tg / 2))
    fit_half <- function(h) {
      s <- in_half[[h]]
      cnt <- tapply(w[s] * y[s], factor(g1[s]), sum)
      ok <- s & as.character(g1) %in% names(cnt)[cnt > 0] &
        stats::ave(as.numeric(s), g1, FUN = sum) > 1
      sel <- which(ok)
      f <- switch(method,
        conditional = .fcr_fe_conditional(y[sel], X[sel, , drop = FALSE],
                                          off[sel], w[sel],
                                          factor(g1[sel]), max.iters, tol),
        demean = .fcr_fe_demean(y[sel], X[sel, , drop = FALSE], off[sel],
                                w[sel], lapply(fe, function(g)
                                  factor(g[sel])), max.iters, tol),
        concentrated = .fcr_fe_concentrated(y[sel], X[sel, , drop = FALSE],
                                            off[sel], w[sel],
                                            factor(g1[sel]), family,
                                            max.iters, tol))
      f$estimate
    }
    b1 <- fit_half(1L)
    b2 <- fit_half(2L)
    jk <- list(uncorrected = fit$estimate, half1 = b1, half2 = b2)
    fit$estimate <- 2 * fit$estimate - (b1 + b2) / 2
  }

  # --- Standard errors -----------------------------------------------------
  x_names <- names(fit$estimate)
  B <- .fcr_hessian_vcov(fit)
  dimnames(B) <- list(x_names, x_names)
  # Units of the scores: groups of the first dimension (conditional and
  # concentrated) or rows (demean)
  unit <- if (method == "demean") seq_along(y) else as.integer(fe[[1]])
  if (stderr == "cluster") {
    if (is.null(cluster)) cluster <- stats::reformulate(fe_vars[1])
    cf <- .fcr_cluster_frame(cluster, data, rows_fit)
    U <- fit$gradientObs
    colnames(U) <- x_names
    V <- .fcr_cluster_vcov(U, B, unit, cf)
    n_clusters <- attr(V, "n_clusters")
    attr(V, "n_clusters") <- NULL
  } else {
    V <- B
    n_clusters <- NULL
  }

  # --- Result ----------------------------------------------------------------
  label <- switch(method,
    conditional = "Poisson fixed effects (conditional likelihood)",
    demean = "Poisson fixed effects (IRLS with alternating projections)",
    concentrated = paste0(.fcr_family_label(family), " fixed effects ",
                          "(concentrated dummy-variable likelihood)"))
  n_fe <- sum(vapply(fit$fixef, length, 0L)) - (length(fe_vars) - 1L)
  model <- list(
    estimate = fit$estimate, maximum = fit$loglik, hessian = fit$hessian,
    hessian_accuracy = fit$hessian_accuracy,
    smooth_table = fit$smooth_table, smooth_info = fit$smooth_info,
    smooth_lambda = fit$smooth_lambda,
    smooth_converged = fit$smooth_converged,
    smooth_iterations = fit$smooth_iterations,
    smooth_coefs = fit$smooth_coefs, edf = fit$edf, df_total = fit$df_total,
    gradientObs = fit$gradientObs, vcov = V, se = .fcr_se(V),
    vcov_type = if (stderr == "cluster") "cluster" else
      if (!is.null(smooth)) "bayesian" else "hessian",
    stderr = stderr, n_clusters = n_clusters,
    cluster = if (stderr == "cluster") names(cf),
    iterations = fit$iterations, converged = fit$converged,
    message = if (isTRUE(fit$converged)) "converged" else
      "iteration limit reached",
    code = if (isTRUE(fit$converged)) 0L else 1L,
    family = family, family_spec = family_spec, method = method,
    fe_method_label = label, modelType = "countreg_fe",
    formula = formula, fixef_formula = fixef, fixef_vars = fe_vars,
    fixef = fit$fixef, n_fe = n_fe, dropped = dropped,
    # Size of each group (sum of the frequency weights of its rows) in the
    # estimation sample, in the order of `fixef`
    fixef_n = stats::setNames(lapply(seq_along(fe), function(d)
      as.vector(rowsum(w, as.integer(fe[[d]]), reorder = TRUE))),
      names(fit$fixef)),
    not_identified = not_identified, jackknife = jk,
    blueprint = prep$blueprint, offset = prep$blueprint$offset_label,
    y = y, nobs = length(y), n_groups = nlevels(fe[[1]]),
    panel_index = if (method != "demean") as.integer(fe[[1]]),
    rows_used = rows_fit, weights = w, weight_type = "frequency",
    n_bic = sum(w))
  names(model$se) <- x_names
  if (!isTRUE(fit$converged)) {
    warning("The fixed-effects estimation did not converge in ", max.iters,
            " iterations.", call. = FALSE)
  }
  .createFlexCountReg(model = model, data = data, keep_data = keep_data,
                      call = match.call(), formula = formula)
}

# ---- Poisson: conditional likelihood (one dimension) -------------------------
# Group i, rows t: l_i = sum_t w y eta - S_i log sum_t w exp(eta) (+ const),
# S_i = sum_t w y. With p = w exp(eta) / sum w exp(eta) within the group,
#   score_i = sum_t w (y - S_i p / w) x = sum_t (w y - S_i p) x,
#   Hessian = -sum_i S_i sum_t p (x - xbar_i)(x - xbar_i)',
#   xbar_i = sum_t p x.
.fcr_fe_conditional <- function(y, X, off, w, g, max.iters, tol,
                                smooth = NULL) {
  gi <- as.integer(g)
  S <- as.vector(rowsum(w * y, gi))
  wy <- w * y
  beta <- rep(0, ncol(X))
  names(beta) <- colnames(X)
  parts <- function(beta) {
    eta <- as.vector(X %*% beta) + off
    m <- as.vector(tapply(eta, gi, max))[gi]
    e <- w * exp(eta - m)
    den <- as.vector(rowsum(e, gi))
    p <- e / den[gi]
    xbar <- rowsum(p * X, gi)
    Xc <- X - xbar[gi, , drop = FALSE]
    Ug <- rowsum((wy - S[gi] * p) * X, gi)
    H <- -crossprod(Xc * sqrt(S[gi] * p))
    ll <- sum(wy * eta) - sum(S * (log(den) + as.vector(tapply(eta, gi, max))))
    list(ll = ll, U = Ug, H = H, eta = eta)
  }
  sf <- NULL
  if (!is.null(smooth)) {
    # Penalized conditional likelihood: per-group contributions and scores
    fn <- function(b) {
      eta <- as.vector(X %*% b) + off
      mg <- as.vector(tapply(eta, gi, max))
      den <- as.vector(rowsum(w * exp(eta - mg[gi]), gi))
      as.vector(rowsum(wy * eta, gi)) - S * (log(den) + mg)
    }
    gr <- function(b) unname(parts(b)$U)
    sf <- .fcr_smooth_fit(smooth$ents, fn, gr, beta, "BFGS", max.iters,
                          smooth$sp, list(main = 0), list(main = "mean"),
                          list(main = X))
    beta <- sf$estimate
    max.iters <- 0L
  }
  pr <- parts(beta)
  conv <- !is.null(sf) && sf$code == 0L
  it <- 0L
  while (it < max.iters) {
    it <- it + 1L
    step <- tryCatch(solve(-pr$H, colSums(pr$U)), error = function(e) NULL)
    if (is.null(step)) break
    s <- 1
    repeat {
      nb <- beta + s * step
      np <- parts(nb)
      if (is.finite(np$ll) && np$ll >= pr$ll - 1e-12 * abs(pr$ll)) break
      s <- s / 2
      if (s < 1e-10) break
    }
    done <- abs(np$ll - pr$ll) <= tol * max(1, abs(pr$ll)) &&
      max(abs(s * step)) < 1e-8
    beta <- nb
    pr <- np
    if (done) {
      conv <- TRUE
      break
    }
  }
  # Group effects and the dummy-variable log-likelihood
  alpha <- log(S) - as.vector(log(rowsum(w * exp(pr$eta), gi)))
  mu <- exp(alpha[gi] + pr$eta)
  ll_full <- sum(w * stats::dpois(y, mu, log = TRUE))
  names(alpha) <- levels(g)
  H <- if (!is.null(sf)) sf$hessian else pr$H
  dimnames(H) <- list(names(beta), names(beta))
  .fcr_fe_with_smooth(
    list(estimate = beta, loglik = ll_full, hessian = H, gradientObs = pr$U,
         fixef = stats::setNames(list(alpha), "fe1"),
         iterations = if (!is.null(sf)) sf$iterations else it,
         converged = conv), sf)
}

# ---- Poisson: IRLS with alternating projections (any dimensions) -------------
# Weighted demeaning of the columns of M on the fixed-effect factors.
.fcr_fe_demean_cols <- function(M, fe, W, tol = 1e-12, max_it = 10000L) {
  M <- as.matrix(M)
  for (it in seq_len(max_it)) {
    M_old <- M
    for (g in fe) {
      gi <- as.integer(g)
      sw <- as.vector(rowsum(W, gi))
      mean_g <- rowsum(W * M, gi) / sw
      M <- M - mean_g[gi, , drop = FALSE]
    }
    if (max(abs(M - M_old)) < tol * max(1, max(abs(M)))) break
  }
  M
}

.fcr_fe_demean <- function(y, X, off, w, fe, max.iters, tol, smooth = NULL) {
  # With smooth terms: penalized IRLS, each weighted least-squares step
  # solving (X'WX + S_lambda) b = X'Wz on the demeaned variables, inside
  # Fellner-Schall updates of the smoothing parameters
  P <- ncol(X)
  terms <- if (!is.null(smooth))
    .fcr_pen_terms(smooth$ents, list(main = 0), list(main = "mean"))
  lambda <- if (is.null(smooth)) numeric(0) else
    smooth$sp %||% rep(1, length(terms))
  Sl <- if (is.null(smooth)) matrix(0, P, P) else
    .fcr_pen_matrix(terms, lambda, P)
  fs_it <- 0L
  fs_conv <- is.null(smooth) || !is.null(smooth$sp)
  repeat {
  mu <- (y + mean(y)) / 2
  eta <- log(mu)
  dev_old <- Inf
  conv <- FALSE
  it <- 0L
  beta <- stats::setNames(rep(0, ncol(X)), colnames(X))
  while (it < max.iters) {
    it <- it + 1L
    W <- w * mu
    z <- eta - off + (y - mu) / mu
    Zt <- .fcr_fe_demean_cols(cbind(z, X), fe, W)
    zt <- Zt[, 1]
    Xt <- Zt[, -1, drop = FALSE]
    beta <- stats::setNames(
      as.vector(solve(crossprod(Xt * sqrt(W)) + Sl,
                      crossprod(Xt * W, zt))),
      colnames(X))
    e <- zt - as.vector(Xt %*% beta)
    eta <- z - e + off
    mu <- exp(eta)
    dev <- 2 * sum(w * ifelse(y > 0, y * log(y / mu), 0) - w * (y - mu))
    if (abs(dev - dev_old) <= tol * max(1, abs(dev))) {
      conv <- TRUE
      break
    }
    dev_old <- dev
  }
  W <- w * mu
  Xt <- .fcr_fe_demean_cols(X, fe, W)
  Hn <- crossprod(Xt * sqrt(W))
  if (fs_conv) break
  fs_it <- fs_it + 1L
  new <- .fcr_fs_step(beta, Hn, terms, lambda)
  done <- max(abs(log(new) - log(lambda))) < 1e-3 || fs_it >= 200L
  lambda <- new
  Sl <- .fcr_pen_matrix(terms, lambda, P)
  if (done) fs_conv <- TRUE
  }
  H <- -(Hn + Sl)
  dimnames(H) <- list(names(beta), names(beta))
  U <- (w * (y - mu)) * Xt
  # Effects of each dimension from the total effect eta - X beta - offset
  # (the last dimensions' effects have weighted mean zero)
  r <- eta - as.vector(X %*% beta) - off
  alphas <- lapply(fe, function(g) stats::setNames(rep(0, nlevels(g)),
                                                    levels(g)))
  for (k in seq_len(1000)) {
    old <- unlist(alphas)
    for (d in seq_along(fe)) {
      gi <- as.integer(fe[[d]])
      other <- Reduce(`+`, lapply(seq_along(fe)[-d], function(e)
        alphas[[e]][as.integer(fe[[e]])]), 0)
      alphas[[d]][] <- as.vector(rowsum(w * (r - other), gi) /
                                   rowsum(w, gi))
    }
    if (max(abs(unlist(alphas) - old)) < 1e-12) break
  }
  if (length(fe) > 1L) {
    for (d in seq_along(fe)[-1]) {
      gi <- as.integer(fe[[d]])
      c_d <- sum(w * alphas[[d]][gi]) / sum(w)
      alphas[[d]] <- alphas[[d]] - c_d
      alphas[[1]] <- alphas[[1]] + c_d
    }
  }
  names(alphas) <- paste0("fe", seq_along(fe))
  out <- list(estimate = beta,
              loglik = sum(w * stats::dpois(y, mu, log = TRUE)),
              hessian = H, gradientObs = U, fixef = alphas, iterations = it,
              converged = conv)
  if (is.null(smooth)) return(out)
  Vb <- tryCatch(solve(Hn + Sl), error = function(e) MASS::ginv(Hn + Sl))
  dimnames(Vb) <- dimnames(H)
  gam <- list(theta = beta, Vb = Vb, edf = stats::setNames(diag(Vb %*% Hn),
                                                          names(beta)),
              lambda = lambda, terms = terms, converged = TRUE,
              iterations = fs_it)
  sf <- .fcr_smooth_outputs(gam, smooth$ents, list(main = X))
  .fcr_fe_with_smooth(out, list(
    smooth_table = sf$table, smooth_info = sf$info, smooth_lambda = lambda,
    smooth_converged = fs_it < 200L, smooth_iterations = fs_it,
    smooth_coefs = unique(unlist(lapply(terms, `[[`, "idx"))),
    edf = gam$edf, df_total = sum(gam$edf)))
}

# ---- Any family: concentrated dummy-variable likelihood -----------------------
.fcr_fe_concentrated <- function(y, X, off, w, g, family, max.iters, tol,
                                 smooth = NULL) {
  params <- get_params(family)
  n1 <- !is.null(params[[1]])
  n2 <- !is.null(params[[2]])
  gi <- as.integer(g)
  G <- nlevels(g)
  dens <- get_probFunc(family)
  dens_args <- list(haltons = NULL, normed_haltons = NULL,
                    engine = "quadrature")
  fam <- .fcr_get_family(family)
  # Starting values: Poisson conditional slopes, registry dispersion starts
  pois <- .fcr_fe_conditional(y, X, off, w, g, 100L, 1e-8)
  st <- fam$start %||% 0
  start <- c(pois$estimate,
             if (n1) stats::setNames(rep_len(st, 2)[1], params[[1]]),
             if (n2) stats::setNames(rep_len(st, 2)[2], params[[2]]))
  kx <- ncol(X)
  alpha <- pois$fixef[[1]]
  logf <- function(mu, a1, a2) {
    tryCatch(do.call(dens, c(list(y = y, predicted = mu, alpha = a1,
                                  sigma = a2), dens_args, list(log = TRUE))),
             error = function(e) rep(NA_real_, length(y)))
  }
  score_eta <- function(mu, a1, a2) {
    .fcr_family_score(family, y, mu, a1, a2, dens_args = dens_args)
  }
  # Group effects for given parameters: Newton on sum_t w g_t(alpha) = 0
  solve_alpha <- function(eta0, a1, a2, alpha) {
    for (k in 1:100) {
      mu <- exp(alpha[gi] + eta0)
      s <- score_eta(mu, a1, a2)$eta
      Gs <- as.vector(rowsum(w * s, gi))
      h <- 1e-4
      sp <- score_eta(mu * exp(h), a1, a2)$eta
      sm <- score_eta(mu * exp(-h), a1, a2)$eta
      D <- as.vector(rowsum(w * (sp - sm) / (2 * h), gi))
      if (any(!is.finite(Gs)) || any(!is.finite(D))) return(NULL)
      step <- -Gs / pmin(D, -1e-10)
      step <- pmax(pmin(step, 2), -2)
      alpha <- alpha + step
      if (max(abs(step)) < 1e-10) break
    }
    alpha
  }
  cache <- new.env()
  cache$alpha <- alpha
  state <- function(p) {
    if (identical(cache$p, p)) return(cache$st)
    beta <- p[seq_len(kx)]
    a1 <- if (n1) rep(exp(p[[kx + 1]]), length(y))
    a2 <- if (n2) rep(exp(p[[kx + 1 + n1]]), length(y))
    eta0 <- as.vector(X %*% beta) + off
    al <- solve_alpha(eta0, a1, a2, cache$alpha)
    st <- if (is.null(al)) NULL else {
      cache$alpha <- al
      mu <- exp(al[gi] + eta0)
      lf <- logf(mu, a1, a2)
      sc <- score_eta(mu, a1, a2)
      list(ll = as.vector(rowsum(w * lf, gi)), mu = mu, sc = sc, alpha = al)
    }
    cache$p <- p
    cache$st <- st
    st
  }
  fn <- function(p) {
    s <- state(p)
    if (is.null(s) || any(!is.finite(s$ll))) return(rep(-Inf, G))
    s$ll
  }
  gr <- function(p) {
    s <- state(p)
    if (is.null(s)) return(matrix(NA_real_, G, length(p)))
    cols <- list(rowsum(w * s$sc$eta * X, gi))
    if (n1) cols <- c(cols, list(rowsum(w * s$sc$alpha, gi)))
    if (n2) cols <- c(cols, list(rowsum(w * s$sc$sigma, gi)))
    out <- do.call(cbind, cols)
    dimnames(out) <- NULL
    out
  }
  sf <- NULL
  if (!is.null(smooth)) {
    ml <- sf <- .fcr_smooth_fit(smooth$ents, fn, gr, start, "BFGS",
                                max.iters, smooth$sp, list(main = 0),
                                list(main = "mean"), list(main = X))
  } else {
    ml <- maxLik::maxLik(fn, grad = gr, start = start, method = "BFGS",
                         finalHessian = FALSE,
                         control = list(iterlim = max.iters,
                                        reltol = max(tol, 1e-12)))
    ml <- .fcr_set_final_hessian(ml, fn, gr, "BFGS")
  }
  est <- ml$estimate
  names(est) <- names(start)
  s <- state(est)
  dimnames(ml$hessian) <- list(names(est), names(est))
  alpha <- stats::setNames(s$alpha, levels(g))
  .fcr_fe_with_smooth(
    list(estimate = est, loglik = sum(s$ll), hessian = ml$hessian,
         gradientObs = gr(est), fixef = list(fe1 = alpha),
         iterations = ml$iterations,
         converged = ml$code %in% c(0L, 1L, 2L, 8L),
         hessian_accuracy = ml$hessian_accuracy), sf)
}

# Copy the smooth summaries of a penalized fit `sf` into a FE fit list.
.fcr_fe_with_smooth <- function(out, sf) {
  if (is.null(sf)) return(out)
  keep <- c("smooth_table", "smooth_info", "smooth_lambda",
            "smooth_converged", "smooth_iterations", "smooth_coefs", "edf",
            "df_total")
  out[keep] <- sf[keep]
  out
}

# Predictions of a countreg_fe() fit: exp(x'beta + offset + effects).
.fcr_fe_predict <- function(object, data, newgroup = c("none", "mean_fe")) {
  newgroup <- match.arg(newgroup)
  m <- object$model
  fd <- forge_new_data(.fcr_blueprint(object), data)
  b <- m$estimate
  xb <- names(b)[names(b) %in% colnames(fd$X)]
  eta <- as.vector(fd$X[, xb, drop = FALSE] %*% b[xb]) + fd$offset
  tot <- rep(0, nrow(data))
  unknown <- rep(FALSE, nrow(data))
  for (d in seq_along(m$fixef_vars)) {
    v <- m$fixef_vars[d]
    if (!v %in% names(data)) {
      stop("predict() needs the fixed-effect variable `", v, "` in ",
           "`newdata`.", call. = FALSE)
    }
    a <- m$fixef[[d]]
    val <- a[match(as.character(data[[v]]), names(a))]
    new <- is.na(val) & !is.na(data[[v]])
    if (any(new) && newgroup == "mean_fe") {
      # A fixed mean of the estimated effects, weighted by the groups' sizes
      # in the estimation sample (plain mean for fits that predate the
      # stored sizes), so a new group's prediction does not depend on the
      # other rows of newdata
      n_g <- m$fixef_n[[d]]
      val[new] <- if (length(n_g) == length(a)) stats::weighted.mean(a, n_g)
        else mean(a)
    }
    unknown <- unknown | is.na(val)
    tot <- tot + ifelse(is.na(val), 0, val)
  }
  out <- exp(eta + tot)
  if (any(unknown & !is.na(eta)) && newgroup == "none") {
    message(sum(unknown & !is.na(eta)), " row(s) belong to groups not in ",
            "the estimation data, which have no fixed effect: their ",
            "predictions are NA (use newgroup = \"mean_fe\" to use the mean ",
            "effect).")
  }
  out[unknown] <- NA_real_
  as.vector(.fcr_family_mean(m$family, out,
                             if (!is.null(get_params(m$family)[[1]]))
                               rep(exp(b[[get_params(m$family)[[1]]]]),
                                   length(out)),
                             if (!is.null(get_params(m$family)[[2]]))
                               rep(exp(b[[get_params(m$family)[[2]]]]),
                                   length(out))))
}
