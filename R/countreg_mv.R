# Multivariate count models: countreg_mv().
#
# J count outcomes observed on the same units (crashes by severity or type
# at a site). Each outcome has its own family, mean formula and
# distribution-parameter formulas. Stage 1: the equations are independent
# given the covariates, so the log-likelihood is the sum of the equations'
# log-likelihoods, and the fit equals separate countreg() fits; the
# dependent forms (correlated random parameters, copulas) build on the same
# equation engine.
#
# Parameters, in order: for each equation j, its mean coefficients
# ("<outcome>:<term>") and distribution parameters ("<outcome>:ln(alpha)",
# or "<outcome>:ln(alpha):<term>" with a formula).

#' Multivariate count regression
#'
#' Fits several count outcomes observed on the same units, such as crashes
#' by severity or by type at each site, each with its own count family,
#' mean formula and distribution-parameter formulas.
#'
#' @param formula either a formula with the outcomes on the left as
#'   `cbind(y_1, ..., y_J)` and one right-hand side for every equation, or a
#'   list of two-sided formulas, one per outcome (named by the list names,
#'   or by their responses). An `offset()` term (an exposure) belongs in the
#'   formulas.
#' @param data a data frame.
#' @param family the count family of each outcome (see [count_families]):
#'   one code for every outcome, or one per outcome (named by outcome or in
#'   their order). Families made with [new_count_family()] are accepted.
#' @param dis_param_formula_1,dis_param_formula_2 formulas for the
#'   outcomes' distribution parameters, as in [countreg()]: `NULL`
#'   (constant), one one-sided formula for every outcome that has the
#'   parameter, or a list named by outcome (outcomes not in the list keep a
#'   constant parameter).
#' @param offset optional names of offset columns (log scale) used by every
#'   equation, instead of `offset()` in the formulas.
#' @param weights optional frequency weights: a column name or a numeric
#'   vector.
#' @param group_var optional column naming the units for cluster-robust
#'   standard errors (such as a site observed over several years).
#' @param stderr `"normal"` (from the Hessian, the default), `"robust"`
#'   (sandwich) or `"cluster"` (sandwich by `group_var`).
#' @param method,max.iters,start.vals,verbose optimizer settings, as in
#'   [countreg()]; the default method is `"BFGS"`.
#' @param ndraws,engine integration settings of families without a closed
#'   form, as in [countreg()].
#' @param .setup_only internal use: return the model setup without fitting.
#'
#' @details
#' **Independent equations.** The outcomes are independent given the
#' covariates, so the log-likelihood is the sum of the equations'
#' log-likelihoods and every equation's estimates equal those of a
#' separate [countreg()] fit. Fitting them together gives one object with
#' joint predictions, simulations and tests across equations (for
#' example, a Wald test that a coefficient is the same in two equations).
#' The rows used are those complete for every equation.
#'
#' Different outcomes often need different families: frequent property
#' damage crashes may be close to NB2, while rare fatal crashes are more
#' skewed and zero-heavy. Each equation can use any family of the package.
#'
#' @return An object of class `"countreg_mv"`, with methods `coef()`,
#'   `vcov()`, `logLik()`, `summary()`, `predict()`, `fitted()`,
#'   `residuals()`, `simulate()`, `tidy()`, `glance()` and `confint()`.
#'
#' @seealso [predict.countreg_mv()], [countreg()], [countreg_alloc()] for a
#'   total and its split across categories, [negmultinom()]
#' @examples
#' \donttest{
#' data("washington_roads")
#' w <- washington_roads
#' w$pdo <- w$Total_crashes - w$Injury_crashes
#' m <- countreg_mv(cbind(pdo, Injury_crashes) ~ lnaadt + speed50 +
#'                    offset(lnlength), data = w,
#'                  family = c(pdo = "NB2", Injury_crashes = "PIG"))
#' summary(m)
#' head(predict(m))
#' }
#' @export
countreg_mv <- function(formula, data, family = "NB2",
                        dis_param_formula_1 = NULL,
                        dis_param_formula_2 = NULL, offset = NULL,
                        weights = NULL, group_var = NULL, stderr = "normal",
                        method = "BFGS", max.iters = 1000, start.vals = NULL,
                        ndraws = 1500, engine = "quadrature", verbose = FALSE,
                        .setup_only = FALSE) {
  cl <- match.call()
  data <- as.data.frame(data)
  method <- .fcr_check_method(method)
  .fcr_check_controls(max.iters, ndraws)
  engine <- match.arg(engine, c("quadrature", "halton", "poilog"))
  stderr <- match.arg(stderr, c("normal", "robust", "cluster"))
  if (stderr == "cluster" && is.null(group_var)) {
    stop("`stderr = \"cluster\"` needs `group_var`.", call. = FALSE)
  }

  # ---- Equations ------------------------------------------------------------
  forms <- .fcr_mv_formulas(formula, data)
  outcomes <- names(forms)
  J <- length(forms)
  fams <- .fcr_mv_per_outcome(family, outcomes, "family")
  family_specs <- lapply(fams, function(f)
    if (inherits(f, "fcr_count_family")) f)
  fams <- vapply(fams, .fcr_resolve_family, character(1))
  d1 <- .fcr_mv_dis_formulas(dis_param_formula_1, outcomes,
                             "dis_param_formula_1")
  d2 <- .fcr_mv_dis_formulas(dis_param_formula_2, outcomes,
                             "dis_param_formula_2")
  for (j in seq_len(J)) {
    params <- get_params(fams[j])
    # One formula given for every outcome only applies where the family has
    # the parameter (a formula named for an outcome must fit its family)
    if (is.null(params[[1]]) && inherits(dis_param_formula_1, "formula")) {
      d1[j] <- list(NULL)
    }
    if (is.null(params[[2]]) && inherits(dis_param_formula_2, "formula")) {
      d2[j] <- list(NULL)
    }
    .fcr_check_dis_formulas(fams[j], params, d1[[j]], d2[[j]])
  }
  if (engine == "poilog" && "PLN" %in% fams) {
    rlang::check_installed("poilog", reason = "for `engine = \"poilog\"`.")
  }

  # ---- Data: the rows complete for every equation ----------------------------
  prep_eq <- function(j, dat, w) {
    prepare_model_data(forms[[j]], dat,
                       aux_formulas = list(dis1 = d1[[j]], dis2 = d2[[j]]),
                       offset = offset, weights = w, group_var = group_var)
  }
  rows <- seq_len(nrow(data))
  for (j in seq_len(J)) {
    p1 <- prep_eq(j, data, weights)
    rows <- intersect(rows, p1$rows_used)
  }
  if (!length(rows)) stop("No row is complete for every equation.",
                          call. = FALSE)
  n_dropped <- nrow(data) - length(rows)
  if (n_dropped > 0) {
    message(n_dropped, " row(s) with missing values were removed.")
  }
  dsub <- data[rows, , drop = FALSE]
  wsub <- if (is.numeric(weights)) weights[rows] else weights
  preps <- lapply(seq_len(J), function(j) prep_eq(j, dsub, wsub))
  if (any(vapply(preps, .fcr_prep_has_smooth, logical(1)))) {
    stop("Smooth terms are not yet available in countreg_mv().",
         call. = FALSE)
  }
  w <- preps[[1]]$weights

  # ---- Parameter layout --------------------------------------------------------
  haltons <- randtoolbox::halton(ndraws)
  dens_args <- list(haltons = haltons, normed_haltons = stats::qnorm(haltons),
                    engine = engine)
  eqs <- vector("list", J)
  names_all <- character(0)
  for (j in seq_len(J)) {
    pr <- preps[[j]]
    S <- .fcr_alloc_layout(fams[j], get_params(fams[j]), pr$X, pr$aux$dis1,
                           pr$aux$dis2, matrix(0, nrow(pr$X), 0), NULL, 2L,
                           c("a", "b"), FALSE)
    off <- length(names_all)
    nm <- paste0(outcomes[j], ":", S$names[seq_len(S$n_count)])
    eqs[[j]] <- list(
      outcome = outcomes[j], family = fams[j], y = as.numeric(pr$y),
      X = pr$X, off = pr$offset, W1 = pr$aux$dis1, W2 = pr$aux$dis2,
      S = list(idx_beta = off + S$idx_beta, idx_a = off + S$idx_a,
               idx_s = off + S$idx_s, n1 = S$n1, n2 = S$n2, kx = S$kx),
      dens_args = dens_args)
    names_all <- c(names_all, nm)
  }
  E <- list(eqs = eqs, J = J, n = length(rows), w = w, names = names_all,
            outcomes = outcomes)
  if (.setup_only) return(E)

  # ---- Fit ---------------------------------------------------------------------
  start <- if (is.null(start.vals)) .fcr_mv_start(E) else
    .fcr_check_start(start.vals, names_all)
  names(start) <- names_all
  fn <- function(p) .fcr_mv_ll(p, E, grad = FALSE)$ll
  gr <- function(p) .fcr_mv_ll(p, E, grad = TRUE)$G
  fit <- maxLik::maxLik(fn, grad = gr, start = start, method = method,
                        finalHessian = method == "BHHH",
                        control = list(iterlim = max.iters,
                                       printLevel = if (verbose) 2 else 0))
  # A few Newton steps (numerical Hessian of the analytic gradient) tighten
  # the solution after a quasi-Newton optimizer; kept only if the
  # log-likelihood does not decrease
  if (method != "NR") {
    f2 <- tryCatch(maxLik::maxLik(fn, grad = gr, start = fit$estimate,
                                  method = "NR",
                                  control = list(iterlim = 20)),
                   error = function(e) NULL)
    if (!is.null(f2) && is.finite(f2$maximum) &&
        f2$maximum >= fit$maximum - 1e-8 &&
        isTRUE(.fcr_converged(f2, "NR"))) {
      fit <- f2
      fit$refined_from <- method
    }
  }
  fit <- .fcr_set_final_hessian(fit, fn, gr, method)
  est <- fit$estimate
  names(est) <- names_all
  conv <- isTRUE(.fcr_converged(fit, if (is.null(fit$refined_from)) method
                                 else "NR"))
  if (!conv) {
    warning("The optimizer did not report convergence: ", fit$message,
            call. = FALSE)
  }
  U <- gr(est)
  colnames(U) <- names_all
  Hv <- .fcr_hessian_vcov(list(hessian = fit$hessian, estimate = est))
  V <- switch(stderr,
    normal = Hv,
    robust = {
      keep <- w > 0
      Hv %*% crossprod(U[keep, , drop = FALSE] / sqrt(w[keep])) %*% Hv
    },
    cluster = {
      Ug <- rowsum(U, preps[[1]]$group)
      G <- nrow(Ug)
      Hv %*% (crossprod(Ug) * G / max(G - 1, 1)) %*% Hv
    })
  dimnames(V) <- list(names_all, names_all)

  out <- list(
    call = cl, coefficients = est, vcov = V, loglik = sum(fn(est)),
    gradientObs = U, hessian = fit$hessian, converged = conv,
    iterations = unname(fit$iterations %||% NA), message = fit$message,
    outcomes = outcomes, family = stats::setNames(fams, outcomes),
    family_spec = family_specs, formulas = lapply(preps, `[[`, "formula"),
    dis_param_formula_1 = lapply(preps, function(p) p$aux_formulas$dis1),
    dis_param_formula_2 = lapply(preps, function(p) p$aux_formulas$dis2),
    layout = lapply(eqs, `[[`, "S"),
    blueprints = lapply(preps, `[[`, "blueprint"), stderr = stderr,
    group_var = group_var, y = vapply(eqs, `[[`, numeric(length(rows)), "y"),
    weights = w, nobs = sum(w > 0), n_obs_ic = sum(w),
    rows_used = rows, n_dropped = n_dropped,
    data_used = dsub, engine = engine, ndraws = ndraws)
  colnames(out$y) <- outcomes
  names(out$formulas) <- outcomes
  class(out) <- "countreg_mv"
  out$fitted.values <- predict(out, newdata = dsub)
  out
}

# The equations' formulas, named by outcome: from cbind(...) ~ rhs or a list
# of two-sided formulas.
.fcr_mv_formulas <- function(formula, data) {
  if (is.list(formula) && !inherits(formula, "formula")) {
    forms <- lapply(formula, stats::as.formula)
    if (length(forms) < 2L) {
      stop("`formula` needs two or more outcomes.", call. = FALSE)
    }
    if (any(vapply(forms, length, 1L) != 3L)) {
      stop("Every formula in the list must have an outcome on the left.",
           call. = FALSE)
    }
    nm <- names(forms)
    resp <- vapply(forms, function(f) paste(deparse(f[[2L]]), collapse = ""),
                   character(1))
    if (is.null(nm)) nm <- resp
    nm[!nzchar(nm)] <- resp[!nzchar(nm)]
  } else {
    formula <- stats::as.formula(formula)
    lhs <- if (length(formula) == 3L) formula[[2L]]
    if (is.null(lhs) || !is.call(lhs) || !identical(lhs[[1L]], as.name("cbind")) ||
        length(lhs) < 3L) {
      stop("`formula` needs two or more outcomes, as cbind(y_1, y_2) ~ x, ",
           "or a list of formulas, one per outcome.", call. = FALSE)
    }
    parts <- as.list(lhs)[-1L]
    nm <- names(parts)
    resp <- vapply(parts, function(e) paste(deparse(e), collapse = ""),
                   character(1))
    if (is.null(nm)) nm <- resp
    nm[!nzchar(nm)] <- resp[!nzchar(nm)]
    forms <- lapply(parts, function(e) {
      f <- formula
      f[[2L]] <- e
      f
    })
  }
  if (anyDuplicated(nm)) {
    stop("The outcomes must have different names: ", toString(nm[duplicated(nm)]),
         ".", call. = FALSE)
  }
  stats::setNames(forms, nm)
}

# One value per outcome from a single value or a vector/list named by
# outcome (or in their order).
.fcr_mv_per_outcome <- function(x, outcomes, what) {
  J <- length(outcomes)
  if (inherits(x, "fcr_count_family")) x <- list(x)
  x <- as.list(x)
  if (length(x) == 1L && is.null(names(x))) return(rep(x, J))
  if (!is.null(names(x)) && all(nzchar(names(x)))) {
    bad <- setdiff(names(x), outcomes)
    miss <- setdiff(outcomes, names(x))
    if (length(bad) || length(miss)) {
      stop("Names of `", what, "` must be the outcomes (",
           toString(outcomes), ").",
           if (length(miss)) paste0(" Missing: ", toString(miss), "."),
           if (length(bad)) paste0(" Unknown: ", toString(bad), "."),
           call. = FALSE)
    }
    return(unname(x[outcomes]))
  }
  if (length(x) != J) {
    stop("`", what, "` needs one value, or one per outcome (", J, ").",
         call. = FALSE)
  }
  unname(x)
}

# Distribution-parameter formulas per outcome: NULL, one formula for all, or
# a list named by outcome.
.fcr_mv_dis_formulas <- function(f, outcomes, what) {
  J <- length(outcomes)
  if (is.null(f)) return(vector("list", J))
  if (inherits(f, "formula")) return(rep(list(f), J))
  if (!is.list(f) || is.null(names(f)) || any(!nzchar(names(f)))) {
    stop("`", what, "` must be a formula or a list of formulas named by ",
         "outcome.", call. = FALSE)
  }
  bad <- setdiff(names(f), outcomes)
  if (length(bad)) {
    stop("Unknown outcome(s) in `", what, "`: ", toString(bad), ".",
         call. = FALSE)
  }
  out <- vector("list", J)
  for (j in seq_len(J)) {
    if (!is.null(f[[outcomes[j]]])) out[[j]] <- stats::as.formula(f[[outcomes[j]]])
  }
  out
}

# Row log-likelihoods of equation `eq` (unweighted) and, with grad, the
# derivatives with respect to its log mean and log distribution parameters.
# `add_eta` is added to the log mean (a random term).
.fcr_mv_eq_parts <- function(p, eq, add_eta = NULL, grad = TRUE) {
  S <- eq$S
  n <- length(eq$y)
  fam <- .fcr_get_family(eq$family)
  lp <- as.vector(eq$X %*% p[S$idx_beta]) + eq$off
  if (!is.null(add_eta)) lp <- lp + add_eta
  mu <- exp(lp)
  dp <- .fcr_alloc_dispar(p, S, eq$W1, eq$W2, n)
  ll <- tryCatch(.fcr_quiet_dens(as.vector(do.call(fam$dens, c(
    list(y = eq$y, predicted = mu, alpha = dp$alpha, sigma = dp$sigma),
    eq$dens_args, list(log = TRUE)))), log = TRUE),
    error = function(e) rep(NA_real_, n))
  ll[!.fcr_family_valid(eq$family, mu, dp$alpha, dp$sigma)] <- -Inf
  if (!grad) return(list(ll = ll, mu = mu, alpha = dp$alpha, sigma = dp$sigma))
  sc <- .fcr_family_score(eq$family, eq$y, mu, dp$alpha, dp$sigma,
                          dens_args = eq$dens_args)
  list(ll = ll, mu = mu, alpha = dp$alpha, sigma = dp$sigma, s_eta = sc$eta,
       s_a = sc$alpha, s_s = sc$sigma)
}

# Gradient rows of equation `eq`'s own parameters, from its parts and the
# row weights w.
.fcr_mv_eq_grad <- function(pt, eq, w, G) {
  S <- eq$S
  G[, S$idx_beta] <- (w * pt$s_eta) * eq$X
  if (S$n1) G[, S$idx_a] <- (w * pt$s_a) * (if (is.null(eq$W1)) 1 else eq$W1)
  if (S$n2) G[, S$idx_s] <- (w * pt$s_s) * (if (is.null(eq$W2)) 1 else eq$W2)
  G
}

# Log-likelihood of independent equations: one contribution per row (the
# sum over equations, times the row weight), and the gradient rows.
.fcr_mv_ll <- function(p, E, grad = FALSE) {
  w <- E$w
  ll <- numeric(E$n)
  G <- if (grad) matrix(0, E$n, length(p))
  for (eq in E$eqs) {
    pt <- .fcr_mv_eq_parts(p, eq, grad = grad)
    ll <- ll + pt$ll
    if (grad) G <- .fcr_mv_eq_grad(pt, eq, w, G)
  }
  ll <- w * ll
  ll[w == 0] <- 0
  if (!grad) {
    ll[!is.finite(ll)] <- NA_real_
    return(list(ll = ll))
  }
  G[!is.finite(G)] <- 0
  list(ll = ll, G = G)
}

# Starting values: a Poisson fit of each equation and the family's
# dispersion start.
.fcr_mv_start <- function(E) {
  st <- numeric(length(E$names))
  for (eq in E$eqs) {
    S <- eq$S
    pm <- suppressWarnings(stats::glm.fit(eq$X, eq$y, weights = E$w,
                                          offset = eq$off,
                                          family = stats::poisson()))
    b <- pm$coefficients
    b[is.na(b)] <- 0
    st[S$idx_beta] <- b
    fs <- .fcr_get_family(eq$family)$start
    if (S$n1) {
      a0 <- .fcr_moment_disp_start(eq$family, eq$y, pm$fitted.values, E$w) %||%
        fs[1L]
      pos <- if (is.null(eq$W1)) 1L else match("(Intercept)", colnames(eq$W1))
      if (!is.na(pos)) st[S$idx_a[pos]] <- a0
    }
    if (S$n2) {
      pos <- if (is.null(eq$W2)) 1L else match("(Intercept)", colnames(eq$W2))
      if (!is.na(pos)) st[S$idx_s[pos]] <- rep_len(fs, 2L)[2L]
    }
  }
  st
}

# The equations of a fit for new data.
.fcr_mv_design <- function(object, newdata) {
  haltons <- randtoolbox::halton(object$ndraws)
  dens_args <- list(haltons = haltons, normed_haltons = stats::qnorm(haltons),
                    engine = object$engine)
  eqs <- lapply(seq_along(object$outcomes), function(j) {
    fd <- forge_new_data(object$blueprints[[j]], newdata)
    list(outcome = object$outcomes[j], family = object$family[[j]],
         y = rep(0, nrow(fd$X)), X = fd$X, off = fd$offset,
         W1 = fd$aux$dis1, W2 = fd$aux$dis2, S = object$layout[[j]],
         dens_args = dens_args)
  })
  list(eqs = eqs, J = length(eqs), n = nrow(newdata),
       outcomes = object$outcomes)
}

# Mean, distribution parameters and variance of every equation for the rows
# of `newdata`.
.fcr_mv_values <- function(object, newdata) {
  E <- .fcr_mv_design(object, newdata)
  p <- object$coefficients
  lapply(E$eqs, function(eq) {
    pt <- .fcr_mv_eq_parts(p, eq, grad = FALSE)
    list(eq = eq, mu = pt$mu, alpha = pt$alpha, sigma = pt$sigma,
         mean = .fcr_family_mean(eq$family, pt$mu, pt$alpha, pt$sigma))
  })
}

# ---------------------------------------------------------------------------
# Methods
# ---------------------------------------------------------------------------

#' Predictions from multivariate count models
#'
#' @param object a [countreg_mv()] fit.
#' @param newdata optional new data; the estimation data by default.
#' @param type `"response"` (default): a matrix of the expected counts, one
#'   column per outcome; `"link"`: the linear predictors (log of the mean
#'   parameter); `"variance"`: the variances; `"prob"`: a list with one
#'   matrix per outcome of the probabilities of the counts `0:max_count`.
#' @param max_count for `type = "prob"`: the largest count (default: the
#'   largest observed count of each outcome).
#' @param ... unused.
#' @return A matrix with one column per outcome, or a list of matrices
#'   (`"prob"`).
#' @export
predict.countreg_mv <- function(object, newdata = NULL,
                                type = c("response", "link", "variance",
                                         "prob"),
                                max_count = NULL, ...) {
  type <- match.arg(type)
  if (is.null(newdata)) newdata <- object$data_used
  vals <- .fcr_mv_values(object, as.data.frame(newdata))
  outs <- object$outcomes
  if (type == "prob") {
    out <- lapply(seq_along(vals), function(j) {
      v <- vals[[j]]
      mc <- max_count %||% max(object$y[, j])
      fam <- .fcr_get_family(v$eq$family)
      P <- vapply(0:mc, function(k) exp(.fcr_quiet_dens(as.vector(do.call(
        fam$dens, c(list(y = rep(k, length(v$mu)), predicted = v$mu,
                         alpha = v$alpha, sigma = v$sigma), v$eq$dens_args,
                    list(log = TRUE)))), log = TRUE)), numeric(length(v$mu)))
      P <- matrix(P, nrow = length(v$mu))
      colnames(P) <- 0:mc
      P
    })
    return(stats::setNames(out, outs))
  }
  M <- vapply(vals, function(v) switch(type,
    response = v$mean,
    link = log(v$mu),
    variance = .fcr_family_variance(v$eq$family, v$mu, v$alpha, v$sigma)),
    numeric(nrow(as.data.frame(newdata))))
  M <- matrix(M, ncol = length(outs))
  colnames(M) <- outs
  M
}

#' @export
fitted.countreg_mv <- function(object, ...) object$fitted.values

#' @export
residuals.countreg_mv <- function(object, type = c("response", "pearson"),
                                  ...) {
  type <- match.arg(type)
  r <- object$y - object$fitted.values
  if (type == "pearson") {
    r <- r / sqrt(predict(object, newdata = object$data_used,
                          type = "variance"))
  }
  r
}

#' @export
coef.countreg_mv <- function(object, ...) object$coefficients

#' @export
vcov.countreg_mv <- function(object, ...) object$vcov

#' @export
nobs.countreg_mv <- function(object, ...) object$nobs

#' @export
logLik.countreg_mv <- function(object, ...) {
  structure(object$loglik, df = length(object$coefficients),
            nobs = object$n_obs_ic, class = "logLik")
}

#' @export
print.countreg_mv <- function(x, ...) {
  cat("Multivariate count model:", length(x$outcomes), "outcomes\n")
  cat("Families:", paste0(x$outcomes, " = ",
                          vapply(x$family, .fcr_family_label, ""),
                          collapse = ", "), "\n\n")
  print(round(x$coefficients, 4))
  invisible(x)
}

# The outcome (equation) of each coefficient.
.fcr_mv_coef_outcome <- function(object) {
  out <- character(length(object$coefficients))
  for (j in seq_along(object$outcomes)) {
    S <- object$layout[[j]]
    out[c(S$idx_beta, S$idx_a, S$idx_s)] <- object$outcomes[j]
  }
  out
}

#' @export
summary.countreg_mv <- function(object, ...) {
  est <- object$coefficients
  se <- sqrt(pmax(diag(object$vcov), 0))
  z <- est / se
  tab <- data.frame(Estimate = est, `Std. Error` = se, `z value` = z,
                    `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)),
                    check.names = FALSE)
  structure(list(table = tab, part = .fcr_mv_coef_outcome(object),
                 object = object), class = "summary.countreg_mv")
}

#' @export
print.summary.countreg_mv <- function(x, digits = 4, ...) {
  o <- x$object
  cat("Multivariate count model\n\n")
  cat("Rows:       ", o$nobs, "\n")
  cat("Std. errors:", o$stderr, "\n")
  cat("Equations:   independent given the covariates\n\n")
  for (j in seq_along(o$outcomes)) {
    oc <- o$outcomes[j]
    cat(oc, " (", .fcr_family_label(o$family[[j]]), "):\n", sep = "")
    tb <- x$table[x$part == oc, , drop = FALSE]
    rownames(tb) <- sub(paste0("^", .fcr_regex_escape(oc), ":"), "",
                        rownames(tb))
    stats::printCoefmat(tb, digits = digits, signif.legend = FALSE)
    cat("\n")
  }
  k <- length(o$coefficients)
  cat("Log-likelihood:", format(o$loglik, digits = 8), " AIC:",
      format(-2 * o$loglik + 2 * k, digits = 8), " BIC:",
      format(-2 * o$loglik + log(o$n_obs_ic) * k, digits = 8), "\n")
  if (!o$converged) cat("Warning: the optimizer did not converge.\n")
  invisible(x)
}

.fcr_regex_escape <- function(x) gsub("([][{}()+*^$|\\\\?.])", "\\\\\\1", x)

#' @export
tidy.countreg_mv <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  s <- summary(x)
  oc <- s$part
  term <- mapply(function(t, o) sub(paste0("^", .fcr_regex_escape(o), ":"),
                                    "", t), rownames(s$table), oc,
                 USE.NAMES = FALSE)
  comp <- rep("mean", length(term))
  for (j in seq_along(x$outcomes)) {
    S <- x$layout[[j]]
    comp[c(S$idx_a, S$idx_s)] <- "distribution"
  }
  out <- data.frame(outcome = oc, term = term, component = comp,
                    estimate = s$table[, 1], std.error = s$table[, 2],
                    statistic = s$table[, 3], p.value = s$table[, 4],
                    stringsAsFactors = FALSE, row.names = NULL)
  if (conf.int) {
    zq <- stats::qnorm(1 - (1 - conf.level) / 2)
    out$conf.low <- out$estimate - zq * out$std.error
    out$conf.high <- out$estimate + zq * out$std.error
  }
  tibble::as_tibble(out)
}

#' @export
glance.countreg_mv <- function(x, ...) {
  k <- length(x$coefficients)
  tibble::tibble(logLik = x$loglik, AIC = -2 * x$loglik + 2 * k,
                 BIC = -2 * x$loglik + log(x$n_obs_ic) * k, nobs = x$nobs,
                 npar = k, outcomes = length(x$outcomes),
                 converged = x$converged)
}

#' @export
simulate.countreg_mv <- function(object, nsim = 1, seed = NULL,
                                 newdata = NULL, ...) {
  if (!is.null(seed)) set.seed(seed)
  if (is.null(newdata)) newdata <- object$data_used
  vals <- .fcr_mv_values(object, as.data.frame(newdata))
  n <- nrow(as.data.frame(newdata))
  one <- function() {
    Y <- vapply(vals, function(v) .fcr_family_r(v$eq$family, n, v$mu,
                                                v$alpha, v$sigma), numeric(n))
    Y <- matrix(Y, nrow = n)
    colnames(Y) <- object$outcomes
    Y
  }
  stats::setNames(replicate(nsim, one(), simplify = FALSE),
                  paste0("sim_", seq_len(nsim)))
}
