# Confidence, prediction and tolerance intervals for countreg() models.
#
# Shared machinery:
#   .fcr_parts_pmf()       P(Y = k), k = 0..K, of the observed (reported)
#                          count for each row of a parts object, including
#                          zero-inflated/hurdle components;
#   .fcr_pmf_resolved()    the same with K grown until the omitted upper tail
#                          is negligible (or a count limit is reached);
#   .fcr_param_draws()     parameter vectors representing estimation
#                          uncertainty: quasi-random (Halton) draws from the
#                          asymptotic normal distribution of the estimates,
#                          or parametric-bootstrap refits.

# ---- probabilities of the observed count -------------------------------------

.fcr_parts_pmf <- function(object, parts, K) {
  m <- object$model
  n <- length(parts$mu_link)
  P <- matrix(NA_real_, n, K + 1)
  mean_par <- parts$mu_link * parts$p_report
  ok <- .fcr_parts_valid(parts)
  if (any(ok)) {
    no <- sum(ok)
    rk <- function(v) if (!is.null(v)) rep(v[ok], K + 1)
    haltons <- randtoolbox::halton(m$ndraws %||% 1500)
    p <- get_probFunc(parts$family)(
      y = rep(0:K, each = no), predicted = rep(mean_par[ok], K + 1),
      alpha = rk(parts$alpha), sigma = rk(parts$sigma), haltons = haltons,
      normed_haltons = stats::qnorm(haltons),
      engine = m$engine %||% "halton", log = FALSE)
    P[ok, ] <- matrix(p, no, K + 1)
  }
  if (!is.null(parts$zero_prob)) {
    zp <- parts$zero_prob
    f0 <- P[, 1]
    if (parts$zero_type == "inflated") {
      P <- P * (1 - zp)
      P[, 1] <- zp + (1 - zp) * f0
    } else {
      lf0 <- parts$lf0 %||% log(f0)
      P <- P * (zp / -expm1(lf0))
      P[, 1] <- 1 - zp
    }
  }
  P
}

# Rows of a parts object at which the count distribution is defined: no
# missing mean, distribution, reporting or zero-model parameter, and
# parameters inside the family's domain (.fcr_family_valid()). Other rows
# have no distribution, so their probabilities are NA rather than the zeros
# the likelihood's quiet densities would give.
.fcr_parts_valid <- function(parts) {
  mean_par <- parts$mu_link * parts$p_report
  ok <- !is.na(mean_par)
  for (nm in c("alpha", "sigma", "zero_prob")) {
    if (!is.null(parts[[nm]])) ok <- ok & !is.na(parts[[nm]])
  }
  ok & .fcr_family_valid(parts$family, mean_par, parts$alpha, parts$sigma)
}

# PMF on 0..K with K large enough that the omitted mass P(Y > K) is below
# `tol` for every row, starting from mean + 20 SD and doubling up to
# `max_count`. Returns list(P, K, resolved).
.fcr_pmf_resolved <- function(object, parts, max_count = 1e5, tol = 1e-10) {
  mean_par <- parts$mu_link * parts$p_report
  v <- tryCatch(.fcr_family_variance(parts$family, mean_par, parts$alpha,
                                     parts$sigma),
                error = function(e) NULL)
  m <- parts$response
  K <- if (is.null(v)) NA else max(m + 20 * sqrt(pmax(v, 0)), na.rm = TRUE)
  if (!is.finite(K)) K <- 200
  K <- min(max(ceiling(K) + 30, 30), max_count)
  repeat {
    P <- .fcr_parts_pmf(object, parts, K)
    miss <- 1 - rowSums(P)
    unresolved <- !is.na(miss) & miss > tol
    if (!any(unresolved) || K >= max_count) break
    K <- min(2 * K, max_count)
  }
  list(P = P, K = K, resolved = !unresolved)
}

# ---- parameter uncertainty -----------------------------------------------------

.fcr_is_pd <- function(V) {
  all(is.finite(V)) &&
    min(eigen(V, symmetric = TRUE, only.values = TRUE)$values) > 0
}

# Parameter vectors representing estimation uncertainty (one per row).
#   "normal":    theta_hat + R'z, R'R = vcov, z = qnorm of ndraws Halton
#                points (deterministic; scrambled with a fixed seed when
#                requested). Requires a positive definite covariance matrix;
#                otherwise the parametric bootstrap is used, with a message.
#   "bootstrap": parametric bootstrap: B response vectors simulated from the
#                fitted model on the training rows, each refit with the
#                model's frozen specification; refits that fail, do not
#                converge or return non-finite or differently named
#                coefficients are rejected and counted.
.fcr_param_draws <- function(object, uncertainty, ndraws, bootstraps, seed,
                             scrambled = FALSE) {
  theta <- object$model$estimate
  V <- stats::vcov(object)
  if (uncertainty == "normal" && !.fcr_is_pd(V)) {
    message("The covariance matrix of the estimates is not positive ",
            "definite (e.g. a parameter at a boundary), so the normal ",
            "approximation is unavailable; using the parametric bootstrap ",
            "(", bootstraps, " refits).")
    uncertainty <- "bootstrap"
  }
  if (uncertainty == "normal") {
    # One point set (no per-unit blocks), scrambled to break the
    # correlation between high dimensions of the standard sequence
    z <- stats::qnorm(.fcr_halton(ndraws, length(theta), scrambled,
                                  seed = seed, scramble = "mixed"))
    draws <- sweep(z %*% chol(V), 2, theta, `+`)
    colnames(draws) <- names(theta)
    return(list(draws = draws, uncertainty = "normal", requested = ndraws,
                accepted = ndraws, failures = character(0)))
  }
  .fcr_param_bootstrap(object, bootstraps, seed)
}

.fcr_param_bootstrap <- function(object, B, seed) {
  w <- object$model[["weights"]]
  if (!is.null(w) && any(w != 1)) {
    stop("The parametric bootstrap is not available for models with ",
         "weights.", call. = FALSE)
  }
  resp <- object$formula[[2]]
  if (!is.name(resp)) {
    stop("The parametric bootstrap needs the response to be a variable ",
         "name in the formula.", call. = FALSE)
  }
  d <- .fcr_training_rows(object)
  ysim <- stats::simulate(object, nsim = B, seed = seed)
  theta <- object$model$estimate
  draws <- matrix(NA_real_, B, length(theta),
                  dimnames = list(NULL, names(theta)))
  why <- character(B)
  for (b in seq_len(B)) {
    db <- d
    db[[as.character(resp)]] <- ysim[[b]]
    fit <- tryCatch(suppressWarnings(suppressMessages(
      .margEff_refit(object, db, start = theta))), error = function(e) e)
    why[b] <- if (inherits(fit, "error")) "error" else
      if (!isTRUE(fit$model$converged)) "not converged" else
        if (!identical(names(fit$model$estimate), names(theta))) "coefficients" else
          if (any(!is.finite(fit$model$estimate))) "non-finite" else ""
    if (why[b] == "") draws[b, ] <- fit$model$estimate
  }
  ok <- why == ""
  if (mean(ok) < 0.9) {
    warning(sum(!ok), " of ", B, " bootstrap refits were rejected (",
            toString(names(table(why[!ok]))), "); the interval may be ",
            "unreliable.", call. = FALSE)
  }
  if (sum(ok) < 20) {
    stop("Too few bootstrap refits succeeded (", sum(ok), " of ", B, ").",
         call. = FALSE)
  }
  list(draws = draws[ok, , drop = FALSE], uncertainty = "bootstrap",
       requested = B, accepted = sum(ok), failures = why[!ok], seed = seed)
}

# ---- confidence and prediction intervals for predict() -----------------------

.fcr_predict_interval <- function(object, data, type, interval, level,
                                  uncertainty = "normal", ndraws = 500,
                                  bootstraps = 200, seed = NULL,
                                  scrambled = FALSE, max_count = 1e5) {
  if (!inherits(object, "countreg")) {
    stop("Intervals are available for countreg() models.", call. = FALSE)
  }
  .fcr_check_prob(level, "level")
  if (interval == "prediction" && type != "response") {
    stop("Prediction intervals are for the observed count ",
         "(type = \"response\").", call. = FALSE)
  }
  if (!type %in% c("response", "true")) {
    stop("Intervals are available for type = \"response\" or \"true\".",
         call. = FALSE)
  }
  fd <- .fcr_countreg_fd(object, data)
  est <- .fcr_parts_from_lp(object, .fcr_countreg_lp(object, fd))
  fit <- if (type == "response") est$response else est$mu_true
  pd <- .fcr_param_draws(object, uncertainty, ndraws, bootstraps, seed,
                         scrambled)
  M <- nrow(pd$draws)
  n <- length(fit)
  a <- (1 - level) / 2
  # Each draw is used only for the rows where its parameters define a
  # distribution (for example GP2's coupled domain): other draws are
  # rejected for that row and counted, never treated as zero mass
  parts_j <- function(j) suppressWarnings(.fcr_parts_from_lp(
    object, .fcr_countreg_lp(object, fd, pd$draws[j, ])))
  accepted <- numeric(n)
  unresolved <- rep(FALSE, n)
  if (interval == "confidence") {
    vals <- matrix(NA_real_, n, M)
    for (j in seq_len(M)) {
      pj <- parts_j(j)
      okj <- .fcr_parts_valid(pj)
      v <- if (type == "response") pj$response else pj$mu_true
      v[!okj] <- NA
      vals[, j] <- v
      accepted <- accepted + (okj & !is.na(v))
    }
    lwr <- apply(vals, 1, stats::quantile, probs = a, na.rm = TRUE,
                 names = FALSE)
    upr <- apply(vals, 1, stats::quantile, probs = 1 - a, na.rm = TRUE,
                 names = FALSE)
  } else {
    # Predictive distribution: the observed-count PMF averaged over the
    # draws accepted for each row
    S <- matrix(0, n, 1)
    for (j in seq_len(M)) {
      pj <- parts_j(j)
      res <- .fcr_pmf_resolved(object, pj, max_count)
      Pj <- res$P
      okj <- !is.na(Pj[, 1])
      Pj[!okj, ] <- 0
      accepted <- accepted + okj
      unresolved <- unresolved | (okj & !res$resolved)
      if (ncol(Pj) > ncol(S)) S <- cbind(S, matrix(0, n, ncol(Pj) - ncol(S)))
      S[, seq_len(ncol(Pj))] <- S[, seq_len(ncol(Pj))] + Pj
    }
    Fbar <- t(apply(S / accepted, 1, cumsum))
    if (n == 1) Fbar <- matrix(Fbar, 1)
    q <- function(prob) apply(Fbar, 1, function(f) {
      k <- which(f >= prob - 1e-12)
      if (length(k)) k[1] - 1 else NA_real_
    })
    lwr <- q(a)
    upr <- q(1 - a)
    # An upper limit beyond the largest count evaluated is unknown, not
    # absent: say so
    cut <- !is.na(fit) & accepted > 0 & is.na(upr)
    if (any(cut)) {
      warning("The upper prediction limit of ", sum(cut), " row(s) lies ",
              "above max_count = ", max_count, " (the predictive ",
              "distribution was not resolved there); it is NA. Increase ",
              "`max_count`.", call. = FALSE)
    }
  }
  # Rows without a distribution at the estimates (missing parameters) are NA
  # for that reason, not because of rejected draws
  low_acc <- !is.na(fit) & .fcr_parts_valid(est) & accepted < 0.9 * M
  if (any(low_acc)) {
    warning("For ", sum(low_acc), " row(s), fewer than 90% of the ", M,
            " parameter draws give a valid distribution (the others fall ",
            "outside the family's parameter space and were rejected); the ",
            "interval may be unreliable.", call. = FALSE)
  }
  lwr[accepted == 0] <- NA
  upr[accepted == 0] <- NA
  na <- is.na(fit)
  lwr[na] <- NA
  upr[na] <- NA
  out <- data.frame(fit = fit, lwr = lwr, upr = upr)
  attr(out, "interval") <- list(
    interval = interval, level = level, type = type,
    uncertainty = pd$uncertainty, draws = M, requested = pd$requested,
    failures = pd$failures, accepted_by_row = accepted,
    unresolved = if (interval == "prediction") unresolved)
  out
}

.fcr_check_prob <- function(x, nm) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || x <= 0 || x >= 1) {
    stop("`", nm, "` must be a single number between 0 and 1.",
         call. = FALSE)
  }
}

# ---- tolerance intervals -----------------------------------------------------

#' Tolerance Intervals for Count Regression Models
#'
#' Content tolerance intervals for the observed count of a [countreg()]
#' model: integer intervals intended to contain at least a fraction
#' `content` of the count distribution at given covariate values, with
#' confidence `confidence` over repeated samples.
#'
#' @param object a model fitted with [countreg()].
#' @param newdata a data frame of covariate profiles (scenarios). Defaults
#'   to the training rows.
#' @param content the required fraction \eqn{p} of the count distribution
#'   (e.g. 0.95).
#' @param confidence the confidence level \eqn{\gamma} (e.g. 0.95).
#' @param side `"two-sided"` (default), `"upper"` (interval `[0, U]`) or
#'   `"lower"` (interval `[L, Inf)`).
#' @param scope `"pointwise"` (default): each row's interval has the stated
#'   confidence on its own; `"simultaneous"`: the confidence holds jointly
#'   for all covariate profiles (based on a confidence set for all
#'   parameters, so wider).
#' @param method `"auto"` (default): `"exact"` when it applies, otherwise
#'   `"joint"`. `"joint"`: the narrowest integer interval whose content is
#'   at least `content` for every parameter value in the confidence set.
#'   `"envelope"`: equal-tail quantiles taken at the worst-case parameter
#'   values of each tail (simpler, more conservative). `"exact"`: the
#'   common-rate Poisson construction, available for intercept-only Poisson
#'   models without auxiliary components.
#' @param uncertainty how the confidence set is formed: `"normal"`
#'   (default), the asymptotic normal (Wald) confidence ellipsoid of the
#'   estimates, or `"bootstrap"`, a parametric-bootstrap confidence region.
#'   `"normal"` switches to the bootstrap, with a message, when the
#'   covariance matrix is not positive definite.
#' @param npoints number of points of the confidence set at which the
#'   distribution is evaluated (normal method). They are deterministic
#'   (based on Halton sequences), so results are reproducible.
#' @param bootstraps number of parametric-bootstrap refits.
#' @param seed random seed for the parametric bootstrap (the caller's
#'   random-number stream is not changed). `NULL` uses the current stream.
#' @param max_count largest count used when resolving the upper tail.
#'
#' @details
#' \strong{What the interval means.} For covariates \eqn{x}, an integer
#' interval \eqn{[L, U]} computed from the estimation sample \eqn{D} is a
#' content tolerance interval if
#' \deqn{\Pr_D\{F_\theta(U \mid x) - F_\theta(L - 1 \mid x) \ge p\} \ge
#'   \gamma,}{P_D{ F(U | x) - F(L - 1 | x) >= p } >= gamma,}
#' where \eqn{F_\theta} is the distribution function of the observed count
#' under the true parameters, \eqn{p} the content and \eqn{\gamma} the
#' confidence. This differs from a prediction interval (which covers one
#' future count on average; see [predict.flexCountReg()]) and from a
#' confidence interval for the mean. The statement is conditional on the
#' model being correct.
#'
#' \strong{Construction.} The distribution of the count at \eqn{x} depends
#' on the parameters only through the linear predictors of the model
#' components at \eqn{x} (mean, distribution-parameter formulas, reporting
#' and zero models). A confidence set for those linear predictors is formed
#' (Wald: an ellipsoid with \eqn{\chi^2} cutoff, with as many degrees of
#' freedom as linear predictors for pointwise intervals, or as parameters
#' for simultaneous intervals; bootstrap: the bootstrap draws closest to
#' the estimate in Mahalanobis distance, a fraction \eqn{\gamma} of them).
#' The distribution is evaluated at `npoints` points of the set (its
#' boundary, interior and center; deterministic Halton-based points) and
#' the interval is chosen so that its content is at least \eqn{p} at every
#' evaluated point. When the true parameter lies in the confidence set the
#' interval has the required content, so the confidence is approximately at
#' least \eqn{\gamma}; it is approximate because the confidence set is
#' asymptotic (or bootstrap-based) and the extrema are found on a finite
#' set of points. Discrete counts generally give more than the requested
#' content; endpoints are never interpolated.
#'
#' \strong{Exact Poisson method.} For an intercept-only Poisson model (a
#' common rate \eqn{\lambda}, with exposures from an offset), an exact
#' confidence interval for \eqn{\lambda} from [stats::poisson.test()] is
#' mapped to Poisson quantiles at the scenario's exposure, giving an
#' interval with at least the stated content and confidence (conservative
#' because of discreteness).
#'
#' In a Monte Carlo study of Poisson and NB2 regressions (n = 50 and 200,
#' 1000 samples each; content 0.90, confidence 0.95), the joint method
#' achieved confidence 0.99 to 1.00 and plug-in quantile intervals (which
#' ignore estimation uncertainty) as little as 0.54. The interval targets
#' the observed (reported) count, including zero-inflated, hurdle and
#' underreporting components. Random-parameter and panel models are not
#' yet supported.
#'
#' @returns A data frame with one row per row of `newdata`: `fit` (the
#'   expected count), `lower` and `upper` (integer limits; `upper` is `Inf`
#'   for `side = "lower"`), `min_content` (the smallest content of the
#'   interval over the evaluated points of the confidence set) and `status`
#'   (`"ok"`, or the reason a limit is unavailable). The settings are
#'   stored in the attribute `"tolerance"`.
#'
#' @references
#' Krishnamoorthy, K. and Mathew, T. (2009). \emph{Statistical Tolerance
#' Regions: Theory, Applications, and Computation}. Wiley.
#'
#' Young, D. S. (2010). tolerance: An R package for estimating tolerance
#' intervals. \emph{Journal of Statistical Software}, 36(5), 1-39.
#'
#' @seealso [predict.flexCountReg()] for confidence and prediction
#'   intervals.
#'
#' @examples
#' data("washington_roads")
#' fit <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
#'                 family = "NB2", method = "BFGS")
#' scen <- data.frame(lnaadt = log(c(2000, 10000)), lnlength = log(1))
#' tolerance_interval(fit, scen, content = 0.90, confidence = 0.95)
#'
#' # Exact common-rate Poisson method
#' set.seed(1)
#' d <- data.frame(y = rpois(30, 2))
#' p0 <- countreg(y ~ 1, data = d, family = "POISSON", method = "BFGS")
#' tolerance_interval(p0, data.frame(x = 1), content = 0.95)
#' @export
tolerance_interval <- function(object, newdata = NULL, content = 0.95,
                               confidence = 0.95,
                               side = c("two-sided", "upper", "lower"),
                               scope = c("pointwise", "simultaneous"),
                               method = c("auto", "joint", "envelope",
                                          "exact"),
                               uncertainty = c("normal", "bootstrap"),
                               npoints = 400, bootstraps = 200, seed = NULL,
                               max_count = 1e5) {
  if (!inherits(object, "countreg")) {
    stop("tolerance_interval() is available for countreg() models.",
         call. = FALSE)
  }
  .fcr_check_prob(content, "content")
  .fcr_check_prob(confidence, "confidence")
  side <- match.arg(side)
  scope <- match.arg(scope)
  method <- match.arg(method)
  uncertainty <- match.arg(uncertainty)
  .fcr_check_controls(max.iters = 1, ndraws = npoints)
  data <- if (is.null(newdata)) .fcr_training_rows(object) else
    as.data.frame(newdata)
  fd <- .fcr_countreg_fd(object, data)
  est <- .fcr_parts_from_lp(object, .fcr_countreg_lp(object, fd))
  exact_ok <- .fcr_ti_exact_applies(object, fd)
  if (method == "exact" && !exact_ok) {
    stop("method = \"exact\" needs an intercept-only Poisson model without ",
         "distribution-parameter, zero or underreporting components.",
         call. = FALSE)
  }
  if (method == "auto") method <- if (exact_ok) "exact" else "joint"
  settings <- list(content = content, confidence = confidence, side = side,
                   scope = scope, method = method)
  if (method == "exact") {
    out <- .fcr_ti_exact(object, fd, est, content, confidence, side)
    settings$uncertainty <- "exact"
    attr(out, "tolerance") <- settings
    return(.fcr_ti_class(out))
  }
  pd <- if (uncertainty == "bootstrap" ||
            !.fcr_is_pd(stats::vcov(object))) {
    .fcr_param_draws(object, "bootstrap", 0, bootstraps, seed)
  }
  settings$uncertainty <- if (is.null(pd)) "normal" else "bootstrap"
  if (!is.null(pd) && uncertainty == "normal") {
    message("The covariance matrix of the estimates is not positive ",
            "definite; using the parametric bootstrap.")
  }
  n <- length(est$mu_link)
  out <- data.frame(fit = est$response, lower = NA_real_, upper = NA_real_,
                    min_content = NA_real_, status = "ok",
                    stringsAsFactors = FALSE)
  comps <- .fcr_ti_components(object, fd)
  theta <- object$model$estimate
  V <- stats::vcov(object)
  k <- length(theta)
  for (i in seq_len(n)) {
    if (is.na(est$mu_link[i])) {
      out$status[i] <- "missing values"
      next
    }
    A <- .fcr_ti_A(comps, i, names(theta))
    off <- comps$offset[i]
    if (is.null(pd)) {
      eta_hat <- as.vector(A %*% theta)
      Ve <- A %*% V %*% t(A)
      ev <- eigen(Ve, symmetric = TRUE)
      # The ellipsoid lives in the span of A V A': its dimension is the rank
      # (a predictor direction with no uncertainty, such as x = 0 in a model
      # without an intercept, is known exactly and adds no degree of freedom)
      df <- if (scope == "pointwise") .fcr_psd_rank(ev$values) else k
      cutoff <- stats::qchisq(confidence, df)
      Lf <- ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), nrow(A))
      u <- .fcr_ball_points(nrow(A), npoints)
      eta <- sweep(sqrt(cutoff) * u %*% t(Lf), 2, eta_hat, `+`)
    } else {
      eta_b <- pd$draws %*% t(A)
      eta_hat <- as.vector(A %*% theta)
      dist <- if (scope == "pointwise") {
        .fcr_mahalanobis_psd(eta_b, eta_hat, stats::cov(eta_b))
      } else {
        .fcr_mahalanobis_psd(pd$draws, theta, stats::cov(pd$draws))
      }
      keep <- dist <= stats::quantile(dist, confidence, names = FALSE)
      eta <- rbind(eta_b[keep, , drop = FALSE], eta_hat)
    }
    lp <- .fcr_ti_lp(comps, eta, off)
    parts <- suppressWarnings(.fcr_parts_from_lp(object, lp))
    res <- .fcr_pmf_resolved(object, parts, max_count)
    # Points of the confidence set outside the family's parameter space
    # (NA rows) define no distribution and are dropped; the estimate itself
    # (the last row) must be valid
    valid <- !is.na(res$P[, 1])
    if (!valid[length(valid)]) {
      out$status[i] <- "invalid parameters at the estimate"
      next
    }
    Pv <- res$P[valid, , drop = FALSE]
    C <- t(apply(Pv, 1, cumsum))
    if (nrow(Pv) == 1) C <- matrix(C, 1)
    lu <- .fcr_ti_search(C, content, side, method, centre = nrow(C))
    if (any(!res$resolved[valid]) && side != "lower") {
      out$status[i] <- paste0("upper tail not resolved within max_count = ",
                              max_count)
      next
    }
    out$lower[i] <- lu[1]
    out$upper[i] <- lu[2]
    out$min_content[i] <- lu[3]
  }
  if (!is.null(pd)) {
    settings$bootstraps <- pd$requested
    settings$accepted <- pd$accepted
    settings$seed <- seed
  } else {
    settings$npoints <- npoints
  }
  attr(out, "tolerance") <- settings
  .fcr_ti_class(out)
}

# Rank of a positive semidefinite matrix from its eigenvalues
.fcr_psd_rank <- function(values) {
  sum(values > max(values, 0) * 1e-10)
}

# Mahalanobis distances with a possibly singular covariance matrix S: the
# distance in the span of S (its positive-eigenvalue directions). With rank
# 0 every point is at the centre (distance 0). stats::mahalanobis() fails on
# a singular S, as happens for bootstrap draws of a linear predictor that
# does not vary in some direction.
.fcr_mahalanobis_psd <- function(X, center, S) {
  X <- sweep(as.matrix(X), 2, center)
  e <- eigen(S, symmetric = TRUE)
  keep <- e$values > max(e$values, 0) * 1e-10
  if (!any(keep)) return(rep(0, nrow(X)))
  Z <- X %*% e$vectors[, keep, drop = FALSE]
  rowSums(sweep(Z^2, 2, e$values[keep], `/`))
}

.fcr_ti_class <- function(out) {
  class(out) <- c("flexCountReg_tolerance", "data.frame")
  out
}

#' @export
print.flexCountReg_tolerance <- function(x, ...) {
  s <- attr(x, "tolerance")
  cat("Tolerance intervals: content ", s$content, ", confidence ",
      s$confidence, ", ", s$side, ", ", s$scope, "\n", sep = "")
  cat("Method: ", s$method, " (", s$uncertainty, ")", sep = "")
  if (!is.null(s$accepted)) {
    cat(", ", s$accepted, " of ", s$bootstraps, " bootstrap refits",
        sep = "")
  }
  cat("\n")
  y <- x
  class(y) <- "data.frame"
  attr(y, "tolerance") <- NULL
  print(y)
  invisible(x)
}

# Is the exact common-rate Poisson method applicable?
.fcr_ti_exact_applies <- function(object, fd) {
  m <- object$model
  identical(m$family, "POISSON") && identical(colnames(fd$X), "(Intercept)") &&
    is.null(fd$aux$dis1) && is.null(fd$aux$dis2) &&
    is.null(fd$aux$underreport) && is.null(fd$aux$zero)
}

.fcr_ti_exact <- function(object, fd, est, content, confidence, side) {
  m <- object$model
  y <- m[["y"]]
  w <- m[["weights"]] %||% rep(1, length(y))
  if (any(w != round(w))) {
    stop("method = \"exact\" needs whole-number frequency weights.",
         call. = FALSE)
  }
  off_tr <- .fcr_countreg_fd(object, .fcr_training_rows(object))$offset
  S <- sum(w * y)
  E <- sum(w * exp(off_tr))
  alt <- switch(side, "two-sided" = "two.sided", upper = "less",
                lower = "greater")
  ci <- stats::poisson.test(S, T = E, conf.level = confidence,
                            alternative = alt)$conf.int
  e0 <- exp(fd$offset)
  lo <- hi <- rep(NA_real_, length(e0))
  if (side == "two-sided") {
    a <- (1 - content) / 2
    lo <- stats::qpois(a, e0 * ci[1])
    hi <- stats::qpois(1 - a, e0 * ci[2])
  } else if (side == "upper") {
    lo <- rep(0, length(e0))
    hi <- stats::qpois(content, e0 * ci[2])
  } else {
    lo <- stats::qpois(1 - content, e0 * ci[1])
    hi <- rep(Inf, length(e0))
  }
  # Content of a fixed interval is unimodal in the rate: its minimum over
  # the confidence interval is at an endpoint
  cont <- function(lam) stats::ppois(hi, e0 * lam) -
    stats::ppois(lo - 1, e0 * lam)
  rates <- ci[is.finite(ci) & ci > 0]
  minc <- if (length(rates)) do.call(pmin, lapply(rates, cont)) else NA
  data.frame(fit = est$response, lower = lo, upper = hi, min_content = minc,
             status = "ok", stringsAsFactors = FALSE)
}

# Design rows of the model components (mean, dis1, dis2, und, zero) and the
# coefficient names each uses.
.fcr_ti_components <- function(object, fd) {
  params <- get_params(object$model$family)
  n <- nrow(fd$X)
  comp <- list()
  comp$mean <- list(X = fd$X, names = colnames(fd$X))
  for (j in 1:2) {
    pn <- params[[j]]
    if (is.null(pn)) next
    X <- fd$aux[[paste0("dis", j)]]
    comp[[paste0("dis", j)]] <- if (is.null(X)) {
      list(X = matrix(1, n, 1), names = pn)
    } else {
      list(X = X, names = paste0(pn, ":", colnames(X)))
    }
  }
  if (!is.null(fd$aux$underreport)) {
    comp$und <- list(X = fd$aux$underreport,
                     names = paste0("Underreporting:",
                                    colnames(fd$aux$underreport)))
  }
  if (!is.null(fd$aux$zero)) {
    comp$zero <- list(X = fd$aux$zero,
                      names = paste0("Zero:", colnames(fd$aux$zero)))
  }
  list(comp = comp, offset = fd$offset)
}

# Linear map from the coefficients to the linear predictors of row i
# (without the offset)
.fcr_ti_A <- function(comps, i, coef_names) {
  A <- matrix(0, length(comps$comp), length(coef_names),
              dimnames = list(names(comps$comp), coef_names))
  for (nm in names(comps$comp)) {
    cc <- comps$comp[[nm]]
    A[nm, cc$names] <- cc$X[i, ]
  }
  A
}

# Linear predictors (one element per point) from points of the linear-
# predictor space; the offset is added to the mean
.fcr_ti_lp <- function(comps, eta, off) {
  nm <- names(comps$comp)
  col <- function(x) if (x %in% nm) eta[, match(x, nm)]
  list(mean = col("mean") + off, dis1 = col("dis1"), dis2 = col("dis2"),
       und = col("und"), zero = col("zero"))
}

# Deterministic points of the d-dimensional unit ball: the centre is last;
# half of the others on the boundary sphere, half inside (Halton-based)
.fcr_ball_points <- function(d, npoints) {
  if (d == 1) return(matrix(c(seq(-1, 1, length.out = npoints), 0), ncol = 1))
  nb <- ceiling(npoints / 2)
  ni <- npoints - nb
  if (d == 2) {
    t <- 2 * pi * (seq_len(nb) - 1) / nb
    bnd <- cbind(cos(t), sin(t))
  } else {
    z <- stats::qnorm(.fcr_halton(nb, d))
    bnd <- z / sqrt(rowSums(z^2))
  }
  h <- .fcr_halton(ni, d + 1)
  z <- stats::qnorm(h[, seq_len(d), drop = FALSE])
  inner <- z / sqrt(rowSums(z^2)) * h[, d + 1]^(1 / d)
  rbind(bnd, inner, rep(0, d))
}

# Interval search on cumulative probabilities C (points x counts 0..K).
# Returns c(L, U, min content over the points).
.fcr_ti_search <- function(C, p, side, method, centre) {
  K <- ncol(C) - 1
  Fm1 <- function(L) if (L == 0) rep(0, nrow(C)) else C[, L]
  colmin <- function(M) if (is.null(dim(M))) min(M) else apply(M, 2, min)
  qidx <- function(prob) apply(C, 1, function(f) {
    k <- which(f >= prob - 1e-12)
    if (length(k)) k[1] - 1 else NA_real_
  })
  if (side == "upper") {
    U <- max(qidx(p))
    return(c(0, U, if (is.na(U)) NA else min(C[, U + 1])))
  }
  if (side == "lower") {
    L <- min(qidx(1 - p))
    return(c(L, Inf, min(1 - Fm1(L))))
  }
  a <- (1 - p) / 2
  env <- c(min(qidx(a)), max(qidx(1 - a)))
  if (method == "envelope" || anyNA(env)) {
    if (anyNA(env)) return(c(NA, NA, NA))
    return(c(env, min(C[, env[2] + 1] - Fm1(env[1]))))
  }
  best <- c(env, min(C[, env[2] + 1] - Fm1(env[1])))
  L <- 0
  while (L <= K && Fm1(L)[centre] <= 1 - p) {
    D <- C[, (L + 1):(K + 1), drop = FALSE] - Fm1(L)
    cm <- colmin(D)
    hit <- which(cm >= p - 1e-12)
    if (length(hit)) {
      U <- L + hit[1] - 1
      c_min <- cm[hit[1]]
      if (U - L < best[2] - best[1] ||
          (U - L == best[2] - best[1] && c_min > best[3])) {
        best <- c(L, U, c_min)
      }
    }
    L <- L + 1
  }
  best
}
