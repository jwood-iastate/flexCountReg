# Diagnostics and further plots for smooth terms.

# The training rows of a fit (NULL with keep_data = FALSE).
.fcr_smooth_rows <- function(object) {
  d <- object$data
  rows <- object$model$rows_used
  if (!is.null(d) && !is.null(rows)) d <- d[rows, , drop = FALSE]
  d
}

.fcr_smooth_info_of <- function(object, term = NULL) {
  info <- object$model$smooth_info
  if (!length(info)) stop("The model has no smooth terms.", call. = FALSE)
  if (!is.null(term)) {
    info <- Filter(function(s) s$term %in% term, info)
    if (!length(info)) {
      stop("No smooth term ", toString(term), " in the model.",
           call. = FALSE)
    }
  }
  info
}

# Residuals for diagnostics: randomized quantile residuals where the model
# provides them (N(0, 1) for a correct count model), else Pearson, else
# response residuals. Returns list(r, type).
.fcr_check_residuals <- function(object, seed) {
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv())
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
      assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  for (type in c("quantile", "pearson", "response")) {
    r <- tryCatch(suppressWarnings(stats::residuals(object, type = type)),
                  error = function(e) NULL)
    if (!is.null(r) && sum(is.finite(r)) > 2) return(list(r = r, type = type))
  }
  stop("Residuals are not available for this model.", call. = FALSE)
}

#' Check the basis dimension of smooth terms
#'
#' For each smooth, compares the residuals of neighboring observations
#' (ordered by the smooth's covariates) with the residual variance. If the
#' basis dimension \code{k} is too small, the residuals keep a pattern in
#' the covariate, so neighbors are more alike than the overall variance
#' suggests and the k-index falls below 1.
#'
#' @param object a fitted model with smooth terms.
#' @param term optional term labels; by default every smooth.
#' @param nperm number of random permutations of the residuals for the
#'   p-value.
#' @param seed seed for the residual randomization and the permutations
#'   (the caller's random-number stream is left unchanged).
#'
#' @details The k-index is \eqn{\frac{1}{2}\overline{(r_i - r_{(i)})^2} /
#'   \mathrm{var}(r)}, where \eqn{r_{(i)}} is the residual of the nearest
#'   neighbor of observation \eqn{i} in the smooth's (scaled) covariates;
#'   for one covariate, the next observation in its order. The residuals are
#'   the randomized quantile residuals of the model, which are N(0, 1) for
#'   a correct count model whatever the family (Pearson residuals where they
#'   are not available). The p-value is the share of residual permutations
#'   with a smaller k-index: a small p-value with a k-index below 1 suggests
#'   increasing \code{k} (refit and compare the EDF and the curve). Smooths
#'   of factors only (random effects, Markov random fields) are not checked.
#'   This follows \code{mgcv::k.check()}, adapted to count residuals.
#'
#' @returns A data frame with \code{parameter}, \code{term}, \code{k} (the
#'   basis dimension after the identifiability constraint), \code{edf},
#'   \code{k_index}, \code{p_value} and the residual type.
#' @examples
#' \donttest{
#' data("washington_roads")
#' m <- countreg(Total_crashes ~ s(lnaadt, k = 5) + offset(lnlength),
#'               data = washington_roads, family = "NB2")
#' smooth_check(m)
#' }
#' @export
smooth_check <- function(object, term = NULL, nperm = 400, seed = 1) {
  info <- .fcr_smooth_info_of(object, term)
  d <- .fcr_smooth_rows(object)
  if (is.null(d)) {
    stop("The model was fitted with `keep_data = FALSE`.", call. = FALSE)
  }
  res <- .fcr_check_residuals(object, seed)
  r <- as.numeric(res$r)
  st <- object$model$smooth_table
  kidx <- function(rr, ord, nn) {
    if (!is.null(ord)) {
      rs <- rr[ord]
      mean(diff(rs)^2) / 2 / stats::var(rr)
    } else {
      mean((rr - rr[nn])^2) / 2 / stats::var(rr)
    }
  }
  out <- lapply(info, function(s) {
    vars <- s$vars
    num <- vars[vapply(vars, function(v) is.numeric(d[[v]]), logical(1))]
    edf <- st$edf[st$term == s$term & st$parameter == s$parameter][1]
    row <- data.frame(parameter = s$parameter, term = s$term,
                      k = length(s$idx), edf = edf, k_index = NA_real_,
                      p_value = NA_real_, residuals = res$type,
                      stringsAsFactors = FALSE)
    if (!length(num)) return(row)
    ok <- is.finite(r) & stats::complete.cases(d[num])
    rr <- r[ok]
    Xc <- scale(as.matrix(d[ok, num, drop = FALSE]))
    Xc[!is.finite(Xc)] <- 0
    ord <- nn <- NULL
    if (length(num) == 1L) {
      ord <- order(Xc[, 1])
    } else {
      # nearest neighbor in the scaled covariates (blocks of rows)
      nn <- integer(nrow(Xc))
      for (st0 in seq(1, nrow(Xc), by = 500)) {
        i <- st0:min(nrow(Xc), st0 + 499)
        D <- outer(rowSums(Xc[i, , drop = FALSE]^2), rowSums(Xc^2), "+") -
          2 * Xc[i, , drop = FALSE] %*% t(Xc)
        D[cbind(seq_along(i), i)] <- Inf
        nn[i] <- max.col(-D, ties.method = "first")
      }
    }
    k0 <- kidx(rr, ord, nn)
    perm <- .fcr_with_seed(seed, vapply(seq_len(nperm), function(b)
      kidx(sample(rr), ord, nn), 0))
    row$k_index <- k0
    row$p_value <- mean(perm <= k0)
    row
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

#' Concurvity of smooth terms
#'
#' Measures how much each smooth can be reproduced by the other terms of
#' its linear predictor: the nonlinear analogue of collinearity. High
#' concurvity makes a smooth's shape poorly determined and its standard
#' errors large.
#'
#' @param object a fitted model with smooth terms.
#'
#' @details For a smooth \eqn{f_j = X_j b_j} at the data rows,
#'   \code{observed} is the share of \eqn{\|f_j\|^2} explained by the
#'   other columns of the same linear predictor (the parametric terms and
#'   the other smooths), and \code{worst} is the largest squared canonical
#'   correlation between the smooth's basis and those columns, the
#'   worst case over all possible estimates. Both are between 0 and 1;
#'   values near 1 indicate a smooth that the rest of the model can mimic.
#'   This follows \code{mgcv::concurvity()} in spirit.
#'
#' @returns A data frame with \code{parameter}, \code{term},
#'   \code{observed} and \code{worst}.
#' @export
concurvity_fcr <- function(object) {
  info <- .fcr_smooth_info_of(object)
  d <- .fcr_smooth_rows(object)
  if (is.null(d)) {
    stop("The model was fitted with `keep_data = FALSE`.", call. = FALSE)
  }
  fd <- forge_new_data(.fcr_blueprint(object), d)
  b <- stats::coef(object)
  out <- lapply(info, function(s) {
    M <- if (s$base == "main") fd$X else fd$aux[[s$base]]
    if (is.null(M)) return(NULL)
    if (identical(s$base, "mixing") && !is.null(object$model$unit)) {
      M <- M[!duplicated(object$model$unit), , drop = FALSE]
    }
    lab <- paste0(s$term, ".", seq_along(s$idx))
    cols <- match(lab, colnames(M))
    ok <- stats::complete.cases(M)
    Xj <- M[ok, cols, drop = FALSE]
    Xo <- M[ok, -cols, drop = FALSE]
    fj <- drop(Xj %*% b[s$idx])
    obs <- if (ncol(Xo)) {
      fit <- stats::lm.fit(Xo, fj)$fitted.values
      sum(fit^2) / sum(fj^2)
    } else 0
    worst <- if (ncol(Xo)) {
      Qj <- qr.Q(qr(Xj))
      Qo <- qr.Q(qr(Xo))
      max(svd(crossprod(Qj, Qo), nu = 0, nv = 0)$d)^2
    } else 0
    data.frame(parameter = s$parameter, term = s$term, observed = obs,
               worst = min(worst, 1), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), out))
  rownames(out) <- NULL
  out
}

#' Derivatives of smooth functions
#'
#' The first derivative \eqn{df/dx} of each one-variable smooth, with
#' pointwise and simultaneous confidence bands: where the band excludes 0,
#' the smooth is increasing (or decreasing).
#'
#' @inheritParams smooth_estimates
#' @param h step for the central differences of the basis (default
#'   \code{1e-4} times the range of the variable).
#'
#' @details The derivative is exact for the linear predictor up to the
#'   central-difference error of the basis: with
#'   \eqn{D = (X(x + h) - X(x - h)) / 2h}, \eqn{df/dx = D b} and its
#'   standard error is \eqn{\sqrt{D V D'}}. Simultaneous bands are computed
#'   as in \code{\link{smooth_estimates}()}.
#'
#' @returns A data frame like that of \code{smooth_estimates()}, with the
#'   derivative in \code{fit}.
#' @seealso \code{\link{plot_derivatives}()}
#' @export
smooth_derivatives <- function(object, term = NULL, n = 100, level = 0.95,
                               nsim = 10000, seed = NULL, h = NULL) {
  info <- .fcr_smooth_info_of(object, term)
  d <- .fcr_smooth_rows(object)
  if (is.null(d)) {
    stop("The model was fitted with `keep_data = FALSE`.", call. = FALSE)
  }
  V <- stats::vcov(object)
  b <- stats::coef(object)
  z <- stats::qnorm((1 + level) / 2)
  out <- lapply(info, function(s) {
    grid <- .fcr_smooth_grid(s, d, n)
    num <- s$vars[vapply(s$vars, function(v) is.numeric(d[[v]]),
                         logical(1))]
    if (is.null(grid) || length(num) != 1L || length(s$vars) != 1L) {
      return(NULL)
    }
    v <- num
    hh <- h %||% (1e-4 * diff(range(d[[v]], na.rm = TRUE)))
    gp <- grid
    gm <- grid
    gp[[v]] <- gp[[v]] + hh
    gm[[v]] <- gm[[v]] - hh
    D <- (mgcv::PredictMat(s$sm, gp) - mgcv::PredictMat(s$sm, gm)) /
      (2 * hh)
    bi <- b[s$idx]
    Vi <- V[s$idx, s$idx, drop = FALSE]
    fit <- drop(D %*% bi)
    se <- sqrt(pmax(rowSums((D %*% Vi) * D), 0))
    e <- eigen((Vi + t(Vi)) / 2, symmetric = TRUE)
    L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values))
    crit <- .fcr_with_seed(seed %||% 1, {
      sims <- L %*% matrix(stats::rnorm(length(bi) * nsim), length(bi))
      stats::quantile(apply(abs(D %*% sims) / pmax(se, 1e-300), 2, max),
                      level, names = FALSE)
    })
    data.frame(parameter = s$parameter, term = s$term,
               grid[intersect(c(s$vars, s$by), names(grid))],
               fit = fit, se = se, lwr = fit - z * se, upr = fit + z * se,
               slwr = fit - crit * se, supr = fit + crit * se,
               stringsAsFactors = FALSE, check.names = FALSE)
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) {
    stop("No one-variable smooth to differentiate.", call. = FALSE)
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

#' Plot derivatives of smooth functions
#'
#' @inheritParams smooth_derivatives
#' @param ... further arguments passed to \code{plot()}.
#' @returns The data of \code{smooth_derivatives()}, invisibly.
#' @export
plot_derivatives <- function(object, term = NULL, n = 100, level = 0.95,
                             nsim = 10000, seed = 1, ...) {
  est <- smooth_derivatives(object, term, n, level, nsim, seed)
  keys <- unique(est[c("parameter", "term")])
  for (r in seq_len(nrow(keys))) {
    e <- est[est$parameter == keys$parameter[r] & est$term == keys$term[r], ]
    v <- setdiff(names(e), c("parameter", "term", "fit", "se", "lwr", "upr",
                             "slwr", "supr"))[1]
    x <- e[[v]]
    main <- if (keys$parameter[r] == "mean") keys$term[r] else
      paste0(keys$parameter[r], ": ", keys$term[r])
    graphics::plot(x, e$fit, type = "n",
                   ylim = range(e$slwr, e$supr, 0, finite = TRUE),
                   xlab = v, ylab = paste0("d f / d ", v),
                   main = paste("Derivative of", main), ...)
    graphics::polygon(c(x, rev(x)), c(e$lwr, rev(e$upr)),
                      col = grDevices::adjustcolor("steelblue", 0.3),
                      border = NA)
    graphics::lines(x, e$slwr, lty = 2)
    graphics::lines(x, e$supr, lty = 2)
    graphics::lines(x, e$fit, lwd = 2)
    graphics::abline(h = 0, lty = 3)
  }
  invisible(est)
}

#' Plot a two-variable smooth as a surface
#'
#' Contours (over a color image) of a smooth of two variables, such as
#' \code{te(lnaadt, lnlength)}, masked where there are no data nearby.
#'
#' @param object a fitted model with a smooth of two variables.
#' @param term the smooth's label (default: the first two-variable smooth).
#' @param parameter optional parameter label, when the same smooth appears
#'   in several components.
#' @param n grid size (as in \code{\link{smooth_estimates}()}).
#' @param too_far grid points farther than this (on the unit square of the
#'   two variables) from every observation are not shown
#'   (\code{mgcv::exclude.too.far()}); \code{0} shows the whole grid.
#' @param est optional output of \code{smooth_estimates()} for the term.
#' @param ... further arguments passed to \code{image()}.
#' @returns The plotted grid, invisibly.
#' @export
plot_surface <- function(object, term = NULL, parameter = NULL, n = 400,
                         too_far = 0.1, est = NULL, ...) {
  info <- object$model$smooth_info
  d <- .fcr_smooth_rows(object)
  two <- Filter(function(s) sum(vapply(s$vars, function(v)
    is.numeric(d[[v]]), logical(1))) == 2L, info)
  if (!is.null(term)) two <- Filter(function(s) s$term == term, two)
  if (!is.null(parameter)) {
    two <- Filter(function(s) s$parameter == parameter, two)
  }
  if (!length(two)) {
    stop("No smooth of two variables in the model.", call. = FALSE)
  }
  s <- two[[1]]
  if (is.null(est)) {
    est <- smooth_estimates(object, term = s$term, n = n, nsim = 200,
                            seed = 1)
    est <- est[est$parameter == s$parameter, ]
  }
  v <- s$vars
  u1 <- sort(unique(est[[v[1]]]))
  u2 <- sort(unique(est[[v[2]]]))
  z <- matrix(est$fit, length(u1), length(u2))
  if (too_far > 0) {
    far <- mgcv::exclude.too.far(est[[v[1]]], est[[v[2]]], d[[v[1]]],
                                 d[[v[2]]], dist = too_far)
    z[matrix(far, length(u1), length(u2))] <- NA
  }
  main <- if (s$parameter == "mean") s$term else
    paste0(s$parameter, ": ", s$term)
  graphics::image(u1, u2, z, col = grDevices::hcl.colors(30, "Blue-Red 3"),
                  xlab = v[1], ylab = v[2], main = main, ...)
  graphics::contour(u1, u2, z, add = TRUE)
  graphics::points(d[[v[1]]], d[[v[2]]], pch = ".", col = "grey30")
  invisible(est)
}

#' Calibration of a count model by a covariate
#'
#' Observed and predicted mean counts, probabilities of zero and variances
#' within groups (quantiles) of a covariate. A smooth can fit the mean well
#' and still miss the zeros or the dispersion; this table shows where.
#'
#' @param object a fitted model.
#' @param var the name of a numeric covariate (a column of the data).
#' @param groups number of quantile groups.
#'
#' @returns A data frame with one row per group: the covariate's range,
#'   \code{n}, observed and predicted mean, observed and predicted
#'   \eqn{P(Y = 0)} and observed and predicted variance. Predicted
#'   probabilities and variances are averaged over the rows of each group
#'   (the predicted variance includes the spread of the means within the
#'   group) and are \code{NA} for models without probability predictions.
#' @export
calibration_by <- function(object, var, groups = 10) {
  d <- .fcr_smooth_rows(object) %||% .fcr_training_rows(object)
  if (!var %in% names(d) || !is.numeric(d[[var]])) {
    stop("`var` must name a numeric column of the model's data.",
         call. = FALSE)
  }
  y <- as.numeric(object$model$y)
  mu <- as.numeric(stats::predict(object, newdata = d))
  P <- tryCatch(suppressWarnings(stats::predict(
    object, newdata = d, type = "prob", max_count = max(y) + 50,
    overflow = TRUE)), error = function(e) NULL)
  p0 <- if (!is.null(P)) P[, 1] else rep(NA_real_, length(y))
  v <- if (!is.null(P)) {
    k <- c(0:(ncol(P) - 2), max(y) + 51)
    rowSums(P * outer(rep(1, nrow(P)), k^2)) - mu^2
  } else rep(NA_real_, length(y))
  br <- unique(stats::quantile(d[[var]], seq(0, 1, length.out = groups + 1),
                               na.rm = TRUE))
  g <- cut(d[[var]], br, include.lowest = TRUE)
  out <- do.call(rbind, lapply(levels(g), function(l) {
    i <- which(g == l)
    data.frame(group = l, n = length(i), observed_mean = mean(y[i]),
               predicted_mean = mean(mu[i]), observed_p0 = mean(y[i] == 0),
               predicted_p0 = mean(p0[i]), observed_var = stats::var(y[i]),
               predicted_var = mean(v[i]) + stats::var(mu[i]),
               stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}
