#' Estimated smooth functions
#'
#' Evaluates the smooth terms (\code{s()}, \code{te()}, ...) of a model fitted
#' by \code{\link{countreg}()} or \code{\link{countreg_rp}()} on a grid over
#' the observed range of their variables, with pointwise and simultaneous
#' confidence bands.
#'
#' @param object a fitted model with smooth terms.
#' @param term optional term labels (as in \code{summary()}, e.g.
#'   \code{"s(lnaadt)"}); by default every smooth.
#' @param n number of grid points for one-variable smooths (the grid of a
#'   two-variable smooth has \code{round(sqrt(n) * 3)} points per variable).
#' @param level confidence level of the bands.
#' @param nsim number of draws for the simultaneous bands.
#' @param seed optional seed for those draws (the caller's random-number
#'   stream is left unchanged).
#'
#' @details Each smooth is shown on the scale of its linear predictor
#'   (log for the mean, the parameter's link for a distribution parameter),
#'   centered as fitted (the constant is in the intercept). Pointwise bands
#'   are \eqn{\hat f \pm z\,SE}, with \eqn{SE^2 = x' V x} from the fit's
#'   covariance matrix (Bayesian by default). Simultaneous bands use the
#'   95th (or \code{level}) percentile of
#'   \eqn{\max |f^* - \hat f| / SE} over the grid, for coefficient draws
#'   \eqn{f^*} from \eqn{N(\hat b, V)} (Ruppert, Wand and Carroll 2003), so
#'   they cover the whole curve with the stated probability.
#'
#' @returns A data frame with columns \code{parameter}, \code{term}, the
#'   smooth's variable(s), \code{fit}, \code{se}, \code{lwr}, \code{upr}
#'   (pointwise) and \code{slwr}, \code{supr} (simultaneous).
#' @references Ruppert, D., Wand, M.P. and Carroll, R.J. (2003).
#'   \emph{Semiparametric Regression}. Cambridge University Press.
#' @seealso \code{\link{plot_smooths}()}
#' @examples
#' \donttest{
#' data("washington_roads")
#' m <- countreg(Total_crashes ~ s(lnaadt) + offset(lnlength),
#'               data = washington_roads, family = "NB2")
#' head(smooth_estimates(m))
#' }
#' @export
smooth_estimates <- function(object, term = NULL, n = 100, level = 0.95,
                             nsim = 10000, seed = NULL) {
  info <- object$model$smooth_info
  if (!length(info)) stop("The model has no smooth terms.", call. = FALSE)
  if (!is.null(term)) {
    info <- Filter(function(s) s$term %in% term, info)
    if (!length(info)) {
      stop("No smooth term ", toString(term), " in the model.",
           call. = FALSE)
    }
  }
  d <- object$data
  rows <- object$model$rows_used
  if (!is.null(d) && !is.null(rows)) d <- d[rows, , drop = FALSE]
  V <- stats::vcov(object)
  b <- stats::coef(object)
  z <- stats::qnorm((1 + level) / 2)
  if (!is.null(seed)) {
    old <- if (exists(".Random.seed", envir = globalenv()))
      get(".Random.seed", envir = globalenv())
    on.exit(if (is.null(old)) rm(".Random.seed", envir = globalenv()) else
      assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  out <- lapply(info, function(s) {
    if (is.null(d)) {
      stop("The model was fitted with `keep_data = FALSE`; the observed ",
           "range of the smooth variables is not available.", call. = FALSE)
    }
    grid <- .fcr_smooth_grid(s, d, n)
    if (is.null(grid)) {
      warning("Skipping ", s$term, ": smooths of more than two numeric ",
              "variables are not evaluated.", call. = FALSE)
      return(NULL)
    }
    Xg <- mgcv::PredictMat(s$sm, grid)
    bi <- b[s$idx]
    Vi <- V[s$idx, s$idx, drop = FALSE]
    fit <- drop(Xg %*% bi)
    se <- sqrt(pmax(rowSums((Xg %*% Vi) * Xg), 0))
    e <- eigen((Vi + t(Vi)) / 2, symmetric = TRUE)
    L <- e$vectors %*% diag(sqrt(pmax(e$values, 0)), length(e$values))
    sims <- L %*% matrix(stats::rnorm(length(bi) * nsim), length(bi))
    dev <- abs(Xg %*% sims) / pmax(se, .Machine$double.eps)
    crit <- stats::quantile(apply(dev, 2, max), level, names = FALSE)
    keep <- intersect(c(s$vars, s$by), names(grid))
    data.frame(parameter = s$parameter, term = s$term,
               grid[keep], fit = fit, se = se,
               lwr = fit - z * se, upr = fit + z * se,
               slwr = fit - crit * se, supr = fit + crit * se,
               stringsAsFactors = FALSE, check.names = FALSE)
  })
  out <- Filter(Negate(is.null), out)
  if (!length(out)) return(data.frame())
  cols <- unique(unlist(lapply(out, names)))
  out <- lapply(out, function(o) {
    o[setdiff(cols, names(o))] <- NA
    o[cols]
  })
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

#' Plot estimated smooth functions
#'
#' Plots the smooth terms of a model with \code{\link{smooth_estimates}()}:
#' one panel per smooth, with the pointwise band shaded, the simultaneous
#' band dashed and the observed values as a rug. Two-variable smooths are
#' drawn as contour plots of the estimate.
#'
#' @inheritParams smooth_estimates
#' @param ... further arguments passed to \code{plot()}.
#'
#' @returns The data of \code{smooth_estimates()}, invisibly.
#' @examples
#' \donttest{
#' data("washington_roads")
#' m <- countreg(Total_crashes ~ s(lnaadt) + offset(lnlength),
#'               data = washington_roads, family = "NB2")
#' plot_smooths(m)
#' }
#' @export
plot_smooths <- function(object, term = NULL, n = 100, level = 0.95,
                         nsim = 10000, seed = 1, ...) {
  est <- smooth_estimates(object, term, n, level, nsim, seed)
  if (!nrow(est)) return(invisible(est))
  d <- object$data
  rows <- object$model$rows_used
  if (!is.null(d) && !is.null(rows)) d <- d[rows, , drop = FALSE]
  keys <- unique(est[c("parameter", "term")])
  info <- object$model$smooth_info
  for (r in seq_len(nrow(keys))) {
    e <- est[est$parameter == keys$parameter[r] & est$term == keys$term[r], ]
    s <- Filter(function(z) z$term == keys$term[r] &&
                  z$parameter == keys$parameter[r], info)[[1]]
    vars <- s$vars
    is_fac <- vapply(vars, function(v) !is.numeric(e[[v]]), logical(1))
    num <- vars[!is_fac]
    ylab <- if (keys$parameter[r] == "mean") "Effect on log mean" else
      paste("Effect on", keys$parameter[r])
    # A numeric `by` variable: the curve is that variable's coefficient
    if (!is.null(s$by) && is.numeric(d[[s$by]])) {
      ylab <- paste("Coefficient of", s$by)
    }
    main <- if (keys$parameter[r] == "mean") keys$term[r] else
      paste0(keys$parameter[r], ": ", keys$term[r])
    if (length(num) == 0L) {
      # Random effects or a Markov random field: estimates by level
      lev <- as.character(e[[vars[1]]])
      o <- order(e$fit)
      yy <- seq_along(o)
      graphics::plot(e$fit[o], yy, xlim = range(e$lwr, e$upr, finite = TRUE),
                     yaxt = "n", pch = 19, xlab = ylab, ylab = "",
                     main = main, ...)
      graphics::segments(e$lwr[o], yy, e$upr[o], yy)
      graphics::abline(v = 0, lty = 3)
      if (length(lev) <= 40L) {
        graphics::axis(2, at = yy, labels = lev[o], las = 1,
                       cex.axis = 0.7)
      }
    } else if (length(num) == 1L && any(is_fac)) {
      # Factor-smooth interaction: one curve per level
      fv <- vars[is_fac][1]
      levs <- unique(as.character(e[[fv]]))
      cols <- grDevices::hcl.colors(length(levs), "Dark 3")
      graphics::plot(e[[num]], e$fit, type = "n",
                     ylim = range(e$fit, finite = TRUE), xlab = num,
                     ylab = ylab, main = main, ...)
      for (j in seq_along(levs)) {
        sel <- as.character(e[[fv]]) == levs[j]
        graphics::lines(e[[num]][sel], e$fit[sel], col = cols[j], lwd = 2)
      }
      if (length(levs) <= 10L) {
        graphics::legend("topleft", legend = levs, col = cols, lwd = 2,
                         bty = "n", cex = 0.8, title = fv)
      }
    } else if (length(num) == 1L) {
      x <- e[[num]]
      graphics::plot(x, e$fit, type = "n",
                     ylim = range(e$slwr, e$supr, finite = TRUE),
                     xlab = num, ylab = ylab, main = main, ...)
      graphics::polygon(c(x, rev(x)), c(e$lwr, rev(e$upr)),
                        col = grDevices::adjustcolor("steelblue", 0.3),
                        border = NA)
      graphics::lines(x, e$slwr, lty = 2)
      graphics::lines(x, e$supr, lty = 2)
      graphics::lines(x, e$fit, lwd = 2)
      if (!is.null(d)) graphics::rug(d[[num]])
    } else {
      plot_surface(object, term = keys$term[r],
                   parameter = keys$parameter[r], est = e, ...)
    }
  }
  invisible(est)
}

# Evaluation grid of one smooth over the observed data `d`: numeric
# variables on a sequence over their range (one variable: n points; two:
# about 3 sqrt(n) per variable), factor variables at each level (random
# effects, Markov random fields, factor-smooth interactions). A numeric
# `by` variable is set to 1, so the term is the varying coefficient itself;
# a factor `by` is set to the smooth's own level. NULL for smooths of more
# than two numeric variables.
.fcr_smooth_grid <- function(s, d, n) {
  vars <- s$vars
  is_fac <- vapply(vars, function(v) is.factor(d[[v]]) ||
                     is.character(d[[v]]), logical(1))
  num <- vars[!is_fac]
  if (length(num) > 2L) return(NULL)
  vals <- list()
  m <- if (length(num) == 2L) round(sqrt(n) * 3) else n
  for (v in num) {
    r <- range(d[[v]], na.rm = TRUE)
    vals[[v]] <- seq(r[1], r[2], length.out = m)
  }
  for (v in vars[is_fac]) {
    lev <- if (is.factor(d[[v]])) levels(d[[v]]) else
      sort(unique(as.character(d[[v]])))
    vals[[v]] <- factor(lev, levels = lev)
  }
  grid <- expand.grid(vals[vars], KEEP.OUT.ATTRS = FALSE,
                      stringsAsFactors = FALSE)
  for (v in vars[is_fac]) {
    if (is.character(d[[v]])) grid[[v]] <- as.character(grid[[v]])
  }
  if (!is.null(s$by)) {
    bv <- d[[s$by]]
    if (is.numeric(bv)) {
      grid[[s$by]] <- 1
    } else {
      lev <- s$sm$by.level
      grid[[s$by]] <- if (is.factor(bv)) factor(lev, levels = levels(bv)) else
        lev
    }
  }
  grid
}
