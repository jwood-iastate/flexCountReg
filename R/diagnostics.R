# Probability predictions, randomized quantile residuals and rootograms for
# countreg() models.

# P(Y = k | x_i) for k = 0..max_count and each row of `data`, evaluated with
# the fitted family exactly as in the likelihood (mean parameter
# exp(X beta + offset) * p, the same distribution parameters and the same
# Halton draws). Returns an n x (max_count + 1) matrix.
.fcr_prob_matrix <- function(object, data, max_count) {
  if (inherits(object, "countreg_mix")) {
    return(.fcr_fm_prob_matrix(object, data, max_count))
  }
  if (!inherits(object, "countreg")) {
    stop("Probability predictions are available for countreg() models.",
         call. = FALSE)
  }
  parts <- .fcr_countreg_parts(object, data)
  P <- vapply(0:max_count, function(k) .fcr_parts_density(object, parts, k),
              numeric(length(parts$mu_link)))
  P <- matrix(P, ncol = max_count + 1,
              dimnames = list(NULL, as.character(0:max_count)))
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

# log P(Y <= q) and log P(Y > q) of the fitted count distribution (before
# any zero component) for each row of a .fcr_countreg_parts() result, at
# whole-number q (q < 0 gives log P(Y <= q) = -Inf). The distribution is the
# one of the likelihood: the same family, reporting transformation and
# parameters. Families with a closed-form CDF use it, so the cost does not
# depend on the counts; the others sum the PMF over 0..q row by row in one
# pass, with memory independent of the counts and a work budget checked
# before any evaluation.
.fcr_parts_cdf <- function(object, parts, q, max_work = 5e7) {
  n <- length(parts$mu_link)
  q <- rep_len(as.numeric(q), n)
  low <- rep(NA_real_, n)   # log P(Y <= q)
  upp <- rep(NA_real_, n)   # log P(Y > q)
  neg <- !is.na(q) & q < 0
  low[neg] <- -Inf
  upp[neg] <- 0
  todo <- which(!is.na(q) & q >= 0 & !is.na(parts$mu_link))
  if (!length(todo)) return(list(lower = low, upper = upp))
  fam <- .fcr_family_registry[[parts$family]]
  sub <- function(nm, i) if (!is.null(parts[[nm]])) parts[[nm]][i]
  mean_par <- parts$mu_link * parts$p_report
  if (!is.null(fam$cdf)) {
    args <- list(q[todo], mean_par[todo], sub("alpha", todo),
                 sub("sigma", todo))
    low[todo] <- suppressWarnings(do.call(fam$cdf, c(args, list(TRUE))))
    upp[todo] <- suppressWarnings(do.call(fam$cdf, c(args, list(FALSE))))
    return(list(lower = low, upper = upp))
  }
  work <- sum(q[todo] + 1)
  if (work > max_work) {
    stop("The ", parts$family, " fitted CDF has no closed form and would ",
         "need ", format(work, big.mark = ","), " probability evaluations ",
         "(the budget is ", format(max_work, big.mark = ","), ") for these ",
         "counts. Quantile residuals are not available for this model and ",
         "data.", call. = FALSE)
  }
  # One pass over k = 0..max(q): rows whose q is reached drop out
  acc <- rep(-Inf, n)
  active <- todo
  k <- 0
  while (length(active)) {
    sp <- parts
    for (nm in c("mu_link", "p_report", "alpha", "sigma")) {
      if (!is.null(parts[[nm]])) sp[[nm]] <- parts[[nm]][active]
    }
    lp <- .fcr_parts_density(object, sp, k, log = TRUE)
    acc[active] <- .fcr_logadd(acc[active], lp)
    done <- q[active] <= k
    active <- active[!done]
    k <- k + 1
  }
  low[todo] <- pmin(acc[todo], 0)
  upp[todo] <- .fcr_log1mexp(low[todo])
  list(lower = low, upper = upp)
}

# log P(Y <= y) of the fitted model, including a zero component: for
# zero-inflated models (structural zero probability pi)
#   F(y) = pi + (1 - pi) F_base(y),
# and for hurdle models (P(Y > 0) = z, base zero probability f0)
#   F(y) = 1 - z P_base(Y > y) / (1 - f0),  y >= 0,
# computed from the upper tail so it stays accurate when f0 is near 1.
.fcr_fitted_log_cdf <- function(object, parts, q) {
  cd <- .fcr_parts_cdf(object, parts, q)
  zp <- parts$zero_prob
  if (is.null(zp)) return(cd$lower)
  q <- rep_len(q, length(cd$lower))
  if (parts$zero_type == "inflated") {
    out <- .fcr_logadd(log(zp), log1p(-zp) + cd$lower)
    out[!is.na(q) & q < 0] <- -Inf
    return(out)
  }
  lf0 <- parts$lf0 %||% .fcr_parts_density(object, parts, 0, log = TRUE)
  log_s <- log(zp) + cd$upper - .fcr_log1mexp(lf0)
  out <- .fcr_log1mexp(pmin(log_s, 0))
  out[!is.na(q) & q < 0] <- -Inf
  out
}

# Randomized quantile residuals (Dunn and Smyth, 1996): for each training
# row, u ~ Uniform(F(y - 1), F(y)) under the fitted distribution, and the
# residual is qnorm(u). Only F(y - 1) and F(y) are evaluated, so the cost
# does not grow with a probability grid up to the largest count. `u` (one
# uniform per row) can be given to reproduce the randomization.
.fcr_quantile_residuals <- function(object, u = NULL) {
  d <- .fcr_training_rows(object)
  y <- as.numeric(object$model[["y"]])
  n <- length(y)
  if (is.null(u)) u <- stats::runif(n)
  parts <- .fcr_countreg_parts(object, d)
  upper <- exp(.fcr_fitted_log_cdf(object, parts, y))
  lower <- exp(.fcr_fitted_log_cdf(object, parts, y - 1))
  bad <- is.na(lower) | is.na(upper) | lower > upper + 1e-12 |
    upper > 1 + 1e-12
  if (any(bad)) {
    warning(sum(bad), " row(s) have invalid fitted probabilities; their ",
            "quantile residuals are NA.", call. = FALSE)
  }
  upper <- pmin(upper, 1)
  lower <- pmin(lower, upper)
  r <- stats::qnorm(lower + u * (upper - lower))
  r[bad] <- NA
  r
}

#' Rootogram for a count regression model
#'
#' Compares observed and expected count frequencies on a square-root scale
#' (Kleiber and Zeileis, 2016). Each bar shows the square root of the observed
#' frequency of a count, hanging from the curve of the square root of its
#' expected frequency under the fitted model, so bars that end above or below
#' zero show counts that the model over- or under-predicts. Rootograms are a
#' standard way to compare the count distributions this package offers.
#'
#' @param object a \code{countreg()} model.
#' @param newdata optional data (with the outcome) to evaluate the model on;
#'   defaults to the rows used in estimation.
#'   Without \code{newdata}, the observed and expected frequencies count each
#'   row with its frequency weight.
#' @param max_count the largest count shown. Defaults to the largest observed
#'   count. The cutoff only limits the display: the expected frequency of
#'   each shown count sums the fitted probabilities of every row, including
#'   rows whose outcome is above \code{max_count}.
#' @param overflow logical; add a last bin for \code{Y > max_count}, so the
#'   observed and expected frequencies both sum to the total. Without it,
#'   the tail is reported in the \code{"tail"} attribute and, when it is not
#'   negligible, under the plot.
#' @param style \code{"hanging"} (bars hang from the expected curve) or
#'   \code{"standing"} (bars stand on zero).
#' @param plot logical; draw the plot?
#' @param main,xlab,ylab plot labels.
#' @param ... further arguments passed to \code{\link[graphics]{plot}}.
#'
#' @returns Invisibly, a data frame with columns \code{count},
#'   \code{observed} and \code{expected} (frequencies), of class
#'   \code{flexCountReg_rootogram}; with \code{overflow = TRUE} its last row
#'   (count \code{max_count + 1}) is the bin \code{Y > max_count}.
#'   Attribute \code{"tail"} holds the observed and expected frequencies of
#'   \code{Y > max_count}. Rows with non-finite fitted means or probabilities are
#'   left out of both frequencies, with a warning.
#'
#' @references
#' Kleiber, C., & Zeileis, A. (2016). Visualizing count data regressions using
#' rootograms. \emph{The American Statistician}, 70(3), 296-303.
#'
#' @examples
#' \donttest{
#' data("washington_roads")
#' pois <- countreg(Total_crashes ~ lnaadt + lnlength,
#'                  data = washington_roads, family = "Poisson",
#'                  method = "BFGS")
#' nb2 <- countreg(Total_crashes ~ lnaadt + lnlength,
#'                 data = washington_roads, family = "NB2", method = "BFGS")
#' oldpar <- par(mfrow = c(1, 2))
#' rootogram(pois, main = "Poisson")
#' rootogram(nb2, main = "NB2")
#' par(oldpar)
#' }
#' @importFrom graphics plot lines points rect abline
#' @export
rootogram <- function(object, newdata = NULL, max_count = NULL,
                      style = c("hanging", "standing"), overflow = FALSE,
                      plot = TRUE, main = NULL, xlab = "Count",
                      ylab = "sqrt(Frequency)", ...) {
  style <- match.arg(style)
  d <- if (is.null(newdata)) .fcr_training_rows(object) else
    as.data.frame(newdata)
  y <- as.numeric(stats::model.response(stats::model.frame(
    object$formula, d, na.action = stats::na.pass)))
  if (is.null(max_count)) max_count <- max(y, na.rm = TRUE)
  if (!is.numeric(max_count) || length(max_count) != 1L ||
      !is.finite(max_count) || max_count < 0 ||
      max_count != round(max_count)) {
    stop("`max_count` must be a single whole number of at least 0.",
         call. = FALSE)
  }
  P <- .fcr_prob_matrix(object, d, max_count)
  # Frequency weights of the training rows count each row w times
  w <- if (is.null(newdata) && length(object$model$weights) == nrow(d))
    object$model$weights else rep(1, nrow(d))
  # The analysis rows: an observed outcome and finite fitted probabilities.
  # The same rows enter the observed and the expected frequencies, whatever
  # the display cutoff: a row whose outcome is above max_count still adds
  # its fitted probabilities to every displayed count.
  mu_hat <- suppressWarnings(stats::predict(object, newdata = d,
                                            type = "response"))
  bad_pred <- !is.na(y) & (!stats::complete.cases(P) | !is.finite(mu_hat))
  if (any(bad_pred)) {
    warning(sum(bad_pred), " row(s) with non-finite fitted means or ",
            "probabilities are left out of the observed and expected frequencies.",
            call. = FALSE)
  }
  valid <- !is.na(y) & !bad_pred
  shown <- valid & y <= max_count
  counts <- 0:max_count
  obs <- as.vector(tapply(w[shown], factor(y[shown], levels = counts), sum))
  obs[is.na(obs)] <- 0
  expd <- colSums(P[valid, , drop = FALSE] * w[valid])
  # The tail Y > max_count, which the displayed counts leave out
  tail_obs <- sum(w[valid & y > max_count])
  tail_exp <- max(0, sum(w[valid] * (1 - rowSums(P[valid, , drop = FALSE]))))
  out <- data.frame(count = counts, observed = obs, expected = unname(expd))
  if (overflow) {
    out <- rbind(out, data.frame(count = max_count + 1, observed = tail_obs,
                                 expected = tail_exp))
  }
  attr(out, "tail") <- c(observed = tail_obs, expected = tail_exp)
  attr(out, "overflow") <- overflow
  class(out) <- c("flexCountReg_rootogram", "data.frame")
  counts <- out$count
  obs <- out$observed
  expd <- out$expected

  if (plot) {
    so <- sqrt(obs)
    se <- sqrt(expd)
    if (style == "hanging") {
      top <- se; bottom <- se - so
    } else {
      top <- so; bottom <- rep(0, length(so))
    }
    if (is.null(main)) {
      main <- paste("Rootogram:", paste(unique(object$model$family %||%
                                                  class(object)[1]),
                                           collapse = " + "))
    }
    graphics::plot(range(counts) + c(-0.5, 0.5),
                   range(c(0, top, bottom, se)), type = "n", main = main,
                   xlab = xlab, ylab = ylab, ...)
    graphics::rect(counts - 0.4, bottom, counts + 0.4, top, col = "grey80",
                   border = "grey40")
    graphics::abline(h = 0, lty = 2)
    graphics::lines(counts, se, col = "firebrick", lwd = 2)
    graphics::points(counts, se, col = "firebrick", pch = 19)
    if (overflow) {
      graphics::mtext(paste0("last bar: Y > ", max_count), side = 1,
                      line = 4, cex = 0.8)
    } else if (tail_obs > 0 || tail_exp > 0.005 * sum(w[valid])) {
      graphics::mtext(sprintf("Not shown: Y > %d (observed %s, expected %s)",
                              max_count, format(tail_obs),
                              format(signif(tail_exp, 3))),
                      side = 1, line = 4, cex = 0.8)
    }
  }
  invisible(out)
}
