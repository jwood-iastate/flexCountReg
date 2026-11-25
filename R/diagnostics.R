# Probability predictions, randomized quantile residuals and rootograms for
# countreg() models.

# P(Y = k | x_i) for the counts and each row of `data`, evaluated with
# the fitted family exactly as in the likelihood (mean parameter
# exp(X beta + offset) * p, the same distribution parameters and the same
# Halton draws). Returns an n x length(counts) matrix. Other models:
# .fcr_count_pmf().
.fcr_prob_matrix <- function(object, data, max_count,
                             counts = seq.int(0, max_count)) {
  if (!inherits(object, "countreg")) {
    P <- .fcr_count_pmf(object, data, counts)
    dimnames(P) <- list(NULL, as.character(counts))
    return(P)
  }
  parts <- .fcr_countreg_parts(object, data)
  P <- vapply(counts, function(k) .fcr_parts_density(object, parts, k),
              numeric(length(parts$mu_link)))
  P <- matrix(P, ncol = length(counts),
              dimnames = list(NULL, as.character(counts)))
  # A row with a missing value in any component (mean, distribution
  # parameters, reporting or zero model) has no distribution: its
  # probabilities are NA, not 0 (which made the overflow column 1)
  miss <- is.na(parts$mu_link)
  for (nm in c("p_report", "alpha", "sigma", "zero_prob")) {
    if (!is.null(parts[[nm]])) miss <- miss | is.na(parts[[nm]])
  }
  if (!is.null(parts$zero_prob)) {
    zp <- parts$zero_prob
    k0 <- which(counts == 0)
    if (parts$zero_type == "inflated") {
      f0 <- P[, k0]
      P <- P * (1 - zp)
      P[, k0] <- zp + (1 - zp) * f0
    } else {
      lf0 <- parts$lf0 %||% .fcr_parts_density(object, parts, 0, log = TRUE)
      P <- P * (zp / -expm1(lf0))
      P[, k0] <- 1 - zp
    }
  }
  P[miss, ] <- NA_real_
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
  fam <- .fcr_family_or_null(parts$family)
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
  .fcr_fitted_log_tails(object, parts, q)$lower
}

# log P(Y <= q) (lower) and log P(Y > q) (upper) of the fitted distribution,
# including any zero component. Both tails are kept on the log scale, so a
# probability near 1 in one tail is not lost to rounding in the other.
.fcr_fitted_log_tails <- function(object, parts, q) {
  cd <- .fcr_parts_cdf(object, parts, q)
  zp <- parts$zero_prob
  if (is.null(zp)) return(cd)
  q <- rep_len(q, length(cd$lower))
  neg <- !is.na(q) & q < 0
  if (parts$zero_type == "inflated") {
    low <- .fcr_logadd(log(zp), log1p(-zp) + cd$lower)
    upp <- log1p(-zp) + cd$upper
  } else {
    lf0 <- parts$lf0 %||% .fcr_parts_density(object, parts, 0, log = TRUE)
    upp <- pmin(log(zp) + cd$upper - .fcr_log1mexp(lf0), 0)
    low <- .fcr_log1mexp(upp)
  }
  low[neg] <- -Inf
  upp[neg] <- 0
  list(lower = low, upper = upp)
}

# Randomized quantile residuals (Dunn and Smyth, 1996): for each training
# row, u ~ Uniform(F(y - 1), F(y)) under the fitted distribution, and the
# residual is qnorm(u). Only F(y - 1) and F(y) are evaluated, so the cost
# does not grow with a probability grid up to the largest count. `u` (one
# uniform per row) can be given to reproduce the randomization.
#
# Everything stays on the log scale. With lp = log P(Y = y),
#   log U     = logadd(log F(y - 1), log(u) + lp)        (lower half)
#   log(1 - U) = logadd(log S(y), log1p(-u) + lp)        (upper half)
# and the residual is qnorm() of that log probability in the matching tail,
# so an outcome far in either tail gives a large finite residual (for
# Poisson(1) and y = 25, F(24) and F(25) both round to 1 and the
# exponentiated form gave Inf).
.fcr_quantile_residuals <- function(object, u = NULL) {
  d <- .fcr_training_rows(object)
  y <- as.numeric(object$model[["y"]])
  n <- length(y)
  if (is.null(u)) u <- stats::runif(n)
  parts <- .fcr_countreg_parts(object, d)
  t_hi <- .fcr_fitted_log_tails(object, parts, y)       # F(y), S(y)
  t_lo <- .fcr_fitted_log_tails(object, parts, y - 1)   # F(y - 1), S(y - 1)
  .fcr_qres_from_log_tails(t_lo$lower, t_hi$lower, t_lo$upper, t_hi$upper, u)
}

# Randomized quantile residuals from log F(y - 1), log F(y), log S(y - 1)
# and log S(y) (see .fcr_quantile_residuals()).
.fcr_qres_from_log_tails <- function(lF_lo, lF_hi, lS_lo, lS_hi, u) {
  bad <- is.na(lF_lo) | is.na(lF_hi) | lF_lo > lF_hi + 1e-12 |
    lF_hi > 1e-12
  if (any(bad)) {
    warning(sum(bad), " row(s) have invalid fitted probabilities; their ",
            "quantile residuals are NA.", call. = FALSE)
  }
  r <- rep(NA_real_, length(u))
  up <- !bad & lF_lo > log(0.5) & !is.na(lS_lo) & !is.na(lS_hi)
  lo <- !bad & !up
  if (any(lo)) {
    a <- pmin(lF_lo[lo], lF_hi[lo])
    lp <- lF_hi[lo] + .fcr_log1mexp(a - lF_hi[lo])
    logU <- .fcr_logadd(a, log(u[lo]) + lp)
    r[lo] <- stats::qnorm(pmin(logU, 0), log.p = TRUE)
  }
  if (any(up)) {
    b <- pmin(lS_hi[up], lS_lo[up])
    lp <- lS_lo[up] + .fcr_log1mexp(b - lS_lo[up])
    log1mU <- .fcr_logadd(b, log1p(-u[up]) + lp)
    r[up] <- stats::qnorm(pmin(log1mU, 0), lower.tail = FALSE, log.p = TRUE)
  }
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
#' @param object a fitted model of any estimator except
#'   \code{countreg_alloc()}. For a multivariate model
#'   (\code{countreg_mv()}), one rootogram is made for each outcome.
#' @param newdata optional data (with the outcome) to evaluate the model on;
#'   defaults to the rows used in estimation.
#'   Without \code{newdata}, the observed and expected frequencies count each
#'   row with its frequency weight.
#' @param max_count the largest count shown. Defaults to the largest observed
#'   count. The cutoff only limits the display: the expected frequency of
#'   each shown count sums the fitted probabilities of every row, including
#'   rows whose outcome is above \code{max_count}. For a multivariate
#'   model, one number for every outcome or a vector named by outcome.
#' @param overflow logical; add a last bin for \code{Y > max_count}, so the
#'   observed and expected frequencies both sum to the total. Without it,
#'   the tail is reported in the \code{"tail"} attribute and, when it is not
#'   negligible, under the plot.
#' @param style \code{"hanging"} (bars hang from the expected curve) or
#'   \code{"standing"} (bars stand on zero).
#' @param plot logical; draw the plot?
#' @param main,xlab,ylab plot labels. For a multivariate model, \code{main}
#'   can have one title per outcome.
#' @param outcome for \code{countreg_mv()} models: the name(s) of the
#'   outcomes to show (default: all).
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
#'   For a multivariate model, a list of these data frames named by outcome
#'   (the data frame itself when \code{outcome} names a single one).
#'
#' @details
#' The rootogram of an outcome of a multivariate model compares its observed
#' counts with its own (marginal) fitted distribution, with any random
#' parameters integrated out; it does not assess the dependence between the
#' outcomes. The plots are drawn side by side when the graphics device shows
#' a single plot, and otherwise fill the layout set with \code{par(mfrow)}.
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
#'
#' # One rootogram per outcome of a multivariate model
#' mv <- countreg_mv(cbind(Animal, Rollover) ~ lnaadt + lnlength,
#'                   data = washington_roads)
#' rootogram(mv)
#' rootogram(mv, outcome = "Animal")
#' }
#' @importFrom graphics plot lines points rect abline
#' @export
rootogram <- function(object, newdata = NULL, max_count = NULL,
                      style = c("hanging", "standing"), overflow = FALSE,
                      plot = TRUE, main = NULL, xlab = "Count",
                      ylab = "sqrt(Frequency)", outcome = NULL, ...) {
  style <- match.arg(style)
  if (inherits(object, "countreg_mv")) {
    return(.fcr_mv_rootogram(object, newdata, max_count, style, overflow,
                             plot, main, xlab, ylab, outcome, ...))
  }
  if (inherits(object, "countreg_alloc")) {
    stop("rootogram() is not available for countreg_alloc() models.",
         call. = FALSE)
  }
  if (!is.null(outcome)) {
    stop("`outcome` is for countreg_mv() models.", call. = FALSE)
  }
  d <- if (is.null(newdata)) .fcr_training_rows(object) else
    as.data.frame(newdata)
  y <- .fcr_response(object, d)
  if (is.null(y)) {
    stop("`newdata` must contain the model's outcome.", call. = FALSE)
  }
  # Frequency weights of the training rows count each row w times
  w <- if (is.null(newdata) && length(object$model$weights) == nrow(d))
    object$model$weights else rep(1, nrow(d))
  mu_hat <- suppressWarnings(stats::predict(object, newdata = d,
                                            type = "response"))
  if (is.null(main)) {
    main <- paste("Rootogram:", paste(unique(object$model$family %||%
                                                class(object)[1]),
                                         collapse = " + "))
  }
  invisible(.fcr_rootogram_one(
    y, function(mc) .fcr_prob_matrix(object, d, mc), w, mu_hat, max_count,
    style, overflow, plot, main, xlab, ylab, ...))
}

# One rootogram per outcome of a countreg_mv() model.
.fcr_mv_rootogram <- function(object, newdata, max_count, style, overflow,
                              plot, main, xlab, ylab, outcome, ...) {
  outs <- outcome %||% object$outcomes
  bad <- setdiff(outs, object$outcomes)
  if (length(bad) || !length(outs)) {
    stop("`outcome` must be among ", toString(object$outcomes), ".",
         call. = FALSE)
  }
  outs <- unique(outs)
  if (!is.null(names(max_count))) {
    miss <- setdiff(outs, names(max_count))
    if (length(miss)) {
      stop("`max_count` has no value for ", toString(miss), ".",
           call. = FALSE)
    }
  } else if (length(max_count) > 1L) {
    stop("`max_count` must be a single number or a vector named by outcome.",
         call. = FALSE)
  }
  if (!is.null(main) && !length(main) %in% c(1L, length(outs))) {
    stop("`main` must have one title, or one per outcome.", call. = FALSE)
  }
  d <- if (is.null(newdata)) object$data_used else as.data.frame(newdata)
  w <- if (is.null(newdata) && length(object$weights) == nrow(d))
    object$weights else rep(1, nrow(d))
  mu <- suppressWarnings(stats::predict(object, newdata = d))
  # The counts of every outcome, before anything is drawn
  ys <- lapply(outs, function(o) {
    j <- match(o, object$outcomes)
    if (is.null(newdata)) return(as.numeric(object$y[, j]))
    y <- .fcr_response(list(formula = object$formulas[[j]]), d)
    if (is.null(y)) {
      stop("`newdata` must contain the outcome ", o, ".", call. = FALSE)
    }
    y
  })
  if (plot && length(outs) > 1L && all(graphics::par("mfrow") == 1)) {
    op <- graphics::par(mfrow = rev(grDevices::n2mfrow(length(outs))))
    on.exit(graphics::par(op))
  }
  E <- .fcr_mv_design(object, d)
  res <- lapply(seq_along(outs), function(i) {
    j <- match(outs[i], object$outcomes)
    pmf <- function(mc) .fcr_mv_integrate(object, E, d, function(adds)
      .fcr_mv_pmf(object$coefficients, E$eqs[[j]], adds[[j]], 0:mc))
    title <- if (is.null(main)) {
      paste0("Rootogram: ", outs[i], " (",
             .fcr_family_label(object$family[[j]]), ")")
    } else rep_len(main, length(outs))[i]
    mc <- if (is.null(names(max_count))) max_count else
      max_count[[outs[i]]]
    .fcr_rootogram_one(ys[[i]], pmf, w, mu[, j], mc, style, overflow, plot,
                       title, xlab, ylab, ...)
  })
  names(res) <- outs
  invisible(if (length(outs) == 1L) res[[1L]] else res)
}

# Observed and expected frequencies of one outcome, and its plot. `pmf(k)`
# gives the probabilities of the counts 0..k (one row per observation).
.fcr_rootogram_one <- function(y, pmf, w, mu_hat, max_count, style, overflow,
                               plot, main, xlab, ylab, ...) {
  if (is.null(max_count)) max_count <- max(y, na.rm = TRUE)
  if (!is.numeric(max_count) || length(max_count) != 1L ||
      !is.finite(max_count) || max_count < 0 ||
      max_count != round(max_count)) {
    stop("`max_count` must be a single whole number of at least 0.",
         call. = FALSE)
  }
  P <- pmf(max_count)
  # The analysis rows: an observed outcome and finite fitted probabilities.
  # The same rows enter the observed and the expected frequencies, whatever
  # the display cutoff: a row whose outcome is above max_count still adds
  # its fitted probabilities to every displayed count.
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
  out
}
