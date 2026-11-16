# Probability predictions, randomized quantile residuals and rootograms for
# countreg() models.

# P(Y = k | x_i) for k = 0..max_count and each row of `data`, evaluated with
# the fitted family exactly as in the likelihood (mean parameter
# exp(X beta + offset) * p, the same distribution parameters and the same
# Halton draws). Returns an n x (max_count + 1) matrix.
.fcr_prob_matrix <- function(object, data, max_count) {
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

# Randomized quantile residuals (Dunn and Smyth, 1996): for each training
# row, u ~ Uniform(F(y - 1), F(y)) under the fitted distribution, and the
# residual is qnorm(u).
.fcr_quantile_residuals <- function(object, chunk = 1000L) {
  d <- .fcr_training_rows(object)
  y <- as.numeric(object$model[["y"]])
  n <- length(y)
  lower <- upper <- numeric(n)
  # Rows sorted by y and processed in chunks, so memory is bounded and a
  # single large count does not widen the probability grid of every row
  ord <- order(y)
  for (idx in split(ord, ceiling(seq_len(n) / chunk))) {
    yc <- y[idx]
    P <- .fcr_prob_matrix(object, d[idx, , drop = FALSE],
                          max_count = max(yc))
    cdf <- P
    if (ncol(P) > 1) for (j in 2:ncol(P)) cdf[, j] <- cdf[, j - 1] + P[, j]
    r <- seq_along(idx)
    upper[idx] <- cdf[cbind(r, yc + 1)]
    lower[idx] <- ifelse(yc > 0, cdf[cbind(r, pmax(yc, 1))], 0)
  }
  upper <- pmin(upper, 1)
  lower <- pmin(lower, upper)
  stats::qnorm(stats::runif(n, lower, upper))
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
#'   count.
#' @param style \code{"hanging"} (bars hang from the expected curve) or
#'   \code{"standing"} (bars stand on zero).
#' @param plot logical; draw the plot?
#' @param main,xlab,ylab plot labels.
#' @param ... further arguments passed to \code{\link[graphics]{plot}}.
#'
#' @returns Invisibly, a data frame with columns \code{count},
#'   \code{observed} and \code{expected} (frequencies), of class
#'   \code{flexCountReg_rootogram}.
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
                      style = c("hanging", "standing"),
                      plot = TRUE, main = NULL, xlab = "Count",
                      ylab = "sqrt(Frequency)", ...) {
  style <- match.arg(style)
  d <- if (is.null(newdata)) .fcr_training_rows(object) else
    as.data.frame(newdata)
  y <- as.numeric(stats::model.response(stats::model.frame(
    object$formula, d, na.action = stats::na.pass)))
  if (is.null(max_count)) max_count <- max(y, na.rm = TRUE)
  P <- .fcr_prob_matrix(object, d, max_count)
  # Frequency weights of the training rows count each row w times
  w <- if (is.null(newdata) && length(object$model$weights) == nrow(d))
    object$model$weights else rep(1, nrow(d))
  keep <- !is.na(y) & stats::complete.cases(P) & y <= max_count
  counts <- 0:max_count
  obs <- as.vector(tapply(w[keep], factor(y[keep], levels = counts), sum))
  obs[is.na(obs)] <- 0
  expd <- colSums(P[keep, , drop = FALSE] * w[keep])
  out <- data.frame(count = counts, observed = obs, expected = unname(expd))
  class(out) <- c("flexCountReg_rootogram", "data.frame")

  if (plot) {
    so <- sqrt(obs)
    se <- sqrt(expd)
    if (style == "hanging") {
      top <- se; bottom <- se - so
    } else {
      top <- so; bottom <- rep(0, length(so))
    }
    if (is.null(main)) {
      main <- paste("Rootogram:", object$model$family %||% class(object)[1])
    }
    graphics::plot(range(counts) + c(-0.5, 0.5),
                   range(c(0, top, bottom, se)), type = "n", main = main,
                   xlab = xlab, ylab = ylab, ...)
    graphics::rect(counts - 0.4, bottom, counts + 0.4, top, col = "grey80",
                   border = "grey40")
    graphics::abline(h = 0, lty = 2)
    graphics::lines(counts, se, col = "firebrick", lwd = 2)
    graphics::points(counts, se, col = "firebrick", pch = 19)
  }
  invisible(out)
}
