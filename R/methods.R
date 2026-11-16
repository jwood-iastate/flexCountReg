#' Standard methods for flexCountReg models
#'
#' Methods that make flexCountReg models work like other fitted models in R,
#' including with \code{AIC()}, \code{BIC()}, \code{confint()},
#' \code{update()} and packages that build on these generics.
#'
#' @param object,x a model estimated with this package.
#' @param parm a specification of which parameters are to be given confidence
#'   intervals, either a vector of numbers or a vector of names. If missing,
#'   all parameters are considered.
#' @param level the confidence level required.
#' @param type the type of residuals: \code{"response"} (observed minus
#'   fitted), \code{"pearson"} (divided by the standard deviation implied by
#'   the fitted distribution) or \code{"quantile"} (randomized quantile
#'   residuals of Dunn and Smyth, 1996, which are standard normal when the
#'   model is correct). Pearson and quantile residuals are available for
#'   \code{countreg()} models.
#' @param nsim number of response vectors to simulate.
#' @param seed an optional random seed, as in \code{\link[stats]{simulate}}.
#' @param digits the number of significant digits to print.
#' @param ... further arguments passed to or from other methods. For
#'   \code{fitted()} and \code{residuals()} of \code{countreg.rp()} models,
#'   this can include \code{method} (see \code{\link{predict.flexCountReg}}).
#'
#' @details
#' \code{coef()} returns all estimated parameters, including distribution,
#' underreporting, random-parameter and heterogeneity parameters, named as in
#' \code{summary()}. \code{vcov()} returns the covariance matrix chosen at
#' fit time: Hessian-based, robust (sandwich) or bootstrapped (see the
#' \code{stderr} argument of each estimator).
#'
#' \code{nobs()} is the number of rows used in estimation. For panel models
#' (\code{renb()}, \code{poisLind.re()}, \code{repois()}, \code{negmultinom()}
#' and \code{countreg.rp()} with \code{panel_id}), the \code{nobs} attribute
#' of \code{logLik()}, which \code{BIC()} uses, is the number of panels,
#' because the panels are the independent units.
#'
#' \code{fitted()} and \code{residuals()} are for the rows used in
#' estimation (rows with missing values are not included). \code{simulate()}
#' draws new responses from the fitted distribution for the same rows; it is
#' available for \code{countreg()} models.
#'
#' @returns \code{coef()}: named numeric vector. \code{vcov()}: matrix.
#'   \code{logLik()}: an object of class \code{logLik}. \code{nobs()}:
#'   integer. \code{fitted()}, \code{residuals()}: numeric vectors.
#'   \code{confint()}: matrix with lower and upper limits. \code{formula()}:
#'   the main model formula. \code{terms()}: its terms object.
#'   \code{simulate()}: a data frame with \code{nsim} columns.
#'
#' @examples
#' \donttest{
#' data("washington_roads")
#' m <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
#'               family = "NB2", method = "BFGS")
#' coef(m)
#' confint(m)
#' AIC(m)
#' BIC(m)
#' head(residuals(m, type = "pearson"))
#' sims <- simulate(m, nsim = 2, seed = 1)
#' }
#' @name flexCountReg-methods
NULL

#' @rdname flexCountReg-methods
#' @export
coef.flexCountReg <- function(object, ...) {
  object$model[["estimate"]]
}

#' @rdname flexCountReg-methods
#' @export
vcov.flexCountReg <- function(object, ...) {
  object$model[["vcov"]]
}

#' @rdname flexCountReg-methods
#' @importFrom stats nobs
#' @export
nobs.flexCountReg <- function(object, ...) {
  n <- object$model[["nobs"]]
  if (is.null(n) && !is.null(object$data)) n <- nrow(object$data)
  as.integer(n)
}

#' @rdname flexCountReg-methods
#' @export
logLik.flexCountReg <- function(object, ...) {
  n_indep <- object$model[["n_groups"]]
  if (is.null(n_indep)) n_indep <- stats::nobs(object)
  structure(object$model[["maximum"]], df = length(stats::coef(object)),
            nobs = as.integer(n_indep), class = "logLik")
}

#' @rdname flexCountReg-methods
#' @export
formula.flexCountReg <- function(x, ...) {
  x$formula
}

#' @rdname flexCountReg-methods
#' @export
terms.flexCountReg <- function(x, ...) {
  if (!is.null(x$model[["blueprint"]])) return(x$model[["blueprint"]]$terms)
  stats::terms(x$formula)
}

#' @rdname flexCountReg-methods
#' @export
confint.flexCountReg <- function(object, parm, level = 0.95, ...) {
  cf <- stats::coef(object)
  se <- .fcr_se(stats::vcov(object), length(cf))
  z <- stats::qnorm((1 + level) / 2)
  ci <- cbind(cf - z * se, cf + z * se)
  pct <- paste(format(100 * c((1 - level) / 2, (1 + level) / 2),
                      trim = TRUE, scientific = FALSE, digits = 3), "%")
  dimnames(ci) <- list(names(cf), pct)
  if (missing(parm)) ci else ci[parm, , drop = FALSE]
}

# Training rows used in estimation.
.fcr_training_rows <- function(object) {
  if (is.null(object$data)) {
    stop("The model was fit with keep_data = FALSE, so the training data ",
         "are not available. Pass `newdata` to predict() instead.",
         call. = FALSE)
  }
  d <- as.data.frame(object$data)
  rows <- object$model[["rows_used"]]
  if (is.null(rows)) rows <- seq_len(nrow(d))
  d[rows, , drop = FALSE]
}

#' @rdname flexCountReg-methods
#' @export
fitted.flexCountReg <- function(object, ...) {
  fv <- object$model[["fitted.values"]]
  if (!is.null(fv) && length(list(...)) == 0) return(as.vector(fv))
  as.vector(stats::predict(object, newdata = .fcr_training_rows(object), ...))
}

#' @rdname flexCountReg-methods
#' @export
residuals.flexCountReg <- function(object,
                                   type = c("response", "pearson", "quantile"),
                                   ...) {
  type <- match.arg(type)
  if (type == "quantile") return(.fcr_quantile_residuals(object))
  y <- object$model[["y"]]
  if (is.null(y)) {
    d <- .fcr_training_rows(object)
    y <- stats::model.response(stats::model.frame(object$formula, d))
  }
  y <- as.numeric(y)
  mu <- stats::fitted(object, ...)
  r <- y - mu
  if (type == "pearson") r <- r / sqrt(.fcr_fitted_variance(object))
  r
}

# Variance of the observed outcome for each training row (countreg only).
# The observed count follows the model's family with mean parameter
# exp(X beta + offset) * p, where p is the reporting probability (p = 1
# without underreporting), as in the likelihood of countreg().
.fcr_fitted_variance <- function(object) {
  if (!inherits(object, "countreg")) {
    stop("Pearson residuals are available for countreg() models.",
         call. = FALSE)
  }
  parts <- .fcr_countreg_parts(object, .fcr_training_rows(object))
  v <- .fcr_family_variance(parts$family, parts$mu_link * parts$p_report,
                            parts$alpha, parts$sigma)
  if (is.null(v)) {
    stop("A variance function is not available for family \"",
         parts$family, "\".", call. = FALSE)
  }
  if (!is.null(parts$zero_prob)) {
    # Mixture moments: E[Y^2] = w (V + m^2) for the count part with weight w
    m <- parts$mu_true
    zp <- parts$zero_prob
    w <- if (parts$zero_type == "inflated") 1 - zp else zp / (1 - parts$f0)
    v <- w * (v + m^2) - (w * m)^2
  }
  v
}

#' @rdname flexCountReg-methods
#' @importFrom stats simulate
#' @export
simulate.flexCountReg <- function(object, nsim = 1, seed = NULL, ...) {
  if (!inherits(object, "countreg")) {
    stop("simulate() is available for countreg() models.", call. = FALSE)
  }
  if (!is.null(seed)) set.seed(seed)
  rng_state <- if (exists(".Random.seed", envir = globalenv()))
    get(".Random.seed", envir = globalenv())
  parts <- .fcr_countreg_parts(object, .fcr_training_rows(object))
  n <- length(parts$mu_link)
  draw <- function(idx) {
    .fcr_family_r(parts$family, length(idx),
                  (parts$mu_link * parts$p_report)[idx],
                  if (!is.null(parts$alpha)) parts$alpha[idx],
                  if (!is.null(parts$sigma)) parts$sigma[idx])
  }
  sims <- lapply(seq_len(nsim), function(i) {
    y <- as.numeric(draw(seq_len(n)))
    zp <- parts$zero_prob
    if (!is.null(zp)) {
      if (parts$zero_type == "inflated") {
        y[stats::runif(n) < zp] <- 0
      } else {
        # Positive counts from the zero-truncated family (redraw zeros)
        need <- which(y == 0)
        for (iter in seq_len(1000)) {
          if (!length(need)) break
          y[need] <- draw(need)
          need <- need[y[need] == 0]
        }
        y[stats::runif(n) >= zp] <- 0
      }
    }
    y
  })
  out <- as.data.frame(sims)
  names(out) <- paste0("sim_", seq_len(nsim))
  attr(out, "seed") <- rng_state
  out
}

#' @rdname flexCountReg-methods
#' @export
print.flexCountReg <- function(x, digits = max(3L, getOption("digits") - 3L),
                               ...) {
  m <- x$model
  cat("flexCountReg model:", class(x)[1])
  if (!is.null(m$family)) cat(" | family:", m$family)
  cat("\n")
  if (!is.null(x$call)) {
    cat("Call: ")
    print(x$call)
  }
  cat("\nCoefficients:\n")
  print.default(format(stats::coef(x), digits = digits), print.gap = 2L,
                quote = FALSE)
  cat("\nLog-likelihood:", format(m[["maximum"]], digits = digits + 3),
      "| N =", stats::nobs(x))
  if (!is.null(m[["n_groups"]]) && m[["n_groups"]] != stats::nobs(x)) {
    cat(" | panels =", m[["n_groups"]])
  }
  cat("\n")
  if (isFALSE(m[["converged"]])) cat("Warning: the optimizer did not converge.\n")
  invisible(x)
}
