#' Tidy, glance and augment flexCountReg models
#'
#' Methods for the \code{tidy()}, \code{glance()} and \code{augment()}
#' generics from the \pkg{generics} package (also used by \pkg{broom}), so
#' flexCountReg models work with \pkg{broom}-style workflows, tables such as
#' \pkg{modelsummary}, and tidymodels.
#'
#' @param x a model estimated with this package.
#' @param conf.int logical; include confidence intervals?
#' @param conf.level the confidence level for \code{conf.int = TRUE}.
#' @param exponentiate logical; exponentiate the mean-model coefficients
#'   (giving incidence rate ratios)? Other components are not transformed.
#' @param data the data to augment. Defaults to the rows of the training data
#'   used in estimation.
#' @param new_data,newdata new data to augment instead of \code{data}; the
#'   outcome is optional.
#' @param type.residuals the residual type for \code{.resid}: \code{"response"}
#'   or \code{"pearson"} (see \code{\link{residuals.flexCountReg}}).
#' @param ... further arguments passed to \code{predict()} (for example
#'   \code{method} for \code{countreg.rp()} models).
#'
#' @details
#' \code{tidy()} has one row per parameter, with a \code{component} column
#' because these models have several equations:
#' \describe{
#'   \item{mean}{coefficients of the mean (count) model;}
#'   \item{distribution}{parameters of the count distribution (for example
#'     \code{ln(alpha)}), including their own regression coefficients;}
#'   \item{underreporting}{coefficients of the reporting model;}
#'   \item{zero}{coefficients of the zero model of zero-inflated and hurdle
#'     models;}
#'   \item{random_mean, random_scale, random_cholesky}{means, scales and
#'     Cholesky elements of random parameters;}
#'   \item{heterogeneity}{heterogeneity-in-means and -variances coefficients.}
#' }
#' Standard errors are those chosen when the model was fit (see
#' \code{vcov()}), and \code{p.value} uses the normal distribution.
#'
#' @returns \code{tidy()}: a tibble with columns \code{term},
#'   \code{component}, \code{estimate}, \code{std.error}, \code{statistic},
#'   \code{p.value} (and \code{conf.low}, \code{conf.high}).
#'   \code{glance()}: a one-row tibble of model statistics.
#'   \code{augment()}: the data as a tibble with \code{.fitted} (and
#'   \code{.resid} when the outcome is available).
#'
#' @examples
#' \donttest{
#' data("washington_roads")
#' m <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
#'               family = "NB2", method = "BFGS")
#' tidy(m, conf.int = TRUE)
#' tidy(m, exponentiate = TRUE)
#' glance(m)
#' head(augment(m))
#' }
#' @name flexCountReg-tidiers
NULL

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics augment
#' @export
generics::augment

# Model component of each parameter name.
.component_of <- function(nm) {
  out <- rep("mean", length(nm))
  dist <- grepl("^ln\\(|^gamma($|:)|^p($|:)|^lnalpha$", nm)
  out[dist] <- "distribution"
  out[grepl("^Underreporting:", nm)] <- "underreporting"
  out[grepl("^Zero:", nm)] <- "zero"
  out[grepl(":Mean$", nm)] <- "random_mean"
  out[grepl(":(St\\.Dev|Spread)$", nm)] <- "random_scale"
  out[grepl("^Chol:", nm)] <- "random_cholesky"
  out[grepl("^Het(Mean|Var):", nm)] <- "heterogeneity"
  out
}

#' @rdname flexCountReg-tidiers
#' @export
tidy.flexCountReg <- function(x, conf.int = FALSE, conf.level = 0.95,
                              exponentiate = FALSE, ...) {
  cf <- stats::coef(x)
  V <- stats::vcov(x)
  se <- if (is.matrix(V) && all(dim(V) == length(cf))) {
    .fcr_se(V)
  } else {
    rep(NA_real_, length(cf))
  }
  out <- tibble::tibble(term = names(cf), component = .component_of(names(cf)),
                        estimate = unname(cf), std.error = unname(se))
  out$statistic <- out$estimate / out$std.error
  out$p.value <- 2 * stats::pnorm(-abs(out$statistic))
  if (conf.int) {
    z <- stats::qnorm((1 + conf.level) / 2)
    out$conf.low <- out$estimate - z * out$std.error
    out$conf.high <- out$estimate + z * out$std.error
  }
  if (exponentiate) {
    i <- out$component == "mean"
    for (col in intersect(c("estimate", "conf.low", "conf.high"),
                          names(out))) {
      out[[col]][i] <- exp(out[[col]][i])
    }
    # std.error on the rate-ratio scale by the delta method
    out$std.error[i] <- out$estimate[i] * out$std.error[i]
  }
  out
}

#' @rdname flexCountReg-tidiers
#' @export
glance.flexCountReg <- function(x, ...) {
  m <- x$model
  tibble::tibble(
    model = class(x)[1],
    family = m$family %||% NA_character_,
    logLik = as.numeric(stats::logLik(x)),
    AIC = stats::AIC(x),
    BIC = stats::BIC(x),
    df = length(stats::coef(x)),
    nobs = stats::nobs(x),
    n_groups = as.integer(m[["n_groups"]] %||% NA_integer_),
    converged = as.logical(m[["converged"]] %||% NA),
    vcov_type = m[["vcov_type"]] %||% NA_character_)
}

#' @rdname flexCountReg-tidiers
#' @export
augment.flexCountReg <- function(x, data = NULL, new_data = NULL,
                                 newdata = NULL,
                                 type.residuals = c("response", "pearson"),
                                 ...) {
  type.residuals <- match.arg(type.residuals)
  if (is.null(new_data)) new_data <- newdata
  if (!is.null(new_data)) {
    nd <- as.data.frame(new_data)
    out <- tibble::as_tibble(nd)
    out$.fitted <- as.vector(stats::predict(x, newdata = nd, ...))
    y_name <- all.vars(x$formula)[1]
    if (y_name %in% names(nd)) {
      if (type.residuals == "pearson") {
        stop("Pearson residuals are only available for the training data.",
             call. = FALSE)
      }
      out$.resid <- as.numeric(nd[[y_name]]) - out$.fitted
    }
    return(out)
  }
  d <- if (is.null(data)) .fcr_training_rows(x) else {
    dd <- as.data.frame(data)
    rows <- x$model[["rows_used"]]
    if (!is.null(rows) && !is.null(x$data) &&
        nrow(dd) == nrow(as.data.frame(x$data))) dd[rows, , drop = FALSE] else dd
  }
  out <- tibble::as_tibble(d)
  out$.fitted <- as.vector(stats::fitted(x, ...))
  out$.resid <- as.vector(stats::residuals(x, type = type.residuals, ...))
  out
}
