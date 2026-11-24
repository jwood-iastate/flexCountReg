#' Summarize a flexCountReg model
#'
#' @param object A flexCountReg model object.
#' @param level the confidence level for the confidence intervals. Default is
#'   0.95.
#' @param x a \code{summary.flexCountReg} object.
#' @param digits the number of digits (decimal places) to round to when
#'   printing. Default is 3.
#' @param ... Optional arguments. The arguments \code{confint_level} and
#'   \code{digits} of earlier versions are still accepted.
#'
#' @details
#' The standard errors are those chosen when the model was fit: Hessian-based,
#' robust (sandwich), or bootstrapped. The z statistics and p-values use the
#' normal distribution.
#'
#' For \code{countreg()} models of families that nest the Poisson model with
#' a regular restriction, a likelihood-ratio test against a Poisson model
#' with the same mean specification, rows, offset and weights is reported.
#' For NB1, NB2, PLN, PIG1, PIG2 and PIG with a scalar dispersion parameter
#' the Poisson model is on the boundary of the parameter space, and the
#' p-value uses a 50:50 mixture of chi-square distributions with 0 and 1
#' degrees of freedom; for GP1, GP2, COM, DPO, HYPERPOIS and GAMMACOUNT the
#' restriction is interior. No test is reported for PL, PLG and PLL (which do
#' not reduce to Poisson), for NBP, SI, PW, PGE, GW, DELAPORTE, HERMITE*, PA
#' and NBGE (Poisson is only a limit at the edge of the parameter space, or
#' not a special case), for DW, for boundary families with a dispersion
#' formula, or for models with zero or underreporting components.
#'
#' @returns An object of class \code{summary.flexCountReg}, a list with
#'   \code{coefficients} (a matrix with the estimate, standard error, z value,
#'   p-value and confidence limits of every parameter), \code{vcov_type},
#'   \code{logLik}, \code{AIC}, \code{BIC}, \code{nobs}, \code{n_groups},
#'   \code{n_dropped}, \code{converged}, \code{offset} and \code{lr_test}.
#'   Printing it shows the results. Use \code{tibble::as_tibble()} on it for
#'   the coefficient table as a tibble (the value returned by
#'   \code{summary()} in earlier versions), or \code{\link{tidy.flexCountReg}}.
#'
#' @examples
#' \donttest{
#' data("washington_roads")
#' nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50,
#'                 data = washington_roads, family = "NB2",
#'                 dis_param_formula_1 = ~ speed50, method = "BFGS")
#' s <- summary(nb2)
#' s
#' s$coefficients
#' tibble::as_tibble(s)
#' }
#'
#' @export
summary.flexCountReg <- function(object, level = 0.95, ...) {
  args <- list(...)
  if (is.numeric(args$confint_level)) level <- args$confint_level
  m <- object$model

  cf <- stats::coef(object)
  V <- stats::vcov(object)
  se <- .fcr_se(V, length(cf))
  se[!is.finite(se)] <- NA_real_
  z <- cf / se
  ci <- stats::confint(object, level = level)
  tab <- cbind(Estimate = cf, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)), ci)
  rownames(tab) <- names(cf)
  # Smooth terms are summarized by term, not by basis coefficient
  if (length(m[["smooth_coefs"]])) {
    tab <- tab[-m[["smooth_coefs"]], , drop = FALSE]
  }

  ll <- stats::logLik(object)
  structure(list(
    call = object$call,
    formula = object$formula,
    model_type = class(object)[1],
    family = m$family,
    method = m[["method"]],
    iterations = m[["iterations"]],
    message = m[["message"]],
    coefficients = tab,
    level = level,
    vcov_type = m[["vcov_type"]],
    n_boot_ok = m[["n_boot_ok"]],
    bootstraps = m[["bootstraps"]],
    n_clusters = m[["n_clusters"]],
    smooth_table = m[["smooth_table"]],
    edf_total = m[["df_total"]],
    fe_label = m[["fe_method_label"]],
    fe_dropped = if (!is.null(m[["fe_method_label"]]))
      sum(lengths(m[["dropped"]] %||% list())),
    logLik = as.numeric(ll),
    AIC = stats::AIC(object),
    BIC = stats::BIC(object),
    nobs = stats::nobs(object),
    n_groups = m[["n_groups"]],
    n_dropped = m[["n_dropped"]],
    converged = m[["converged"]],
    offset = m[["offset"]],
    weight_type = m[["weight_type"]],
    # Not for penalized fits: their smoothing parameters are selected
    # separately, so the likelihood-ratio statistic has no chi-squared basis
    lr_test = if (is.null(m[["smooth_table"]])) .fcr_lr_vs_poisson(object),
    digits = if (is.numeric(args$digits)) args$digits else 3
  ), class = "summary.flexCountReg")
}

#' @rdname summary.flexCountReg
#' @export
print.summary.flexCountReg <- function(x, digits = x$digits, ...) {
  cat("Call:\n")
  if (!is.null(x$call)) print(x$call) else print(x$formula)
  cat("\nModel:", x$model_type)
  if (!is.null(x$family)) cat(" | family:", .fcr_family_label(x$family))
  if (!is.null(x$method)) cat(" | optimizer:", x$method)
  cat("\n")
  if (!is.null(x$fe_label)) {
    cat(x$fe_label, "\n")
    if (isTRUE(x$fe_dropped > 0)) {
      cat("Groups dropped (all-zero counts or one observation):",
          x$fe_dropped, "\n")
    }
  }
  if (!is.null(x$iterations)) cat("Iterations:", x$iterations, "\n")
  if (!is.null(x$message)) cat("Convergence:", trimws(x$message), "\n")
  if (isFALSE(x$converged)) {
    cat("Warning: the optimizer did not report convergence.\n")
  }

  se_label <- switch(x$vcov_type %||% "hessian",
                     hessian = "Hessian-based",
                     robust = "robust (sandwich)",
                     bayesian = "Bayesian (penalized Hessian)",
                     cluster = paste0("cluster-robust (",
                                      paste(x$n_clusters, "clusters of",
                                            names(x$n_clusters),
                                            collapse = "; "), ")"),
                     bootstrap = paste0("bootstrapped (", x$n_boot_ok,
                                        " of ", x$bootstraps,
                                        " replicates)"),
                     x$vcov_type)
  cat("\nCoefficients (standard errors: ", se_label, "):\n", sep = "")
  stats::printCoefmat(x$coefficients[, 1:4, drop = FALSE], digits = digits,
                      P.values = TRUE, has.Pvalue = TRUE, na.print = "NA")
  if (length(x$offset)) {
    cat("Offset:", paste(x$offset, collapse = " + "),
        "(coefficient fixed at 1)\n")
  }
  if (!is.null(x$smooth_table)) {
    st <- x$smooth_table
    cat("\nSmooth terms (approximate tests):\n")
    tb <- data.frame(edf = round(st$edf, 2), k = st$k,
                     lambda = st$lambda,
                     `Chi.sq` = round(st$chisq, 2),
                     `p-value` = format.pval(st$p.value, digits = 3),
                     check.names = FALSE)
    rownames(tb) <- ifelse(st$parameter == "mean", st$term,
                           paste0(st$parameter, ": ", st$term))
    print(tb)
    cat("Total effective degrees of freedom:",
        format(round(x$edf_total, 2)), "\n")
  }

  cat("\nLog-likelihood:", format(round(x$logLik, digits)),
      "| AIC:", format(round(x$AIC, digits)),
      "| BIC:", format(round(x$BIC, digits)), "\n")
  if (identical(x$weight_type, "sampling")) {
    cat("Note: with sampling weights the log-likelihood is a",
        "pseudo-likelihood; LR tests, AIC and BIC are not valid.\n")
  }
  cat("N:", x$nobs)
  if (!is.null(x$n_groups) && x$n_groups != x$nobs) {
    cat(" | panels:", x$n_groups)
  }
  if (!is.null(x$n_dropped) && x$n_dropped > 0) {
    cat(" | rows removed (missing values):", x$n_dropped)
  }
  cat("\n")
  if (!is.null(x$lr_test)) {
    lr <- x$lr_test
    cat("LR test vs Poisson: chi-sq =", format(round(lr$statistic, digits)),
        "on", lr$df, "df, p-value =", format.pval(lr$p.value, digits = 3),
        if (lr$boundary) "(boundary-corrected)", "\n")
  }
  invisible(x)
}

#' @importFrom tibble as_tibble
#' @export
tibble::as_tibble

#' @rdname summary.flexCountReg
#' @export
as_tibble.summary.flexCountReg <- function(x, ...) {
  tab <- x$coefficients
  out <- tibble::tibble(parameter = rownames(tab),
                        coeff = unname(tab[, "Estimate"]),
                        `Std. Err.` = unname(tab[, "Std. Error"]),
                        `z value` = unname(tab[, "z value"]),
                        `p-value` = unname(tab[, "Pr(>|z|)"]),
                        `lower CI` = unname(tab[, 5]),
                        `upper CI` = unname(tab[, 6]))
  num <- vapply(out, is.numeric, logical(1))
  out[num] <- lapply(out[num], round, digits = x$digits)
  out
}

# Likelihood-ratio test of a countreg() model against a Poisson model with
# the same mean specification, rows, offset and weights. NULL when Poisson is
# not nested in the model with a regular restriction (see .fcr_nesting()),
# or when the training data are not available.
.fcr_lr_vs_poisson <- function(object) {
  m <- object$model
  if (!inherits(object, "countreg") || is.null(object$data)) return(NULL)
  if (identical(m$weight_type, "sampling")) return(NULL)
  full <- .fcr_model_spec(object)
  bp <- m$blueprint %||% .fcr_blueprint(object)
  nest <- .fcr_nesting(.fcr_base_spec("POISSON", bp$x_names), full)
  if (!nest$nested || !isTRUE(m$converged)) return(NULL)
  d <- .fcr_training_rows(object)
  fd <- tryCatch(forge_new_data(bp, d), error = function(e) NULL)
  if (is.null(fd)) return(NULL)
  y <- as.numeric(m[["y"]])
  w <- m[["weights"]]
  if (length(w) != length(y)) w <- rep(1, length(y))
  pf <- tryCatch(stats::glm.fit(fd$X, y, weights = w, offset = fd$offset,
                                family = stats::poisson()),
                 error = function(e) NULL)
  if (is.null(pf)) return(NULL)
  ll0 <- sum(w * stats::dpois(y, pf$fitted.values, log = TRUE))
  if (!isTRUE(pf$converged)) return(NULL)
  test <- .fcr_lr_pvalue(m$maximum, ll0, nest)
  if (is.na(test$LR)) return(NULL)
  list(statistic = test$LR, df = nest$df, p.value = test$p,
       logLik_poisson = ll0, boundary = nest$type == "boundary")
}
