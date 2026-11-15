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
#' For \code{countreg()} models of families that nest the Poisson model, a
#' likelihood-ratio test against a Poisson model with the same mean
#' specification, rows, offset and weights is reported. For one-parameter
#' mixtures where the Poisson model lies on the boundary of the parameter
#' space (NB1, NB2, PLN, PIG1, PIG2, PIG, PL), the p-value is halved (a
#' 50:50 mixture of chi-square distributions with 0 and 1 degrees of freedom).
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
  se <- if (is.matrix(V) && all(dim(V) == length(cf))) {
    .fcr_se(V)
  } else {
    rep(NA_real_, length(cf))
  }
  se[!is.finite(se)] <- NA_real_
  z <- cf / se
  ci <- stats::confint(object, level = level)
  tab <- cbind(Estimate = cf, `Std. Error` = se, `z value` = z,
               `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)), ci)
  rownames(tab) <- names(cf)

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
    logLik = as.numeric(ll),
    AIC = stats::AIC(object),
    BIC = stats::BIC(object),
    nobs = stats::nobs(object),
    n_groups = m[["n_groups"]],
    n_dropped = m[["n_dropped"]],
    converged = m[["converged"]],
    offset = m[["offset"]],
    lr_test = .fcr_lr_vs_poisson(object),
    digits = if (is.numeric(args$digits)) args$digits else 3
  ), class = "summary.flexCountReg")
}

#' @rdname summary.flexCountReg
#' @export
print.summary.flexCountReg <- function(x, digits = x$digits, ...) {
  cat("Call:\n")
  if (!is.null(x$call)) print(x$call) else print(x$formula)
  cat("\nModel:", x$model_type)
  if (!is.null(x$family)) cat(" | family:", x$family)
  if (!is.null(x$method)) cat(" | optimizer:", x$method)
  cat("\n")
  if (!is.null(x$iterations)) cat("Iterations:", x$iterations, "\n")
  if (!is.null(x$message)) cat("Convergence:", trimws(x$message), "\n")
  if (isFALSE(x$converged)) {
    cat("Warning: the optimizer did not report convergence.\n")
  }

  se_label <- switch(x$vcov_type %||% "hessian",
                     hessian = "Hessian-based",
                     robust = "robust (sandwich)",
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

  cat("\nLog-likelihood:", format(round(x$logLik, digits)),
      "| AIC:", format(round(x$AIC, digits)),
      "| BIC:", format(round(x$BIC, digits)), "\n")
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
# the same mean specification, rows, offset and weights. NULL when it does
# not apply (Poisson family, underreporting, or no training data).
.fcr_lr_vs_poisson <- function(object) {
  m <- object$model
  if (!inherits(object, "countreg") || identical(m$family, "POISSON") ||
      !is.null(m$underreport_formula) || is.null(object$data)) {
    return(NULL)
  }
  fam <- m$family
  nests <- c("NB1", "NB2", "NBP", "PLN", "PIG1", "PIG2", "PIG", "PL", "GP1",
             "GP2", "COM", "PLG", "PLL", "SI", "GW", "PW", "PGE")
  if (!fam %in% nests) return(NULL)
  d <- .fcr_training_rows(object)
  fd <- tryCatch(forge_new_data(m$blueprint %||% .fcr_blueprint(object), d),
                 error = function(e) NULL)
  if (is.null(fd)) return(NULL)
  y <- as.numeric(m[["y"]])
  w <- m[["weights"]]
  if (length(w) != length(y)) w <- rep(1, length(y))
  pf <- tryCatch(stats::glm.fit(fd$X, y, weights = w, offset = fd$offset,
                                family = stats::poisson()),
                 error = function(e) NULL)
  if (is.null(pf)) return(NULL)
  mu0 <- pf$fitted.values
  ll0 <- sum(w * stats::dpois(y, mu0, log = TRUE))
  df <- length(stats::coef(object)) - ncol(fd$X)
  stat <- max(0, 2 * (m$maximum - ll0))
  p <- stats::pchisq(stat, df = df, lower.tail = FALSE)
  boundary <- fam %in% c("NB1", "NB2", "PLN", "PIG1", "PIG2", "PIG", "PL") &&
    df == 1
  if (boundary) p <- p / 2
  list(statistic = stat, df = df, p.value = p, logLik_poisson = ll0,
       boundary = boundary)
}
