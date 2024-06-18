#' Summary method for random parameter negative binomial model
#' @param object a rpnb model, as produced by the rpnb function.
#' @param confint_level A numeric value between 0 and 1 indicating the confidence level for confidence intervals. Default is 0.95.
#' @param digits Number of digits to round to.
#' @param ...
#' @export
summary.rpnb <- function(object, confint_level = 0.95, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nSummary of Random Parameter Negative Binomial Model\n")
  cat("Call:\n", deparse(object$formula), "\n")
  cat("\n", object$form, "Estimation\n")
  cat("Method: ", object$method, "\n")
  cat("Maximum iterations: ", object$control$iterlim, "\n")
  cat("Convergence: ", ifelse(object$convergence == 0, "Converged", "Did not converge"), "\n")
  cat("Log-likelihood: ", object$logLik, "\n")
  
  if (object$correlated) {
    cat("Random Parameter Correlation:\n")
    print(object$Correlation)
  }
  
  cat("\nFixed Effects:\n")
  print(coef(object))
  
  cat("\nRandom Parameters:\n")
  print(coef(object)$rpr)
  
  cat("\nStandard Errors:\n")
  if (!is.null(object$bootstrapped_se)) {
    print(object$bootstrapped_se)
  } else {
    print(object$sd)
  }
  
  cat("\nDispersion Parameter (alpha):\n")
  cat("Estimate: ", object$alpha, "\n")
  if (!is.null(object$bootstrapped_se)) {
    cat("Bootstrapped SE: ", sd(object$bootstrapped_se[, "ln(alpha)"]), "\n")
    cat("t-value: ", object$alpha / sd(object$bootstrapped_se[, "ln(alpha)"]), "\n")
  } else {
    cat("Standard Error: ", object$sd["ln(alpha)"], "\n")
    cat("t-value: ", object$alpha / object$sd["ln(alpha)"], "\n")
  }
  
  if (object$form == 'nbp') {
    cat("\nP Parameter:\n")
    cat("Estimate: ", object$P, "\n")
    if (!is.null(object$bootstrapped_se)) {
      cat("Bootstrapped SE: ", sd(object$bootstrapped_se[, "P"]), "\n")
      cat("t-value: ", object$P / sd(object$bootstrapped_se[, "P"]), "\n")
    } else {
      cat("Standard Error: ", object$sd["P"], "\n")
      cat("t-value: ", object$P / object$sd["P"], "\n")
    }
  }
  
  cat("\nConfidence Intervals:\n")
  if (!is.null(object$bootstrapped_se)) {
    ci <- boot::boot.ci(object = object$bootstrapped_se, type = "norm", conf = confint_level)
    ci_ln_alpha <- ci$normal[, "CI"]
    cat("ln(alpha) (", confint_level * 100, "%): [", ci_ln_alpha[1], ", ", ci_ln_alpha[2], "]\n")
    if (object$form == 'nbp') {
      ci_P <- ci$normal[, "CI"]
      cat("P (", confint_level * 100, "%): [", ci_P[1], ", ", ci_P[2], "]\n")
    }
  } else {
    ci_ln_alpha <- confint(object, "ln(alpha)", level = confint_level)
    cat("ln(alpha) (", confint_level * 100, "%): [", ci_ln_alpha[1], ", ", ci_ln_alpha[2], "]\n")
    if (object$form == 'nbp') {
      ci_P <- confint(object, "P", level = confint_level)
      cat("P (", confint_level * 100, "%): [", ci_P[1], ", ", ci_P[2], "]\n")
    }
  }
  
  invisible(NULL)
}


#' Custom summary method for maxLik objects
#' @param object An object of class maxLik.
#' @param confint_level A numeric value between 0 and 1 indicating the confidence level for confidence intervals. Default is 0.95.
#' @param digits Number of digits to round to.
#' @param ...
#' @export
summary.maxLik <- function(object, confint_level = 0.95, digits = max(3, getOption("digits") - 3), ...) {
  if (inherits(object, "rpnb")) {
    summary.rpnb(object, confint_level = confint_level, digits = digits, ...)
  } else {
    stats:::summary.glm(object, confint_level = confint_level, digits = digits, ...)
  }
}