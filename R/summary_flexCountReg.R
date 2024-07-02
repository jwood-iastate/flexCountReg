#' Custom summary method for flexCountReg models
#' 
#' @param object A flexCountReg model object.
#' @param confint_level A numeric value between 0 and 1 indicating the confidence level for confidence intervals. Default is 0.95.
#' @param digits Number of digits (decimal places) to round to. Default is 3.
#' 
#' @import tibble dplyr
#' @details
#' This summary method accounts for bootstrapped standard errors (when used). Bootstrapped standard errors are currently only implemented in \code{\link{pwiebreg}}.
#' 
#' @examples
#' # Poisson-Weibull
#' pw_rp <- pwiebreg(Total_crashes ~ lnlength + lnaadt,
#'                                  data = washington_roads,
#'                                  ndraws = 10,
#'                                  bootstraps = 10)
#' summary(pw_rp)
#' 
#' @export
summary.flexCountReg <- function(object, confint_level = 0.95, digits = 3) {
  object <- object$model
  cat("Call:\n", deparse(object$formula), "\n")
  cat("\n", "Method: ", object$type, "\n")
  cat("Iterations: ", object$iterations, "\n")
  cat("Convergence: ", object$message, "\n")
  cat("Log-likelihood: ", object$maximum, "\n")
  
  cat("\nParameter Estimates:\n")
  
  if (!is.null(object$bootstraps)) {
    mod.sum <- tibble::tibble(parameter = names(object$estimate), 
                              coeff = object$estimate, 
                              `Std. Err.` = object$bootstrapped_se)
  } else {
    se <- sqrt(diag(-1/(object$hessian)))
    mod.sum <- tibble::tibble(parameter = names(object$estimate), 
                              coeff = object$estimate, 
                              `Std. Err.` = se)
  }
  
  mod.sum$`t-stat` <- mod.sum$coeff / mod.sum$`Std. Err.`
  mod.sum$`p-value` <- 2 * pnorm(-1 * abs(mod.sum$`t-stat`), mean = 0, sd = 1)
  
  lower.level <- (1 - confint_level) / 2
  upper.level <- 1 - lower.level
  
  ci_z.lower <- qnorm(lower.level)
  ci_z.upper <- qnorm(upper.level)
  
  mod.sum$`lower CI` <- mod.sum$coeff + ci_z.lower * mod.sum$`Std. Err.`
  mod.sum$`upper CI` <- mod.sum$coeff + ci_z.upper * mod.sum$`Std. Err.`
  
  mod.sum <- mod.sum %>%
    dplyr::mutate(across(where(is.numeric), round, digits = digits))
  
  print(mod.sum)
}