#' Custom summary method for flexCountReg models
#' 
#' @param object A flexCountReg model object.
#' @param ... Optional parameters that include `confint_level` and `digits`.
#' 
#' @note Optional parameter `confint_level`: A numeric value between 0 and 1 indicating the confidence level for confidence intervals. Default is 0.95.
#' @note Optional parameter `digits`: Number of digits (decimal places) to round to. Default is 3.
#' 
#' @import tibble 
#' @importFrom dplyr mutate across  %>% where
#' @details
#' This summary method accounts for bootstrapped or robust standard errors (when used).
#' 
#' @examples
#' # NB2 Model
#' data("washington_roads")
#' washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
#' nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'                data = washington_roads, family = "NB2",
#'                dis_param_formula_1 = ~ speed50, method='BFGS')
#' summary(nb2)
#' 
#' @export
summary.flexCountReg <- function(object, ...) {
  # Extract optional parameters from '...'
  args <- list(...)
  if (is.numeric(args$confint_level)) confint_level <- args$confint_level else confint_level <- 0.95  
  if (is.numeric(args$digits)) digits <- args$digits else digits <- 3  

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
  } else if (!is.null(object$se)) {
    mod.sum <- tibble::tibble(parameter = names(object$estimate), 
                              coeff = object$estimate, 
                              `Std. Err.` = object$se)
    
  } # else {
  #   se <- sqrt(diag(-1/(object$hessian)))
  #   mod.sum <- tibble::tibble(parameter = names(object$estimate), 
  #                             coeff = object$estimate, 
  #                             `Std. Err.` = se)
  # }
  # 
  mod.sum$`t-stat` <- mod.sum$coeff / mod.sum$`Std. Err.`
  mod.sum$`p-value` <- 2 * pnorm(-1 * abs(mod.sum$`t-stat`), mean = 0, sd = 1)
  
  lower.level <- (1 - confint_level) / 2
  upper.level <- 1 - lower.level
  
  ci_z.lower <- qnorm(lower.level)
  ci_z.upper <- qnorm(upper.level)
  
  mod.sum$`lower CI` <- mod.sum$coeff + ci_z.lower * mod.sum$`Std. Err.`
  mod.sum$`upper CI` <- mod.sum$coeff + ci_z.upper * mod.sum$`Std. Err.`
  
  mod.sum <- mod.sum %>%
    dplyr::mutate(across(where(is.numeric), \(x) round(x, digits = digits)))
  
  print(digits)
  print(confint_level)
  
  if (!is.null(object$offset)){
    mod.sum <- tibble::add_row(mod.sum, parameter = paste(object$offset, "(Offset variable)"), coeff = 1, `Std. Err.` = NA_real_, `t-stat` = NA_real_, `p-value` = NA_real_, `lower CI` = NA_real_, `upper CI` = NA_real_)
  }
  
  print(mod.sum)
}
