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
#' \donttest{
#' # NB2 Model
#' data("washington_roads")
#' washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
#' nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'                 data = washington_roads, family = "NB2",
#'                 dis_param_formula_1 = ~ speed50, method='BFGS')
#' summary(nb2)
#' }
#' 
#' @export
summary.flexCountReg <- function(object, ...) {
  # Extract optional parameters from '...'
  args <- list(...)
  if (!is.null(args$confint_level) && is.numeric(args$confint_level)) confint_level <- args$confint_level else confint_level <- 0.95  
  if (!is.null(args$digits) && is.numeric(args$digits)) digits <- args$digits else digits <- 3  
  
  model_obj <- object$model
  
  cat("Call:\n", deparse(model_obj$formula), "\n")
  if(!is.null(model_obj$modelType)) cat("\n", "Method: ", model_obj$modelType, "\n")
  if(!is.null(model_obj$iterations)) cat("Iterations: ", model_obj$iterations, "\n")
  if(!is.null(model_obj$message)) cat("Convergence: ", model_obj$message, "\n")
  if(!is.null(model_obj$maximum)) cat("Log-likelihood: ", model_obj$maximum, "\n")
  
  cat("\nParameter Estimates:\n")
  
  # Prepare Coefficients
  params <- names(model_obj$estimate)
  coeffs <- as.numeric(model_obj$estimate)
  se_vec <- rep(NA_real_, length(coeffs))
  
  # Extract Standard Errors Logic
  using_bootstrap <- FALSE
  
  if (!is.null(model_obj$bootstrapped_se)) {
    # Check if it's a data frame (common output from broom logic in tests) or vector
    if (is.data.frame(model_obj$bootstrapped_se)) {
      # Try to match by name if 'term' and 'sd' columns exist
      if(all(c("term", "sd") %in% names(model_obj$bootstrapped_se))){
        match_idx <- match(params, model_obj$bootstrapped_se$term)
        # FIX: Explicitly convert to numeric to prevent type issues (e.g., list columns or attributes)
        se_vec <- as.numeric(model_obj$bootstrapped_se$sd[match_idx])
      } else {
        # Fallback: if lengths match, assume specific column order (2nd col)
        if(nrow(model_obj$bootstrapped_se) == length(coeffs)){
          se_vec <- as.numeric(model_obj$bootstrapped_se[[2]])
        }
      }
    } else {
      # It is a vector
      se_vec <- as.numeric(model_obj$bootstrapped_se)
    }
    using_bootstrap <- TRUE
    cat("(Using bootstrapped standard errors)\n")
    
  } else if (!is.null(model_obj$se)) {
    se_vec <- as.numeric(model_obj$se)
  }
  
  # Final length check: if SE vector doesn't match coeffs, force all to NA to prevent crash
  if(length(se_vec) != length(coeffs)) {
    se_vec <- rep(NA_real_, length(coeffs))
  }
  
  # Standardize invalid SEs (NaN, Inf) to NA to ensure propagation works as requested
  # Note: is.finite() fails on non-numeric types, so we check numeric status implicit in se_vec definition above
  se_vec[!is.finite(se_vec)] <- NA_real_
  
  # Create tibble with explicit numeric columns
  mod.sum <- tibble::tibble(parameter = params, 
                            coeff = coeffs, 
                            `Std. Err.` = se_vec)
  
  # Calculate T-stats
  # Division by NA results in NA. 
  mod.sum$`t-stat` <- mod.sum$coeff / mod.sum$`Std. Err.`
  
  # Calculate P-values
  # pnorm(NA) results in NA.
  mod.sum$`p-value` <- 2 * pnorm(-1 * abs(mod.sum$`t-stat`), mean = 0, sd = 1)
  
  # Calculate Confidence Intervals
  lower.level <- (1 - confint_level) / 2
  upper.level <- 1 - lower.level
  
  ci_z.lower <- qnorm(lower.level)
  ci_z.upper <- qnorm(upper.level)
  
  mod.sum$`lower CI` <- mod.sum$coeff + ci_z.lower * mod.sum$`Std. Err.`
  mod.sum$`upper CI` <- mod.sum$coeff + ci_z.upper * mod.sum$`Std. Err.`
  
  # Rounding
  mod.sum <- mod.sum %>%
    dplyr::mutate(across(where(is.numeric), \(x) round(x, digits = digits)))
  
  # Add Offset row if present
  if (!is.null(model_obj$offset)){
    # Handle multiple offsets case
    offset_name <- if(length(model_obj$offset) > 1) "Offset" else model_obj$offset
    
    mod.sum <- tibble::add_row(mod.sum, 
                               parameter = paste(offset_name, "(Offset variable)"), 
                               coeff = 1, 
                               `Std. Err.` = NA_real_, 
                               `t-stat` = NA_real_, 
                               `p-value` = NA_real_, 
                               `lower CI` = NA_real_, 
                               `upper CI` = NA_real_)
  }
  
  print(mod.sum)
  invisible(mod.sum)
}