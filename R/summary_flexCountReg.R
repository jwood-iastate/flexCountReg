#' Custom summary method for flexCountReg models
#' 
#' @param object A flexCountReg model object.
#' @param ... Optional parameters that include `confint_level` and `digits`.
#' 
#' @note Optional parameter `confint_level`: A numeric value between 0 and 1
#'   indicating the confidence level for confidence intervals. Default is 0.95.
#' @note Optional parameter `digits`: Number of digits (decimal places) to round
#'   to. Default is 3.
#' 
#' @import tibble 
#' @importFrom dplyr mutate across  %>% where
#' @details
#' This summary method accounts for bootstrapped or robust standard errors (when
#' used).
#' 
#' @returns Prints the model formula, method used for estimation, number of 
#'  iterations used, if the model converged, and the log-likelihood. Then, it 
#'  prints a table containing parameter estimates, standard errors, 
#'  t-statistics, p-values, and confidence intervals. Also quietly returns a 
#'  tibble with these values.
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
  if (!is.null(args$confint_level) && is.numeric(args$confint_level)) {
    confint_level <- args$confint_level
    } else {
      confint_level <- 0.95
  }
  if (!is.null(args$digits) && is.numeric(args$digits)) {
    digits <- args$digits
    } else {
      digits <- 3  
  }
  
  model_obj <- object$model
  
  cat("Call:\n", deparse(model_obj$formula), "\n")
  if(!is.null(model_obj$modelType)) 
    cat("\n", "Method: ", model_obj$modelType, "\n")
  if(!is.null(model_obj$iterations)) 
    cat("Iterations: ", model_obj$iterations, "\n")
  if(!is.null(model_obj$message)) 
    cat("Convergence: ", model_obj$message, "\n")
  if(!is.null(model_obj$maximum)) 
    cat("Log-likelihood: ", model_obj$maximum, "\n")
  
  if (!is.null(model_obj$vcov_type)) {
    se_label <- switch(model_obj$vcov_type,
                       hessian = "Hessian-based",
                       robust = "robust (sandwich)",
                       bootstrap = paste0("bootstrapped (",
                                          model_obj$n_boot_ok, " of ",
                                          model_obj$bootstraps,
                                          " replicates)"),
                       model_obj$vcov_type)
    cat("Standard errors: ", se_label, "\n")
  }
  if (!is.null(model_obj$n_dropped) && model_obj$n_dropped > 0) {
    cat("Rows removed due to missing values: ", model_obj$n_dropped, "\n")
  }
  
  cat("\nParameter Estimates:\n")
  
  # Prepare Coefficients
  params <- names(model_obj$estimate)
  coeffs <- as.numeric(model_obj$estimate)
  se_vec <- rep(NA_real_, length(coeffs))
  
  # Extract standard errors, always matched to the parameters by name
  if (!is.null(model_obj$vcov) && is.matrix(model_obj$vcov)) {
    se_all <- sqrt(pmax(diag(model_obj$vcov), 0))
    names(se_all) <- rownames(model_obj$vcov)
  } else if (is.data.frame(model_obj$bootstrapped_se) &&
             all(c("term", "sd") %in% names(model_obj$bootstrapped_se))) {
    se_all <- stats::setNames(as.numeric(model_obj$bootstrapped_se$sd),
                              model_obj$bootstrapped_se$term)
  } else if (!is.null(model_obj$bootstrapped_se)) {
    se_all <- model_obj$bootstrapped_se
  } else {
    se_all <- model_obj$se
  }
  if (!is.null(se_all)) {
    if (!is.null(names(se_all)) && all(params %in% names(se_all))) {
      se_vec <- as.numeric(se_all[params])
    } else if (length(se_all) == length(coeffs)) {
      se_vec <- as.numeric(se_all)
    }
  }
  
  # Final length check: if SE vector doesn't match coeffs, force all to NA to
  # prevent crash
  if(length(se_vec) != length(coeffs)) {
    se_vec <- rep(NA_real_, length(coeffs))
  }
  
  # Standardize invalid SEs (NaN, Inf) to NA to ensure propagation works as
  # requested
  # Note: is.finite() fails on non-numeric types, so we check numeric status
  # implicit in se_vec definition above
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
  if (!is.null(model_obj[["offset"]]) && length(model_obj[["offset"]]) > 0){
    # Handle multiple offsets case
    offset_name <- if (length(model_obj[["offset"]]) > 1) {
      "Offset" 
      } else { 
        model_obj[["offset"]]
      }
    
    mod.sum <- tibble::add_row(
      mod.sum, 
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