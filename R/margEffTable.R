#' Marginal Effects, Elasticities, and Pseudo-Elasticities for flexCountReg Models
#'
#' Compute average marginal effects (AMEs), elasticities, or pseudo-elasticities
#' for fitted `flexCountReg` objects. Standard errors are computed using either
#' the delta method or a bootstrap, and the results can be returned as a tibble,
#' a gt table, or a LaTeX table.
#'
#' Continuous variables are handled by numerical differentiation on the
#' expected-response scale. Indicator variables are handled by discrete changes
#' (or pseudo-elasticities, depending on `indicator`).
#'
#' @param object A fitted `flexCountReg` object.
#' @param data Optional data frame. Defaults to the data used to fit the model.
#' @param vars Optional character vector of variable names to evaluate. If
#'   `NULL`, the function uses the non-response variables appearing in the
#'   model formula that are present in `data`.
#' @param measure Character scalar. One of:
#'   \describe{
#'     \item{"auto"}{Continuous variables receive AMEs; indicator variables
#'     receive pseudo-elasticities (or discrete changes if `indicator =
#'     "discrete_change"`).}
#'     \item{"ame"}{Compute average marginal effects for continuous variables.}
#'     \item{"elasticity"}{Compute elasticities for continuous variables.
#'     Indicator variables receive pseudo-elasticities (or discrete changes if
#'     `indicator = "discrete_change"`).}
#'   }
#' @param indicator Character scalar. For indicator variables, return
#'   `"pseudo_elasticity"` (default) or `"discrete_change"`.
#' @param se Character scalar. Standard error method: `"delta"` (default) or
#'   `"bootstrap"`.
#' @param bootstraps Integer. Number of bootstrap replications when `se =
#'   "bootstrap"`.
#' @param pred_method Optional prediction method forwarded to
#'   `predict.flexCountReg()`. For random-parameters models, this is commonly
#'   `"Exact"` or `"Simulated"`. If `NULL`, the function uses `"Exact"` for
#'   random-parameters models when possible.
#' @param confint_level Numeric scalar between 0 and 1. Confidence level for
#'   the confidence interval. Default is 0.95.
#' @param tableType Character scalar. One of `"tibble"`, `"gt"`, or `"latex"`.
#' @param digits Integer. Number of digits to round in the returned table.
#' @param cluster_var Optional character vector giving a clustering variable for
#'   bootstrap resampling. If supplied, bootstrap samples are drawn at the
#'   cluster level. If omitted, the function attempts row-level bootstrap
#'   resampling.
#' @param ... Additional arguments passed to `predict.flexCountReg()`.
#'
#' @details
#' The function works on the expected-response scale and is intended to be
#' compatible with all fitted `flexCountReg` model types, provided that
#' `predict.flexCountReg()` can generate predictions for the model.
#'
#' For continuous regressors, the average marginal effect is computed using a
#' central finite difference:
#' \deqn{
#' \text{AME}_j = \frac{1}{n}\sum_{i=1}^n \frac{\mu_i(x_{ij}+h)-\mu_i(x_{ij}-h)}{2h}
#' }
#' where \eqn{\mu_i} is the predicted mean response.
#'
#' For elasticities, the effect is computed as:
#' \deqn{
#' E_j = \frac{1}{n}\sum_{i=1}^n \left(\frac{x_{ij}}{\mu_i}\frac{\partial \mu_i}{\partial x_{ij}}\right)
#' }
#'
#' For indicator variables, the function computes a discrete change between the
#' observed baseline and a switched value. If `indicator = "pseudo_elasticity"`,
#' the result is reported as a percent change relative to the baseline.
#'
#' @returns A table object of the type requested by `tableType`. For
#'   `"tibble"`, the returned object also carries an `effect_info` attribute with
#'   metadata about the calculation.
#'
#' @examples
#' \donttest{
#' data("washington_roads")
#' washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
#'
#' nb2 <- countreg(
#'   Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'   data = washington_roads,
#'   family = "NB2"
#' )
#'
#' margEffTable(nb2, tableType = "tibble")
#' margEffTable(nb2, tableType = "gt")
#' }
#'
#' @export
margEffTable <- function(object,
                         data = NULL,
                         vars = NULL,
                         measure = c("auto", "ame", "elasticity"),
                         indicator = c("pseudo_elasticity", "discrete_change"),
                         se = c("delta", "bootstrap"),
                         bootstraps = 200,
                         pred_method = NULL,
                         confint_level = 0.95,
                         tableType = c("tibble", "gt", "latex"),
                         digits = 3,
                         cluster_var = NULL,
                         ...) {
  measure <- match.arg(measure)
  indicator <- match.arg(indicator)
  se <- match.arg(se)
  tableType <- match.arg(tableType)
  
  data <- .margEff_get_data(object, data)
  data <- as.data.frame(data)
  
  if (!is.numeric(confint_level) || length(confint_level) != 1L ||
      confint_level <= 0 || confint_level >= 1) {
    stop("`confint_level` must be a single numeric value between 0 and 1.")
  }
  
  if (!is.numeric(bootstraps) || length(bootstraps) != 1L || bootstraps < 1) {
    stop("`bootstraps` must be a positive integer.")
  }
  
  if (!is.numeric(digits) || length(digits) != 1L || digits < 0) {
    stop("`digits` must be a nonnegative integer.")
  }
  
  if (is.null(pred_method) && identical(.margEff_model_type(object), "countreg.rp")) {
    pred_method <- "Exact"
  }
  
  candidate_vars <- .margEff_default_vars(object, data)
  if (is.null(vars)) {
    vars <- candidate_vars
  } else {
    vars <- intersect(vars, names(data))
    missing_vars <- setdiff(vars, names(data))
    if (length(missing_vars) > 0L) {
      stop("The following variables are not in `data`: ", paste(missing_vars, collapse = ", "))
    }
  }
  
  if (length(vars) == 0L) {
    stop("No eligible variables were found for marginal-effect computation.")
  }
  
  eff_rows <- vector("list", length(vars))
  
  for (i in seq_along(vars)) {
    v <- vars[[i]]
    eff_rows[[i]] <- .margEff_one_term(
      object = object,
      data = data,
      var = v,
      measure = measure,
      indicator = indicator,
      se = se,
      bootstraps = bootstraps,
      pred_method = pred_method,
      confint_level = confint_level,
      cluster_var = cluster_var,
      ...
    )
  }
  
  out <- do.call(rbind, lapply(eff_rows, as.data.frame))
  rownames(out) <- NULL
  
  numeric_cols <- intersect(
    c("estimate", "std_error", "t_value", "p_value", "lower_ci", "upper_ci",
      "n_obs", "successful_bootstraps"),
    names(out)
  )
  
  out[numeric_cols] <- lapply(out[numeric_cols], function(x) round(as.numeric(x), digits))
  
  attr(out, "effect_info") <- list(
    measure = measure,
    indicator = indicator,
    se = se,
    bootstraps = bootstraps,
    confint_level = confint_level,
    pred_method = pred_method,
    model_type = .margEff_model_type(object),
    n_obs = nrow(data)
  )
  
  if (tableType == "tibble") {
    if (requireNamespace("tibble", quietly = TRUE)) {
      out <- tibble::as_tibble(out)
    }
    class(out) <- c("flexEffectTable", class(out))
    return(out)
  }
  
  if (tableType == "gt") {
    if (!requireNamespace("gt", quietly = TRUE)) {
      stop("Package `gt` is required for `tableType = \"gt\"`.")
    }
    gt_tbl <- gt::gt(out) |>
      gt::tab_header(
        title = "Marginal Effects / Elasticities",
        subtitle = paste0(
          "Method: ", se,
          if (!is.null(pred_method)) paste0(" | Prediction method: ", pred_method) else ""
        )
      ) |>
      gt::fmt_number(
        columns = c("estimate", "std_error", "t_value", "p_value", "lower_ci", "upper_ci"),
        decimals = digits
      ) |>
      gt::tab_source_note(
        source_note = "Significance codes: * p <= 0.05, ** p <= 0.01, *** p <= 0.001."
      )
    return(gt_tbl)
  }
  
  if (tableType == "latex") {
    if (!requireNamespace("knitr", quietly = TRUE)) {
      stop("Package `knitr` is required for `tableType = \"latex\"`.")
    }
    latex_tbl <- knitr::kable(
      out,
      format = "latex",
      booktabs = TRUE,
      caption = "Marginal Effects / Elasticities"
    )
    return(latex_tbl)
  }
  
  out
}

#' @export
print.flexEffectTable <- function(x, ...) {
  info <- attr(x, "effect_info")
  if (!is.null(info)) {
    cat("Marginal Effects / Elasticities\n")
    cat("Method:", info$se, "\n")
    if (!is.null(info$pred_method)) {
      cat("Prediction method:", info$pred_method, "\n")
    }
    cat("\n")
  }
  print.data.frame(as.data.frame(x), row.names = FALSE, ...)
  invisible(x)
}

.margEff_get_data <- function(object, data) {
  if (!is.null(data)) {
    return(data)
  }
  if (!is.null(object$data)) {
    return(object$data)
  }
  stop("No data supplied and no data stored in the fitted object.")
}

.margEff_model_type <- function(object) {
  if (!is.null(object$modelType)) {
    return(object$modelType)
  }
  if (!is.null(object$model) && !is.null(object$model$modelType)) {
    return(object$model$modelType)
  }
  "countreg"
}

.margEff_default_vars <- function(object, data) {
  form <- object$formula
  if (is.null(form)) {
    stop("The fitted object does not contain a formula.")
  }
  
  term_vars <- all.vars(stats::delete.response(stats::terms(form)))
  
  special_vars <- character(0)
  if (!is.null(object$model)) {
    if (!is.null(object$model$offset)) {
      special_vars <- c(special_vars, object$model$offset)
    }
    if (!is.null(object$model$panel_id)) {
      special_vars <- c(special_vars, object$model$panel_id)
    }
    if (!is.null(object$model$group_var)) {
      special_vars <- c(special_vars, object$model$group_var)
    }
  }
  
  term_vars <- setdiff(term_vars, special_vars)
  term_vars <- intersect(term_vars, names(data))
  term_vars
}

.margEff_is_indicator <- function(x) {
  if (is.factor(x)) {
    return(nlevels(x) == 2L)
  }
  if (is.logical(x)) {
    return(TRUE)
  }
  if (is.numeric(x) || is.integer(x)) {
    ux <- unique(stats::na.omit(x))
    if (length(ux) <= 2L && all(ux %in% c(0, 1))) {
      return(TRUE)
    }
  }
  FALSE
}

.margEff_indicator_levels <- function(x) {
  if (is.factor(x)) {
    if (nlevels(x) != 2L) {
      stop("Indicator variables must have exactly two levels.")
    }
    return(list(
      low = factor(rep(levels(x)[1], length(x)), levels = levels(x)),
      high = factor(rep(levels(x)[2], length(x)), levels = levels(x))
    ))
  }
  
  if (is.logical(x)) {
    return(list(
      low = rep(FALSE, length(x)),
      high = rep(TRUE, length(x))
    ))
  }
  
  if (is.numeric(x) || is.integer(x)) {
    return(list(
      low = rep(0, length(x)),
      high = rep(1, length(x))
    ))
  }
  
  stop("Unsupported indicator variable type.")
}

.margEff_predict <- function(object, newdata, pred_method = NULL, ...) {
  if (is.null(pred_method)) {
    pred <- predict(object, newdata = newdata, ...)
  } else {
    pred <- predict(object, newdata = newdata, method = pred_method, ...)
  }
  as.numeric(pred)
}

.margEff_set_theta <- function(object, theta) {
  out <- object
  theta <- as.numeric(theta)
  if (!is.null(names(object$model$estimate))) {
    names(theta) <- names(object$model$estimate)
  }
  out$model$estimate <- theta
  out
}

.margEff_theta <- function(object) {
  theta <- object$model$estimate
  theta <- as.numeric(unlist(theta, recursive = TRUE, use.names = TRUE))
  if (is.null(names(theta))) {
    est_names <- names(unlist(object$model$estimate, recursive = TRUE, use.names = TRUE))
    if (!is.null(est_names)) {
      names(theta) <- est_names
    }
  }
  theta
}

.margEff_vcov <- function(object) {
  model <- object$model
  
  if (!is.null(model$vcov) && is.matrix(model$vcov)) {
    return(model$vcov)
  }
  
  if (!is.null(model$vcov_matrix) && is.matrix(model$vcov_matrix)) {
    return(model$vcov_matrix)
  }
  
  if (!is.null(model$hessian) && is.matrix(model$hessian)) {
    vc <- tryCatch(solve(-model$hessian), error = function(e) NULL)
    if (!is.null(vc)) {
      return(vc)
    }
  }
  
  if (!is.null(model$bootstrapped_se)) {
    se <- model$bootstrapped_se
    if (is.data.frame(se) && all(c("term", "sd") %in% names(se))) {
      sd_vec <- as.numeric(se$sd)
      vc <- diag(sd_vec^2)
      dimnames(vc) <- list(se$term, se$term)
      return(vc)
    }
    if (is.numeric(se)) {
      vc <- diag(as.numeric(se)^2)
      if (!is.null(names(se))) {
        dimnames(vc) <- list(names(se), names(se))
      }
      return(vc)
    }
  }
  
  if (!is.null(model$se) && is.numeric(model$se)) {
    vc <- diag(as.numeric(model$se)^2)
    if (!is.null(names(model$se))) {
      dimnames(vc) <- list(names(model$se), names(model$se))
    }
    warning("Using a diagonal covariance approximation because the full covariance matrix was not available.")
    return(vc)
  }
  
  NULL
}

.margEff_effect_metric <- function(measure, is_indicator, indicator) {
  if (is_indicator) {
    if (indicator == "pseudo_elasticity") {
      return("Pseudo-elasticity (%)")
    }
    return("Discrete change")
  }
  
  if (measure == "elasticity") {
    return("Elasticity")
  }
  
  "Average marginal effect"
}

.margEff_one_term <- function(object,
                              data,
                              var,
                              measure,
                              indicator,
                              se,
                              bootstraps,
                              pred_method,
                              confint_level,
                              cluster_var,
                              ...) {
  x <- data[[var]]
  
  if (is.factor(x) && nlevels(x) > 2L) {
    stop(
      "Variable `", var, "` is a factor with more than two levels. ",
      "Create binary indicators or supply a specific contrast before calling `margEffTable()`."
    )
  }
  
  is_indicator <- .margEff_is_indicator(x)
  metric <- .margEff_effect_metric(measure, is_indicator, indicator)
  
  effect_fun <- function(theta) {
    obj2 <- .margEff_set_theta(object, theta)
    .margEff_effect_value(
      object = obj2,
      data = data,
      var = var,
      measure = measure,
      indicator = indicator,
      pred_method = pred_method,
      ...
    )$estimate
  }
  
  eff0 <- .margEff_effect_value(
    object = object,
    data = data,
    var = var,
    measure = measure,
    indicator = indicator,
    pred_method = pred_method,
    ...
  )
  
  est <- eff0$estimate
  n_obs <- nrow(data)
  
  if (se == "delta") {
    theta <- .margEff_theta(object)
    vc <- .margEff_vcov(object)
    if (is.null(vc)) {
      warning("No covariance matrix was available; delta-method standard errors could not be computed.")
      se_val <- NA_real_
    } else {
      if (length(theta) != nrow(vc)) {
        warning("Parameter vector and covariance matrix dimensions do not match; standard errors set to NA.")
        se_val <- NA_real_
      } else {
        grad <- .margEff_num_jacobian(effect_fun, theta)
        se_val <- as.numeric(sqrt(drop(grad %*% vc %*% t(grad))))
      }
    }
    successful_bootstraps <- NA_integer_
  } else {
    boot_out <- .margEff_bootstrap_se(
      object = object,
      data = data,
      var = var,
      measure = measure,
      indicator = indicator,
      bootstraps = bootstraps,
      pred_method = pred_method,
      cluster_var = cluster_var,
      ...
    )
    se_val <- boot_out$se
    successful_bootstraps <- boot_out$successful_bootstraps
  }
  
  t_val <- est / se_val
  p_val <- 2 * stats::pnorm(-abs(t_val))
  alpha <- 1 - confint_level
  z <- stats::qnorm(1 - alpha / 2)
  lower_ci <- est - z * se_val
  upper_ci <- est + z * se_val
  
  sig <- ifelse(
    is.na(p_val), "",
    ifelse(p_val <= 0.001, "***",
           ifelse(p_val <= 0.01, "**",
                  ifelse(p_val <= 0.05, "*", "")))
  )
  
  data.frame(
    term = var,
    variable_type = ifelse(is_indicator, "indicator", "continuous"),
    effect_metric = metric,
    estimate = est,
    std_error = se_val,
    t_value = t_val,
    p_value = p_val,
    sig = sig,
    lower_ci = lower_ci,
    upper_ci = upper_ci,
    n_obs = n_obs,
    successful_bootstraps = successful_bootstraps,
    stringsAsFactors = FALSE
  )
}

.margEff_effect_value <- function(object,
                                  data,
                                  var,
                                  measure,
                                  indicator,
                                  pred_method,
                                  ...) {
  x <- data[[var]]
  mu0 <- .margEff_predict(object, data, pred_method = pred_method, ...)
  
  is_indicator <- .margEff_is_indicator(x)
  
  if (is_indicator) {
    vals <- .margEff_indicator_levels(x)
    data_low <- data
    data_high <- data
    data_low[[var]] <- vals$low
    data_high[[var]] <- vals$high
    
    mu_low <- .margEff_predict(object, data_low, pred_method = pred_method, ...)
    mu_high <- .margEff_predict(object, data_high, pred_method = pred_method, ...)
    
    disc_change <- mu_high - mu_low
    
    if (indicator == "pseudo_elasticity") {
      est <- mean(100 * disc_change / pmax(mu_low, .Machine$double.eps), na.rm = TRUE)
    } else {
      est <- mean(disc_change, na.rm = TRUE)
    }
    
    return(list(
      estimate = est,
      variable_type = "indicator"
    ))
  }
  
  x_num <- as.numeric(x)
  h <- pmax(sqrt(.Machine$double.eps) * pmax(abs(x_num), 1), 1e-7)
  
  data_plus <- data
  data_minus <- data
  data_plus[[var]] <- x_num + h
  data_minus[[var]] <- x_num - h
  
  mu_plus <- .margEff_predict(object, data_plus, pred_method = pred_method, ...)
  mu_minus <- .margEff_predict(object, data_minus, pred_method = pred_method, ...)
  
  dmu_dx <- (mu_plus - mu_minus) / (2 * h)
  
  if (measure == "elasticity") {
    eff_i <- (x_num / pmax(mu0, .Machine$double.eps)) * dmu_dx
    est <- mean(eff_i, na.rm = TRUE)
  } else {
    est <- mean(dmu_dx, na.rm = TRUE)
  }
  
  list(
    estimate = est,
    variable_type = "continuous"
  )
}

.margEff_num_jacobian <- function(f, theta, eps = sqrt(.Machine$double.eps)) {
  theta <- as.numeric(theta)
  p <- length(theta)
  grad <- numeric(p)
  
  for (j in seq_len(p)) {
    step <- eps * max(abs(theta[j]), 1)
    if (!is.finite(step) || step == 0) {
      step <- eps
    }
    
    th_plus <- theta
    th_minus <- theta
    th_plus[j] <- th_plus[j] + step
    th_minus[j] <- th_minus[j] - step
    
    f_plus <- f(th_plus)
    f_minus <- f(th_minus)
    
    grad[j] <- (f_plus - f_minus) / (2 * step)
  }
  
  matrix(grad, nrow = 1)
}

.margEff_bootstrap_se <- function(object,
                                  data,
                                  var,
                                  measure,
                                  indicator,
                                  bootstraps,
                                  pred_method,
                                  cluster_var,
                                  ...) {
  successful <- 0L
  draws <- numeric(0)
  
  for (b in seq_len(bootstraps)) {
    boot_data <- .margEff_resample_data(data, cluster_var = cluster_var)
    
    fit <- tryCatch(
      .margEff_refit(object, boot_data),
      error = function(e) NULL
    )
    
    if (is.null(fit)) {
      next
    }
    
    eff <- tryCatch(
      .margEff_effect_value(
        object = fit,
        data = boot_data,
        var = var,
        measure = measure,
        indicator = indicator,
        pred_method = pred_method,
        ...
      )$estimate,
      error = function(e) NA_real_
    )
    
    if (is.finite(eff)) {
      successful <- successful + 1L
      draws <- c(draws, eff)
    }
  }
  
  if (length(draws) < 2L) {
    warning("Too few successful bootstrap replications to compute a bootstrap standard error.")
    return(list(se = NA_real_, successful_bootstraps = successful))
  }
  
  list(
    se = stats::sd(draws, na.rm = TRUE),
    successful_bootstraps = successful
  )
}

.margEff_resample_data <- function(data, cluster_var = NULL) {
  if (!is.null(cluster_var)) {
    cluster_var <- cluster_var[cluster_var %in% names(data)]
  }
  
  if (!is.null(cluster_var) && length(cluster_var) > 0L) {
    if (length(cluster_var) == 1L) {
      ids <- data[[cluster_var]]
    } else {
      ids <- interaction(data[cluster_var], drop = TRUE, sep = "__")
    }
    
    unique_ids <- unique(ids)
    sampled_ids <- sample(unique_ids, size = length(unique_ids), replace = TRUE)
    
    idx <- unlist(lapply(sampled_ids, function(id) which(ids == id)), use.names = FALSE)
    out <- data[idx, , drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  
  idx <- sample.int(nrow(data), size = nrow(data), replace = TRUE)
  out <- data[idx, , drop = FALSE]
  rownames(out) <- NULL
  out
}

.margEff_refit <- function(object, boot_data) {
  fit_call <- object$call
  if (is.null(fit_call)) {
    stop("No stored call is available for bootstrap refitting.")
  }
  
  fit_call$data <- boot_data
  eval(fit_call, envir = parent.frame())
}