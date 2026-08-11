#' Estimate a Random Effects Poisson Panel Model
#'
#' Estimate a random effects Poisson panel model for clustered count data.
#' The model is also known as a negative multinomial model.
#'
#' @name repois
#' @param formula an R formula.
#' @param group_var the grouping variable(s) for the random effects (e.g.,
#'   individual ID or other panel ID variables).
#' @param data a dataframe that has all of the variables in the \code{formula}.
#' @param offset an optional offset term provided as a string.
#' @param method a method to use for optimization in the maximum likelihood
#'   estimation. For options, see \code{\link[maxLik]{maxLik}}.
#' @param max.iters the maximum number of iterations to allow the optimization
#'   method to perform.
#' @param print.level integer specifying the verbosity of output during
#'   optimization.
#' @param bootstraps optional integer specifying the number of bootstrap samples
#'   to be used for estimating standard errors. If not specified, no
#'   bootstrapping is performed.
#'
#' @import maxLik stats modelr tibble
#' @importFrom purrr map map_df
#' @importFrom dplyr group_by reframe
#'
#' @details
#' This function estimates a random effects Poisson panel model. The likelihood
#' is integrated analytically over a gamma-distributed cluster effect with mean
#' 1 and variance \code{alpha}. The model estimates the regression coefficients
#' and one heterogeneity parameter \code{ln(alpha)}.
#' 
#' The PMF of this distribution is:
#' \deqn{P(y_{it} \mid \mu_{it}) =
#' \frac{\Gamma\left(\sum_t y_{it} + \frac{1}{\alpha}\right)}
#' {\Gamma\left(\frac{1}{\alpha}\right)\prod_t y_{it}!}
#' \left(\frac{1}{\alpha\sum_t \mu_{it} + 1}\right)^{\frac{1}{\alpha}}
#' \prod_t \left(\frac{\alpha\mu_{it}}{\alpha\sum_t \mu_{it} + 1}
#' \right)^{y_{it}}}
#'
#' @returns
#' An object of class \code{countreg} which is a list with the following
#' components:
#' \itemize{
#' \item model: the fitted model object.
#' \item data: the data frame used to fit the model.
#' \item call: the matched call.
#' \item formula: the formula used to fit the model.
#' }
#'
#' @examples
#' \donttest{
#' ## RE Poisson Panel Model
#' data("washington_roads")
#' washington_roads$AADTover10k <-
#'   ifelse(washington_roads$AADT > 10000, 1, 0)
#'
#' repois.mod <- repois(
#'   Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
#'   data = washington_roads,
#'   offset = "lnlength",
#'   group_var = "ID",
#'   method = "BFGS",
#'   max.iters = 1000
#' )
#'
#' summary(repois.mod)
#' }
#'
#' @export
repois <- function(formula,
                   group_var,
                   data,
                   method = "NM",
                   max.iters = 1000,
                   print.level = 0,
                   bootstraps = NULL,
                   offset = NULL) {
  
  call <- match.call()
  
  prep <- .repois_prepare_inputs(
    formula = formula,
    group_var = group_var,
    data = data,
    offset = offset
  )
  
  # Keep local copies to avoid R CMD check notes about X/y.
  X <- prep$X
  y <- prep$y
  
  fit <- .repois_fit_engine(
    y = y,
    X = X,
    offset_vec = prep$offset_vec,
    cluster_list = prep$cluster_list,
    method = method,
    max.iters = max.iters,
    print.level = print.level,
    bootstraps = bootstraps
  )
  
  # Keep the same object fields your summary/predict methods already use.
  fit$formula <- formula
  fit$observed <- y
  fit$predictions <- fit$fitted.values
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum
  fit$modelType <- "REPoisson"
  fit$offset <- offset
  fit$call <- call
  
  obj <- .createFlexCountReg(
    model = fit,
    data = data,
    call = call,
    formula = formula
  )
  
  obj
}

.repois_prepare_inputs <- function(formula, group_var, data, offset) {
  if (!inherits(data, "data.frame")) {
    stop("`data` must be a data frame.")
  }
  if (missing(group_var) || length(group_var) < 1L) {
    stop("`group_var` must be supplied.")
  }
  if (is.null(names(data))) {
    stop("`data` must have named columns.")
  }
  
  group_var <- as.character(group_var)
  missing_group <- setdiff(group_var, names(data))
  if (length(missing_group) > 0L) {
    stop("Grouping variable(s) not found in `data`: ",
         paste(missing_group, collapse = ", "))
  }
  
  data_work <- data
  
  # Keep only rows with complete grouping information before model.frame().
  group_ok <- stats::complete.cases(data_work[, group_var, drop = FALSE])
  data_work <- data_work[group_ok, , drop = FALSE]
  
  # Reset row names so model.frame() row indices are easy to recover.
  rownames(data_work) <- seq_len(nrow(data_work))
  
  formula_aug <- formula
  if (!is.null(offset)) {
    if (!is.character(offset) || length(offset) != 1L) {
      stop("`offset` must be a single character string naming a variable.")
    }
    if (!offset %in% names(data_work)) {
      stop("Offset variable not found in `data`: ", offset)
    }
    if (grepl("offset\\s*\\(", paste(deparse(formula), collapse = " "))) {
      stop("Specify offset either in `formula` or via `offset`, not both.")
    }
    
    data_work$.__repois_offset__ <- data_work[[offset]]
    formula_aug <- stats::as.formula(
      paste(
        paste(deparse(formula), collapse = " "),
        "+ offset(.__repois_offset__)"
      )
    )
  }
  
  mf <- stats::model.frame(formula_aug, data = data_work, na.action = na.omit)
  y <- stats::model.response(mf)
  
  if (!is.numeric(y)) {
    stop("The response variable must be numeric.")
  }
  if (any(y < 0, na.rm = TRUE)) {
    stop("The response variable must be nonnegative.")
  }
  if (any(abs(y - round(y)) > 1e-8, na.rm = TRUE)) {
    warning("The response variable is not integer-valued; Poisson counts are typically integers.")
  }
  
  X <- stats::model.matrix(stats::terms(mf), mf)
  
  offset_vec <- stats::model.offset(mf)
  if (is.null(offset_vec)) {
    offset_vec <- rep(0, length(y))
  }
  
  row_idx <- as.integer(rownames(mf))
  data_used <- data_work[row_idx, , drop = FALSE]
  
  group_vals <- lapply(group_var, function(v) data_used[[v]])
  if (length(group_vals) == 1L) {
    group_fac <- as.factor(group_vals[[1L]])
  } else {
    group_fac <- interaction(group_vals, drop = TRUE, sep = ":")
  }
  
  cluster_list <- split(seq_len(length(y)), group_fac, drop = TRUE)
  cluster_list <- cluster_list[lengths(cluster_list) > 0L]
  
  if (length(cluster_list) < 1L) {
    stop("No non-empty groups remain after removing missing values.")
  }
  
  list(
    y = y,
    X = X,
    offset_vec = offset_vec,
    cluster_list = cluster_list,
    data_used = data_used,
    model_frame = mf,
    terms = stats::terms(mf),
    xlevels = .getXlevels(stats::terms(mf), mf),
    contrasts = attr(X, "contrasts")
  )
}

.repois_map_method <- function(method) {
  if (is.null(method) || length(method) != 1L) {
    stop("`method` must be a single string.")
  }
  
  method <- toupper(method)
  
  switch(
    method,
    "NM" = "NM",
    "NELDER-MEAD" = "NM",
    "BFGS" = "BFGS",
    "NR" = "NR",
    "BHHH" = "BHHH",
    "CG" = "CG",
    "SANN" = "SANN",
    stop("Unsupported optimization method: ", method)
  )
}

.repois_start_values <- function(y, X, offset_vec) {
  fit0 <- try(
    stats::glm.fit(
      x = X,
      y = y,
      family = stats::poisson(),
      offset = offset_vec
    ),
    silent = TRUE
  )
  
  if (inherits(fit0, "try-error")) {
    beta0 <- rep(0, ncol(X))
    names(beta0) <- colnames(X)
  } else {
    beta0 <- fit0$coefficients
    beta0[is.na(beta0)] <- 0
    if (length(beta0) != ncol(X)) {
      beta0 <- rep(0, ncol(X))
      names(beta0) <- colnames(X)
    }
  }
  
  eta0 <- drop(X %*% beta0) + offset_vec
  mu0 <- exp(pmin(eta0, 700))
  mu0[!is.finite(mu0)] <- 1
  
  pearson_phi <- sum((y - mu0)^2 / pmax(mu0, .Machine$double.eps)) /
    max(length(y) - ncol(X), 1L)
  
  alpha0 <- max(1e-6, (pearson_phi - 1) / max(mean(mu0), .Machine$double.eps))
  lnalpha0 <- log(alpha0)
  
  c(beta0, lnalpha0)
}

.repois_loglik <- function(par, y, X, offset_vec, cluster_list) {
  p <- ncol(X)
  beta <- par[seq_len(p)]
  lnalpha <- par[p + 1L]
  alpha <- exp(lnalpha)
  
  if (!is.finite(alpha) || alpha <= 0) {
    return(-Inf)
  }
  
  eta <- drop(X %*% beta) + offset_vec
  if (any(!is.finite(eta)) || any(eta > 700)) {
    return(-Inf)
  }
  
  ll <- 0
  for (idx in cluster_list) {
    y_i <- y[idx]
    eta_i <- eta[idx]
    mu_i <- exp(eta_i)
    
    ysum <- sum(y_i)
    s <- sum(mu_i)
    
    ll_i <- lgamma(ysum + 1 / alpha) -
      lgamma(1 / alpha) +
      ysum * log(alpha) -
      (ysum + 1 / alpha) * log1p(alpha * s) +
      sum(y_i * eta_i) -
      sum(lgamma(y_i + 1))
    
    if (!is.finite(ll_i)) {
      return(-Inf)
    }
    
    ll <- ll + ll_i
  }
  
  ll
}

.repois_fit_engine <- function(y,
                               X,
                               offset_vec,
                               cluster_list,
                               method = "NM",
                               max.iters = 1000,
                               print.level = 0,
                               bootstraps = NULL,
                               compute_vcov = TRUE) {
  if (!requireNamespace("maxLik", quietly = TRUE)) {
    stop("Package `maxLik` is required for `repois()`.")
  }
  
  method_opt <- .repois_map_method(method)
  start <- .repois_start_values(y, X, offset_vec)
  
  control <- maxLik::maxControl(
    iterlim = max.iters,
    printLevel = print.level
  )
  
  fit_ml <- maxLik::maxLik(
    logLik = .repois_loglik,
    start = start,
    method = method_opt,
    control = control,
    y = y,
    X = X,
    offset_vec = offset_vec,
    cluster_list = cluster_list
  )
  
  par <- fit_ml$estimate
  names(par) <- c(colnames(X), "lnalpha")
  
  p <- ncol(X)
  beta <- par[seq_len(p)]
  lnalpha <- par[p + 1L]
  alpha <- exp(lnalpha)
  
  eta <- drop(X %*% beta) + offset_vec
  mu <- exp(pmin(eta, 700))
  
  logLik_value <- as.numeric(fit_ml$maximum)
  k <- length(par)
  nobs <- length(y)
  n_groups <- length(cluster_list)
  
  vcov_mat <- matrix(
    NA_real_,
    nrow = k,
    ncol = k,
    dimnames = list(names(par), names(par))
  )
  
  se <- rep(NA_real_, k)
  z <- rep(NA_real_, k)
  pval <- rep(NA_real_, k)
  names(se) <- names(z) <- names(pval) <- names(par)
  
  if (isTRUE(compute_vcov)) {
    if (!is.null(bootstraps) && length(bootstraps) == 1L && bootstraps > 1L) {
      boot_vcov <- .repois_bootstrap_vcov(
        y = y,
        X = X,
        offset_vec = offset_vec,
        cluster_list = cluster_list,
        method = method,
        max.iters = max.iters,
        print.level = print.level,
        B = bootstraps
      )
      if (!is.null(boot_vcov)) {
        vcov_mat <- boot_vcov
      }
    }
    
    if (all(is.na(vcov_mat))) {
      vc_try <- try(stats::vcov(fit_ml), silent = TRUE)
      if (!inherits(vc_try, "try-error") && all(is.finite(vc_try))) {
        vcov_mat <- vc_try
      }
    }
    
    if (all(is.finite(vcov_mat))) {
      se <- sqrt(pmax(diag(vcov_mat), 0))
      z <- par / se
      pval <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    }
  }
  
  out <- list(
    par = par,
    estimate = par,
    coefficients = par,
    beta = beta,
    lnalpha = lnalpha,
    alpha = alpha,
    vcov = vcov_mat,
    se = se,
    z = z,
    p = pval,
    logLik_value = logLik_value,
    aic = -2 * logLik_value + 2 * k,
    bic = -2 * logLik_value + log(nobs) * k,
    fitted.values = mu,
    linear.predictors = eta,
    residuals = y - mu,
    convergence = fit_ml$code,
    message = fit_ml$message,
    iterations = if (!is.null(fit_ml$iterations)) fit_ml$iterations else NA_integer_,
    counts = fit_ml$counts,
    nobs = nobs,
    n_groups = n_groups,
    method = method_opt,
    hessian = fit_ml$hessian,
    gradient = fit_ml$gradient,
    maximum = fit_ml$maximum,
    maxLik_object = fit_ml
  )
  
  if (isTRUE(compute_vcov) && all(is.finite(vcov_mat))) {
    out$se <- sqrt(pmax(diag(vcov_mat), 0))
  }
  
  if (!is.null(bootstraps) && length(bootstraps) == 1L && bootstraps > 1L) {
    out$bootstrapped_se <- .repois_bootstrap_se(
      y = y,
      X = X,
      offset_vec = offset_vec,
      cluster_list = cluster_list,
      method = method,
      max.iters = max.iters,
      print.level = print.level,
      B = bootstraps,
      start = par
    )
    out$successful_bootstraps <- attr(out$bootstrapped_se, "successful_bootstraps")
  } else {
    out$bootstrapped_se <- NULL
    out$successful_bootstraps <- 0L
  }
  
  out
}

.repois_bootstrap_vcov <- function(y,
                                   X,
                                   offset_vec,
                                   cluster_list,
                                   method,
                                   max.iters,
                                   print.level,
                                   B) {
  boot_se <- .repois_bootstrap_se(
    y = y,
    X = X,
    offset_vec = offset_vec,
    cluster_list = cluster_list,
    method = method,
    max.iters = max.iters,
    print.level = print.level,
    B = B,
    start = NULL
  )
  
  if (is.null(boot_se) || nrow(boot_se) < 2L) {
    return(NULL)
  }
  
  stats::cov(boot_se[, -1L, drop = FALSE])
}

.repois_bootstrap_se <- function(y,
                                 X,
                                 offset_vec,
                                 cluster_list,
                                 method,
                                 max.iters,
                                 print.level,
                                 B,
                                 start = NULL) {
  par_len <- ncol(X) + 1L
  estimates <- matrix(NA_real_, nrow = B, ncol = par_len)
  
  cluster_ids <- seq_along(cluster_list)
  ok <- 0L
  
  if (is.null(start)) {
    start <- .repois_start_values(y, X, offset_vec)
  }
  
  for (b in seq_len(B)) {
    sampled_clusters <- sample(cluster_ids, size = length(cluster_ids), replace = TRUE)
    sampled_rows_list <- cluster_list[sampled_clusters]
    
    boot_idx <- unlist(sampled_rows_list, use.names = FALSE)
    boot_y <- y[boot_idx]
    boot_X <- X[boot_idx, , drop = FALSE]
    boot_offset <- offset_vec[boot_idx]
    
    boot_group <- rep(seq_along(sampled_rows_list), times = lengths(sampled_rows_list))
    boot_cluster_list <- split(seq_len(length(boot_idx)), boot_group)
    
    boot_fit <- try(
      .repois_fit_engine(
        y = boot_y,
        X = boot_X,
        offset_vec = boot_offset,
        cluster_list = boot_cluster_list,
        method = method,
        max.iters = max.iters,
        print.level = print.level,
        bootstraps = NULL,
        compute_vcov = FALSE
      ),
      silent = TRUE
    )
    
    if (!inherits(boot_fit, "try-error") && is.finite(boot_fit$logLik_value)) {
      ok <- ok + 1L
      estimates[ok, ] <- boot_fit$par
    }
  }
  
  if (ok < 2L) {
    warning("Too few successful bootstrap replications; falling back to Hessian-based standard errors.")
    return(NULL)
  }
  
  estimates <- estimates[seq_len(ok), , drop = FALSE]
  se_df <- data.frame(
    term = c(colnames(X), "ln(alpha)"),
    sd = apply(estimates, 2L, stats::sd)
  )
  
  attr(se_df, "successful_bootstraps") <- ok
  se_df
}


.coef.repois <- function(object, ...) {
  object$coefficients
}


.vcov.repois <- function(object, ...) {
  object$vcov
}


.logLik.repois <- function(object, ...) {
  structure(
    object$logLik_value,
    class = "logLik",
    nobs = object$nobs,
    df = length(object$coefficients)
  )
}


.fitted.repois <- function(object, ...) {
  object$fitted.values
}


`%||%` <- function(a, b) if (!is.null(a)) a else b
