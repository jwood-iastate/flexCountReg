#' Estimate a Random Effects Poisson Panel Model
#'
#' Estimate a random effects Poisson panel model for clustered count data.
#' The model is also known as a negative multinomial model. A more flexible 
#' version of this model can be estimated using the `negmultinom()` function.
#'
#' @name repois
#' @param formula an R formula.
#' @param group_var the grouping variable(s) for the random effects (e.g.,
#'   individual ID or other panel ID variables).
#' @param data a dataframe that has all of the variables in the \code{formula}.
#' @param offset an optional offset term provided as a string.
#' @param sp optional fixed smoothing parameters for smooth terms
#'   (\code{s()}, \code{te()}, ...) in the formula(s); by default they are
#'   selected automatically. Smooth terms are fitted by penalized maximum
#'   likelihood with Bayesian standard errors (see the "Smooth terms"
#'   section of \code{\link{countreg}}); bootstrap standard errors are not
#'   available with them.
#' @param weights the name of a column of `data`, or a numeric vector with one
#'   value per row of `data`, of frequency weights of the panels: a panel with
#'   weight w counts as w independent panels in the log-likelihood, the robust
#'   covariance, BIC and the bootstrap. Weights must be constant within each
#'   panel. Bootstrap standard errors need whole-number weights.
#' @param method a method to use for optimization in the maximum likelihood
#'   estimation. For options, see \code{\link[maxLik]{maxLik}}.
#' @param max.iters the maximum number of iterations to allow the optimization
#'   method to perform.
#' @param verbose logical; if \code{TRUE}, print the progress of the
#'   optimization.
#' @param stderr the type of standard errors: \code{"normal"} (Hessian-based,
#'   the default), \code{"robust"} (sandwich, clustered by group), or
#'   \code{"boot"} (bootstrap that resamples whole groups). Giving
#'   \code{bootstraps} without \code{stderr} also requests bootstrapped
#'   standard errors.
#' @param keep_data logical; if \code{FALSE}, the training data are not stored
#'   in the fitted object.
#' @param cre optional one-sided formula of time-varying regressors, e.g.
#'   \code{~ x1 + x2}, for correlated random effects (Mundlak 1978): the
#'   group means of these variables (columns \code{mean_x1}, ...) are added
#'   to the mean model, so the random effect may depend on the regressors
#'   through their means. The means are computed on the estimation rows and
#'   stored with the fit: predictions for a group in the estimation data use
#'   its training means, and a new group uses the means of its own rows in
#'   \code{newdata}. \code{\link{cre_test}()} tests random against fixed
#'   effects.
#' @param print.level deprecated; use \code{verbose}.
#' @param bootstraps optional integer specifying the number of bootstrap samples
#'   to be used for estimating standard errors. If not specified, no
#'   bootstrapping is performed.
#'
#' @importFrom maxLik maxLik
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
                   verbose = FALSE,
                   stderr = "normal",
                   bootstraps = NULL,
                   offset = NULL,
                   weights = NULL,
                   keep_data = TRUE,
                   print.level = NULL,
                   cre = NULL,
                   sp = NULL) {
  
  call <- match.call()
  print.level <- .fcr_print_level("repois", verbose, print.level)
  stderr <- .fcr_resolve_stderr(stderr, missing(stderr), bootstraps)
  if (stderr == "boot" && is.null(bootstraps)) {
    bootstraps <- 200L
    message("`bootstraps` not given; using 200 bootstrap samples.")
  }
  if (stderr != "boot") bootstraps <- NULL
  
  method <- .fcr_check_method(method)
  .fcr_check_controls(max.iters)
  if (missing(group_var) || is.null(group_var)) {
    stop("`group_var` must name the column(s) that identify the groups.",
         call. = FALSE)
  }
  
  # Correlated random effects: group means of time-varying regressors
  cre_info <- NULL
  if (!is.null(cre)) {
    cp <- .fcr_cre_prepare(cre, formula, as.data.frame(data), group_var,
                           extra_vars = c(if (is.character(offset)) offset,
                                          if (is.character(weights)) weights))
    data <- cp$data
    formula <- cp$formula
    cre_info <- cp$info
  }

  # Shared data preparation: one complete-case frame (formula, offsets and
  # groups), training blueprint and the rows used
  prep <- prepare_model_data(formula, data, offset = offset,
                             group_var = group_var, weights = weights)
  formula <- prep$formula   # with `.` expanded
  if (prep$n_dropped > 0) {
    message(prep$n_dropped, " row(s) with missing values were removed.")
  }
  y <- as.numeric(prep$y)
  X <- prep$X
  cluster_list <- split(seq_along(y), prep$group, drop = TRUE)
  # Panel frequency weights (NULL without `weights`)
  w_panel <- .fcr_panel_weights(prep, weights)
  if (stderr == "boot" && !is.null(w_panel)) {
    .fcr_check_boot_freq(w_panel, length(w_panel))
  }
  x_names <- c(colnames(X), "lnalpha")
  
  has_smooth <- .fcr_prep_has_smooth(prep)
  if (has_smooth) .fcr_smooth_no_boot(stderr)
  fit <- .repois_fit_engine(
    y = y,
    X = X,
    offset_vec = prep$offset,
    cluster_list = cluster_list,
    method = method,
    max.iters = max.iters,
    print.level = print.level,
    bootstraps = NULL,
    w = w_panel,
    smooth = if (has_smooth) list(ents = list(main = prep$smooth$main),
                                  sp = sp)
  )
  names(fit$estimate) <- x_names
  
  fit$formula <- formula
  fit$observed <- y
  fit$y <- y
  fit$predictions <- fit$fitted.values
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum
  fit$modelType <- "REPoisson"
  fit$offset <- prep$blueprint$offset_label
  fit$call <- call
  fit$beta_pred <- unname(fit$estimate[seq_len(ncol(X))])
  fit$stderr <- stderr
  fit$blueprint <- prep$blueprint
  fit$rows_used <- prep$rows_used
  fit$n_dropped <- prep$n_dropped
  fit$nobs <- length(y)
  fit$n_groups <- length(cluster_list)
  fit$weights <- prep$weights
  fit$panel_index <- as.integer(prep$group)
  fit$lr_sample <- .fcr_lr_sample(prep, data, group_var)
  fit$panel_weights <- w_panel
  fit$n_bic <- if (is.null(w_panel)) length(cluster_list) else sum(w_panel)
  
  if (stderr == "robust") {
    # One log-likelihood contribution per cluster: cluster-robust sandwich
    V <- .fcr_robust_vcov(fit$maxLik_object, w_panel)
  } else if (stderr == "boot") {
    # Resample whole groups and refit from the full-sample estimates
    bs <- .fcr_cluster_bootstrap(prep, as.integer(bootstraps), x_names,
                                 refit = function(sub) {
      f <- .repois_fit_engine(
        y = as.numeric(sub$y), X = sub$X, offset_vec = sub$offset,
        cluster_list = split(seq_along(sub$y), sub$group, drop = TRUE),
        method = method, max.iters = max.iters, print.level = 0,
        bootstraps = NULL, compute_vcov = FALSE, w = sub$panel_w)
      .fcr_refit_result(stats::setNames(f$par, x_names),
                        .fcr_converged(list(code = f$convergence), method))
    }, freq = w_panel)
    V <- bs$vcov
    fit$boot_failures <- bs$failures
    fit$n_boot_ok <- bs$n_ok
    fit$successful_bootstraps <- bs$n_ok
    fit$bootstraps <- as.integer(bootstraps)
  } else {
    V <- fit$vcov
  }
  dimnames(V) <- list(x_names, x_names)
  fit$vcov <- V
  fit$se <- .fcr_se(V)
  fit$bootstrapped_se <- if (stderr == "boot") fit$se else NULL
  fit$vcov_type <- switch(stderr,
                          normal = if (has_smooth) "bayesian" else "hessian",
                          robust = "robust",
                          boot = "bootstrap")
  fit$cre <- cre_info
  
  obj <- .createFlexCountReg(
    model = fit,
    data = data, keep_data = keep_data,
    call = call,
    formula = formula
  )
  
  obj
}


.repois_start_values <- function(y, X, offset_vec, w_obs = NULL) {
  fit0 <- try(
    stats::glm.fit(
      x = X,
      y = y,
      weights = w_obs,
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

.repois_loglik <- function(par, y, X, offset_vec, cluster_list, w = NULL) {
  # One log-likelihood value per cluster (maxLik sums them); per-cluster
  # values also give cluster-robust sandwich standard errors.
  n_cl <- length(cluster_list)
  p <- ncol(X)
  beta <- par[seq_len(p)]
  lnalpha <- par[p + 1L]
  alpha <- exp(lnalpha)
  
  if (!is.finite(alpha) || alpha <= 0) {
    return(rep(-Inf, n_cl))
  }
  
  eta <- drop(X %*% beta) + offset_vec
  if (any(!is.finite(eta)) || any(eta > 700)) {
    return(rep(-Inf, n_cl))
  }
  
  ll <- vapply(cluster_list, function(idx) {
    y_i <- y[idx]
    eta_i <- eta[idx]
    ysum <- sum(y_i)
    s <- sum(exp(eta_i))
    r <- 1 / alpha
    # lgamma(Y + r) - lgamma(r) + Y log(alpha) without cancellation, which
    # tends to the Poisson log-likelihood as alpha -> 0
    .fcr_lgamma_diff_rel(r, ysum) - (ysum + r) * log1p(alpha * s) +
      sum(y_i * eta_i) - sum(lgamma(y_i + 1))
  }, numeric(1))
  ll <- .fcr_weight_ll(unname(ll), w)
  if (any(!is.finite(ll))) return(rep(-Inf, n_cl))
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
                               compute_vcov = TRUE,
                               w = NULL,
                               smooth = NULL) {
  if (!requireNamespace("maxLik", quietly = TRUE)) {
    stop("Package `maxLik` is required for `repois()`.")
  }
  
  method_opt <- .fcr_check_method(method)
  # Row weights for the starting values: each row has its panel's weight
  w_obs <- if (!is.null(w)) rep(w, lengths(cluster_list))[
    order(unlist(cluster_list, use.names = FALSE))]
  start <- .repois_start_values(y, X, offset_vec, w_obs)
  
  control <- maxLik::maxControl(
    iterlim = max.iters,
    printLevel = print.level
  )
  
  if (!is.null(smooth)) {
    names(start) <- c(colnames(X), "lnalpha")
    fit_ml <- .fcr_smooth_fit(
      smooth$ents,
      function(p) .repois_loglik(p, y = y, X = X, offset_vec = offset_vec,
                                 cluster_list = cluster_list, w = w),
      function(p) .repois_grad(p, y = y, X = X, offset_vec = offset_vec,
                               cluster_list = cluster_list, w = w),
      start, method_opt, max.iters, smooth$sp, list(main = 0),
      list(main = "mean"), list(main = X))
    method_opt <- fit_ml$method
  } else {
    fit_ml <- maxLik::maxLik(
      logLik = .repois_loglik,
      start = start,
      method = method_opt,
      control = control,
      y = y,
      X = X,
      offset_vec = offset_vec,
      cluster_list = cluster_list,
      w = w
    )
  }
  
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
    
    if (all(is.na(vcov_mat))) {
      vc_try <- try(stats::vcov(fit_ml), silent = TRUE)
      if (!inherits(vc_try, "try-error") && all(is.finite(vc_try))) {
        vcov_mat <- vc_try
      }
    }
    
    if (all(is.finite(vcov_mat))) {
      se <- .fcr_se(vcov_mat)
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
    bic = -2 * logLik_value +
      log(if (is.null(w)) n_groups else sum(w)) * k,
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
  if (!is.null(smooth)) {
    keep <- c("smooth_table", "smooth_info", "smooth_lambda",
              "smooth_converged", "smooth_iterations", "smooth_coefs", "edf",
              "df_total", "gradientObs")
    out[keep] <- fit_ml[keep]
  }
  
  if (isTRUE(compute_vcov) && all(is.finite(vcov_mat))) {
    out$se <- .fcr_se(vcov_mat)
  }
  
  
  out
}




`%||%` <- function(a, b) if (!is.null(a)) a else b

# Analytic per-cluster scores of .repois_loglik() (rows: clusters, in the
# order of cluster_list; columns: beta, ln(alpha)). With r = 1 / alpha,
# S = sum(mu) and Y = sum(y) over a cluster,
#   d/d beta      = sum(y x) - (Y + r) alpha / (1 + alpha S) sum(mu x)
#   d/d ln(alpha) = Y - r [psi(Y + r) - psi(r)] + r log1p(alpha S)
#                   - (Y + r) alpha S / (1 + alpha S).
.repois_grad <- function(par, y, X, offset_vec, cluster_list, w = NULL) {
  p <- ncol(X)
  alpha <- exp(par[p + 1L])
  r <- 1 / alpha
  g <- integer(length(y))
  for (k in seq_along(cluster_list)) g[cluster_list[[k]]] <- k
  mu <- exp(drop(X %*% par[seq_len(p)]) + offset_vec)
  S <- as.numeric(rowsum(mu, g, reorder = TRUE))
  Y <- as.numeric(rowsum(y, g, reorder = TRUE))
  f <- (Y + r) * alpha / (1 + alpha * S)
  d_b <- rowsum(y * X, g, reorder = TRUE) - f * rowsum(mu * X, g,
                                                       reorder = TRUE)
  d_a <- Y - r * .fcr_digamma_diff(r, Y) + r * log1p(alpha * S) -
    (Y + r) * alpha * S / (1 + alpha * S)
  G <- unname(cbind(d_b, d_a))
  if (!is.null(w)) G <- G * w
  G
}
