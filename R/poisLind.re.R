#' Function for estimating a Random Effects Poisson-Lindley regression model
#'
#' @name poisLindRE
#' @param formula an R formula.
#' @param group_var the grouping variable(s) indicating random effects
#'   (e.g., individual ID).
#' @param data a dataframe that has all of the variables in the
#'   \code{formula}.
#' @param offset an optional offset term provided as a string.
#' @param method a method to use for optimization in the maximum likelihood
#'   estimation. For options, see \code{\link[maxLik]{maxLik}}. Note that
#'   "BHHH" is not available for this function due to the implementation
#'   for the random effects.
#' @param max.iters the maximum number of iterations to allow the
#'   optimization method to perform.
#' @param print.level Integer specifying the verbosity of output during
#'   optimization.
#' @param bootstraps Optional integer specifying the number of bootstrap
#'   samples to be used for estimating standard errors. If not specified,
#'   no bootstrapping is performed.
#'
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom maxLik maxLik
#' @importFrom modelr model_matrix
#' @importFrom MASS glm.nb
#' @importFrom purrr map map_df
#' @importFrom broom tidy
#' @importFrom dplyr group_by %>% across select mutate all_of reframe
#' @importFrom tibble deframe
#' @include plind.R
#'
#' @details
#' The function \code{poisLindRE} is similar to the \code{poisLind}
#' function, but it includes an additional argument \code{group_var} that
#' specifies the grouping variable for the random effects. The function
#' estimates a Random Effects Poisson-Lindley regression model using
#' maximum likelihood. It is similar to \code{poisLind}, but includes
#' additional terms to account for the random effects.
#'
#' The Random Effects Poisson-Lindley model is useful for panel data and
#' assumes that the random effects follow a gamma distribution. The PDF is
#' \deqn{
#' f(y_{it}|\mu_{it},\theta)=\frac{\theta^2}{\theta+1}
#' \prod_{t=1}^{n_i}\frac{\left(\mu_{it}\frac{\theta(\theta+1)}
#' {\theta+2}\right)^{y_{it}}}{y_{it}!}
#' \cdot
#' \frac{
#' \left(\sum_{t=1}^{n_i}y_{it}\right)!
#' \left(\sum_{t=1}^{n_i}\mu_{it}\frac{\theta(\theta+1)}{\theta+2}
#'       + \theta + \sum_{t=1}^{n_i}y_{it} + 1\right)
#' }{
#' \left(\sum_{t=1}^{n_i}\mu_{it}\frac{\theta(\theta+1)}{\theta+2}
#'       + \theta\right)^{\sum_{t=1}^{n_i}y_{it}+2}
#' }
#' }
#'
#' The log-likelihood function is:
#' \deqn{
#' LL = 2\log(\theta) - \log(\theta+1)
#'   + \sum_{t=1}^{n_i} y_{it}\log(\mu_{it})
#'   + \sum_{t=1}^{n_i} y_{it}\log\!\left(
#'       \frac{\theta(\theta+1)}{\theta+2}
#'     \right)
#'   - \sum_{t=1}^{n_i}\log(y_{it}!)
#'   + \log\!\left(
#'       \left(\sum_{t=1}^{n_i}y_{it}\right)!
#'     \right)
#'   + \log\!\left(
#'       \sum_{t=1}^{n_i}\mu_{it}\frac{\theta(\theta+1)}{\theta+2}
#'       + \theta + \sum_{t=1}^{n_i}y_{it} + 1
#'     \right)
#'   - \left(\sum_{t=1}^{n_i}y_{it} + 2\right)
#'     \log\!\left(
#'       \sum_{t=1}^{n_i}\mu_{it}\frac{\theta(\theta+1)}{\theta+2}
#'       + \theta
#'     \right)
#' }
#'
#' The mean and variance are:
#' \deqn{\mu_{it}=\exp(X_{it} \beta)}
#' \deqn{
#' V(\mu_{it})=\mu_{it}+
#' \left(1-\frac{2}{(\theta+2)^2}\right)\mu_{it}^2
#' }
#' 
#' @returns
#' An object of class `countreg` which is a list with the following components:
#' \itemize{
#'  \item model: the fitted model object.
#'  \item data: the data frame used to fit the model.
#'  \item call: the matched call.
#'  \item formula: the formula used to fit the model.
#' }
#'
#' @examples
#' \donttest{
#' data("washington_roads")
#' washington_roads$AADTover10k <-
#'   ifelse(washington_roads$AADT > 10000, 1, 0)
#'
#' poislind.mod <- poisLind.re(
#'   Animal ~ lnaadt + lnlength + speed50 +
#'     ShouldWidth04 + AADTover10k,
#'   data      = washington_roads,
#'   group_var = "ID",
#'   method    = "NM",
#'   max.iters = 1000
#' )
#' summary(poislind.mod)
#' }
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @export
poisLind.re <- function(formula, group_var, data,
                        method = "NM", max.iters = 1000,
                        print.level = 0, bootstraps = NULL,
                        offset = NULL) {
  
  # if (method == "BHHH" | method == "bhhh") {
  #   msg <- paste(
  #     "The `BHHH` method is not available for poisLindRE.",
  #     "Switching to `NM` method."
  #   )
  #   print(msg)
  #   method <- "NM"
  # }
  
  # Data preparation
  mod_df <- stats::model.frame(formula, data, na.action = stats::na.fail)
  X <- stats::model.matrix(attr(mod_df, "terms"), mod_df)
  y <- as.numeric(stats::model.response(mod_df))
  
  # Generate a panel ID for the model
  if (!("panel_id" %in% names(data))) {
    if (is.null(group_var)) {
      warning("The `group_var` must be defined for this model.")
    } else {
      if (length(group_var) > 1) {
        # Combine multiple columns into a single string column
        data <- data %>%
          tidyr::unite("panel_id", all_of(group_var),
                sep = "_", remove = FALSE)
      } else {
        data <- data %>%
          mutate(panel_id = as.character(.data[[group_var]]))
      }
    }
  }
  
  group <- data$panel_id
  x_names <- colnames(X)
  
  # Negative Binomial as starting values
  p_model <- glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  t <- 116761 / exp(-9.04 * p_model$theta)
  theta <- ifelse(t < 100, t, 1)  # Approximate initial theta
  
  full_start <- append(start, log(theta))
  x_names <- append(x_names, "ln(theta)")
  names(full_start) <- x_names
  
  # Local helpers keep this replacement independent of other repair files.
  .fcr_lse <- function(x) {
    if (!length(x)) return(-Inf)
    if (anyNA(x)) return(NA_real_)
    m <- max(x)
    if (is.infinite(m)) return(m)
    m + log(sum(exp(x - m)))
  }

  .fcr_offset <- function(mf, data, offset_cols = NULL) {
    form_off <- stats::model.offset(mf)
    if (!is.null(form_off) && !is.null(offset_cols))
      stop("Use either a formula offset or offset columns")
    if (!is.null(form_off)) out <- as.numeric(form_off)
    else if (is.null(offset_cols)) out <- rep(0, nrow(data))
    else {
      if (!is.character(offset_cols) || !length(offset_cols) ||
          anyDuplicated(offset_cols) ||
          !all(offset_cols %in% names(data)))
        stop("offset must name distinct columns in data")
      z <- data[, offset_cols, drop = FALSE]
      if (!all(vapply(z, is.numeric, logical(1))))
        stop("Offset columns must be numeric")
      out <- rowSums(z)
    }
    if (length(out) != nrow(data) || any(!is.finite(out)))
      stop("Offsets must be finite and aligned with model rows")
    out
  }

  offset_vec <- .fcr_offset(mod_df, data, offset)

  # Poisson-Lindley log likelihood, evaluated separately for each panel.
  reg.run.RE <- function(beta, y, X, group, offset_vec) {
    n_covariates <- ncol(X)
    if (length(beta) != n_covariates + 1L || nrow(X) != length(y) ||
        length(group) != length(y) || length(offset_vec) != length(y) ||
        anyNA(group)) {
      stop("Poisson-Lindley parameters, rows, groups and offsets must align.")
    }
    beta_cov <- beta[seq_len(n_covariates)]
    theta <- exp(beta[n_covariates + 1L])
    panels <- split(seq_along(y), group, drop = TRUE)
    n_panels <- length(panels)

    if (!is.finite(theta) || theta <= 0) {
      return(rep(-Inf, n_panels))
    }

    # Convert the marginal mean into the conditional Poisson rate multiplier.
    log_theta <- log(theta)
    log_1p_theta <- log1p(theta)
    log_lambda <- drop(X %*% beta_cov) + offset_vec +
      log_theta + log_1p_theta - log(theta + 2)
    if (any(!is.finite(log_lambda))) {
      return(rep(-Inf, n_panels))
    }
    theta_const <- 2 * log_theta - log_1p_theta

    # Integrate over the shared Lindley random effect within each panel.
    vapply(panels, function(idx) {
      y_panel <- y[idx]
      log_lambda_panel <- log_lambda[idx]
      y_sum <- sum(y_panel)
      log_A <- .fcr_lse(log_lambda_panel)
      log_denom <- .fcr_lse(c(log_A, log_theta))
      log_linear <- .fcr_lse(c(log_A, log(theta + y_sum + 1)))
      log_lik_obs <- sum(y_panel * log_lambda_panel - lgamma(y_panel + 1))
      log_lik_gamma <- lgamma(y_sum + 1)

      theta_const + log_lik_obs + log_lik_gamma +
        log_linear - (y_sum + 2) * log_denom
    }, numeric(1))
  }

  fit <- maxLik::maxLik(
    reg.run.RE,
    start  = full_start,
    y      = y,
    X      = X,
    group  = group,
    offset_vec = offset_vec,
    method = method,
    control = list(
      iterlim    = max.iters,
      printLevel = print.level
    )
  )
  
  # Bootstrap SEs (optional)
  plind.boot <- function(data) {
    data <- as.data.frame(data)
    mod1_frame <- stats::model.frame(formula, data, na.action = stats::na.fail)
    X_Fixed <- stats::model.matrix(attr(mod1_frame, "terms"), mod1_frame)
    y <- stats::model.response(mod1_frame)
    
    int_res <- maxLik::maxLik(
      reg.run.RE,
      start  = fit$estimate,
      y      = y,
      X      = X_Fixed,
      group  = data$panel_id,
      offset_vec = .fcr_offset(mod1_frame, data, offset),
      method = method,
      control = list(
        iterlim    = max.iters,
        printLevel = print.level
      )
    )
    
    return(int_res)
  }
  
  if (!is.null(bootstraps) && is.numeric(bootstraps)) {
    bs.data <- modelr::bootstrap(data, n = bootstraps)
    
    mod1_frame <- stats::model.frame(formula, data)
    X_Fixed <- stats::model.matrix(formula, data)
    y <- stats::model.response(mod1_frame)
    
    models <- map(bs.data$strap, ~ plind.boot(data = .x))
    tidied <- map_df(models, broom::tidy, .id = "id")
    
    SE <- tidied %>%
      group_by(term) %>%
      reframe(sd = sd(estimate))
    
    fit$bootstrapped_se <- SE
  }
  
  fit$bootstraps <- if (!is.null(bootstraps)) bootstraps else NULL
  
  beta_est <- fit$estimate
  npars <- length(beta_est) - 1
  beta_pred <- as.vector(unlist(beta_est[1:npars]))
  fit$beta_pred <- beta_pred
  fit$theta <- exp(fit$estimate[length(fit$estimate)])
  
  mu <- exp(drop(X %*% beta_pred) + offset_vec)
  fit$predictions <- mu
  fit$se <- sqrt(diag(vcov(fit)))
  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum
  fit$modelType <- "poisLindRE"
  fit$offset <- offset
  fit$offset_values <- offset_vec
  
  obj <- .createFlexCountReg(
    model = fit,
    data = data,
    call = match.call(),
    formula = formula
  )
  
  return(obj)
}
