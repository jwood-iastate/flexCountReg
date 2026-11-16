#' Function for estimating a Random Effects Poisson-Lindley regression model
#'
#' @name poisLindRE
#' @aliases poislind_re poisLind.re
#' @param formula an R formula.
#' @param group_var the grouping variable(s) indicating random effects
#'   (e.g., individual ID). Required. Several columns are combined into one
#'   panel identifier.
#' @param data a dataframe that has all of the variables in the
#'   \code{formula}.
#' @param offset an optional character vector of offset column names (summed
#'   when several are given; on the log scale). An offset can instead be given
#'   with \code{offset()} in \code{formula}; use only one of the two.
#' @param method a method to use for optimization in the maximum likelihood
#'   estimation. For options, see \code{\link[maxLik]{maxLik}}. Note that
#'   "BHHH" is not available for this function due to the implementation
#'   for the random effects.
#' @param max.iters the maximum number of iterations to allow the
#'   optimization method to perform.
#' @param verbose logical; if \code{TRUE}, print the progress of the
#'   optimization.
#' @param stderr the type of standard errors: \code{"normal"} (Hessian-based,
#'   the default), \code{"robust"} (sandwich; clustered by panel because the
#'   likelihood has one contribution per panel), or \code{"boot"} (bootstrap
#'   that resamples whole panels). Giving \code{bootstraps} without
#'   \code{stderr} also requests bootstrapped standard errors.
#' @param keep_data logical; if \code{FALSE}, the training data are not stored
#'   in the fitted object (saving memory). \code{predict()} then needs
#'   \code{newdata}, and \code{fitted()}/\code{residuals()} are unavailable.
#' @param print.level deprecated; use \code{verbose}.
#' @param bootstraps Optional integer specifying the number of bootstrap
#'   samples to be used for estimating standard errors. Whole panels
#'   (groups) are resampled, so the standard errors account for the
#'   correlation within panels. If not specified, no bootstrapping is
#'   performed.
#'
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom maxLik maxLik
#' @importFrom MASS glm.nb
#' @importFrom tibble deframe
#' @include plind.R
#'
#' @details
#' \code{poislind_re()} is the same function as \code{poisLind.re()}, under a
#' name that tools do not mistake for an S3 method.
#'
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
#' An object of class `flexCountReg` which is a list with the following components:
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
                        verbose = FALSE, stderr = "normal",
                        bootstraps = NULL, offset = NULL,
                        keep_data = TRUE, print.level = NULL) {
  print.level <- .fcr_print_level("poisLind.re", verbose, print.level)
  stderr <- .fcr_resolve_stderr(stderr, missing(stderr), bootstraps)
  
  method <- .fcr_check_method(method)
  if (missing(group_var) || is.null(group_var)) {
    stop("`group_var` must name the column(s) that identify the panels.",
         call. = FALSE)
  }
  
  # Data preparation: one complete-case frame; panels from `group_var`
  prep <- prepare_model_data(formula, data, offset = offset,
                             group_var = group_var)
  formula <- prep$formula   # with `.` expanded
  if (prep$n_dropped > 0) {
    message(prep$n_dropped, " row(s) with missing values were removed.")
  }
  X <- prep$X
  y <- as.numeric(prep$y)
  group <- prep$group
  offset_vec <- prep$offset
  x_names <- c(colnames(X), "ln(theta)")
  
  # Starting values: Poisson GLM with the same rows and offset, and a moment
  # estimate of theta. For the mean-parameterized Poisson-Lindley,
  # Var(Y) = mu + phi mu^2 with phi = 1 - 2 / (theta + 2)^2 in (1/2, 1).
  p_model <- stats::glm.fit(X, y, offset = offset_vec,
                            family = stats::poisson())
  start <- p_model$coefficients
  start[is.na(start)] <- 0
  mu0 <- p_model$fitted.values
  phi <- sum((y - mu0)^2 - y) / sum(mu0^2)
  phi <- min(max(phi, 0.51), 0.99)
  theta <- sqrt(2 / (1 - phi)) - 2
  full_start <- c(start, log(theta))
  names(full_start) <- x_names
  
  .fcr_lse <- function(x) {
    if (!length(x)) return(-Inf)
    if (anyNA(x)) return(NA_real_)
    m <- max(x)
    if (is.infinite(m)) return(m)
    m + log(sum(exp(x - m)))
  }

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
  
  fit$fcr_method <- method
  fit$method <- method
  fit <- .fcr_finish_re(fit, prep, x_names, stderr, bootstraps,
                        refit = function(sub) {
    m <- maxLik::maxLik(reg.run.RE, start = fit$estimate, y = sub$y,
                        X = sub$X, group = sub$group,
                        offset_vec = sub$offset, method = method,
                        control = list(iterlim = max.iters, printLevel = 0))
    stats::setNames(m$estimate, x_names)
  })
  
  beta_est <- fit$estimate
  npars <- length(beta_est) - 1
  beta_pred <- as.vector(unlist(beta_est[1:npars]))
  fit$beta_pred <- beta_pred
  fit$theta <- exp(fit$estimate[length(fit$estimate)])
  
  mu <- exp(drop(X %*% beta_pred) + offset_vec)
  fit$predictions <- mu
  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum
  fit$modelType <- "poisLindRE"
  fit$offset <- prep$blueprint$offset_label
  fit$offset_values <- offset_vec
  
  obj <- .createFlexCountReg(
    model = fit,
    data = data, keep_data = keep_data,
    call = match.call(),
    formula = formula
  )
  
  return(obj)
}

#' @rdname poisLindRE
#' @usage NULL
#' @export
poislind_re <- poisLind.re
