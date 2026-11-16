#' Estimate a Random Effects Negative Binomial regression model
#'
#' @name renb
#' @param formula an R formula.
#' @param group_var the grouping variable(s) for the random effects (e.g.,
#'   individual ID or other panel ID variables). Required. Several columns
#'   are combined into one panel identifier.
#' @param data a dataframe that has all of the variables in the \code{formula}.
#' @param offset an optional character vector of offset column names (summed
#'   when several are given; on the log scale). An offset can instead be given
#'   with \code{offset()} in \code{formula}; use only one of the two.
#' @param method a method to use for optimization in the maximum likelihood
#'   estimation. For options, see \code{\link[maxLik]{maxLik}}. Note that "BHHH"
#'   is not available for this function due to the implementation for the random
#'   effects.
#' @param max.iters the maximum number of iterations to allow the optimization
#'   method to perform.
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
#' @param bootstraps Optional integer specifying the number of bootstrap samples
#'   to be used for estimating standard errors. Whole panels (groups) are
#'   resampled, so the standard errors account for the correlation within
#'   panels. If not specified, no bootstrapping is performed.
#'        
#' @importFrom maxLik maxLik
#' @importFrom MASS glm.nb
#' @include renbLL.R
#' 
#' @details
#' This function estimates a random effects negative binomial (RENB) regression 
#' model. This model is based on the NB-1 model. The PDF for the RENB is:
#' \deqn{f(y_{it}|\lambda_{it}, a, b) = 
#'   \frac{\Gamma(a+b) 
#'     \Gamma(a + \sum_{t = 1}^{n_i} \lambda_{it}) 
#'     \Gamma(b + \sum_{t=1}^{n_i}y_{it})}
#'     {\Gamma(a) \Gamma(b) \Gamma(a + b + 
#'        \sum_{t=1}^{n_i}\lambda_{it} + \sum_{t=1}^{n_i}y_{it})} \prod_{t=1}^{n_i}
#'        \frac{\Gamma(\lambda_{it}+y_{it})}{\Gamma(\lambda_{it})\Gamma(y_{it}+1)}}
#'        
#' Where \eqn{y_{it}} is the count outcome for individual \eqn{i} at time 
#' \eqn{t}, and \eqn{\lambda_{it}} is the latent Poisson mean parameter for 
#' individual \eqn{i} at time \eqn{t}. The parameters \eqn{a} and \eqn{b} are 
#' the shape parameters for the beta distribution that is used to model the 
#' random effects. The RENB model allows for overdispersion in the count data 
#' and accounts for unobserved heterogeneity across individuals by including 
#' random effects in the model. This formulation follows the approach described 
#' in the paper by Hausman, Hall, and Griliches (1984) for modeling panel data
#' with random effects.
#' 
#' The marginal mean and marginal variance of the RENB model are given by:
#' \deqn{E[y_{it}] = \lambda_{it}\frac{b}{a-1}=\mu_it=exp(X_{it}\beta)}
#' 
#' \deqn{Var[y_{it}] = \frac{a+b-1}{a-2}\mu_{it} + \frac{a+b-1}{b(a-2)} \mu_{it}^2}
#' 
#' Thus, the formulation of the model estimated here allows the use of the 
#' estimated coefficients to directly compute the marginal mean.
#' 
#' Note that the RENB model is a panel data model, and the \code{group_var} 
#' argument must be specified to indicate the grouping variable(s) for the 
#' random effects. The model is estimated using maximum likelihood estimation, 
#' and the optimization is performed using the \code{\link[maxLik]{maxLik}} 
#' package. The user can specify the optimization method and maximum iterations.
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
#' ## RENB Model
#' data("washington_roads")
#' washington_roads$AADTover10k <- 
#'   ifelse(washington_roads$AADT > 10000, 1, 0) # create a dummy variable
#' renb.mod <- renb(Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 offset = "lnlength",
#'                                 group_var="ID",
#'                                 method="nm",
#'                                 max.iters = 1000)
#' summary(renb.mod)
#' }
#' 
#' @references
#' Hausman, Jerry A., Bronwyn H. Hall, and Zvi Griliches. "Econometric models 
#' for count data with an application to the patents-R&D relationship." 
#' Econometrica: Journal of the Econometric Society (1984): 909-938.
#' 
#' @export
renb <- function(formula, group_var, data, method = 'NM', max.iters = 1000, 
                 verbose = FALSE, stderr = "normal", bootstraps = NULL,
                 offset = NULL, keep_data = TRUE, print.level = NULL) {
  print.level <- .fcr_print_level("renb", verbose, print.level)
  stderr <- .fcr_resolve_stderr(stderr, missing(stderr), bootstraps)
  
  method <- .fcr_check_method(method)
  if (missing(group_var) || is.null(group_var)) {
    stop("`group_var` must name the column(s) that identify the panels.",
         call. = FALSE)
  }
  
  # Data preparation: one complete-case frame; panels from `group_var`
  prep <- prepare_model_data(formula, data, offset = offset,
                             group_var = group_var)
  if (prep$n_dropped > 0) {
    message(prep$n_dropped, " row(s) with missing values were removed.")
  }
  X <- prep$X
  y <- as.numeric(prep$y)
  group <- prep$group
  offset_vec <- prep$offset
  x_names <- c(colnames(X), "ln(a-1)", "ln(b)")
  
  # Starting values: Poisson GLM with the same rows and offset
  p_model <- stats::glm.fit(X, y, offset = offset_vec,
                            family = stats::poisson())
  start <- p_model$coefficients
  start[is.na(start)] <- 0
  a <- 2
  b <- 1
  full_start <- c(start, log(a - 1), log(b))
  names(full_start) <- x_names

  # Random Effects Negative Binomial log-likelihood, one value per panel.
  reg.run.RE <- function(beta, y, X, group, offset_vec) {
    n_covariates <- ncol(X)
    n_panels <- length(unique(group))
    if (length(beta) != n_covariates + 2L ||
        length(offset_vec) != length(y) || nrow(X) != length(y) ||
        length(group) != length(y) || anyNA(group)) {
      stop("RENB parameters, responses, groups and offsets must align.")
    }

    # The marginal-mean parameterization requires a > 1 and b > 0.
    beta_cov <- beta[seq_len(n_covariates)]
    a <- 1 + exp(beta[n_covariates + 1L])
    b <- exp(beta[n_covariates + 2L])
    if (!is.finite(a) || a <= 1 || !is.finite(b) || b <= 0) {
      return(rep(-Inf, n_panels))
    }

    # lambda is the conditional NB shape used by renb_ll().
    log_mean <- drop(X %*% beta_cov) + offset_vec
    log_lambda <- log_mean + log(a - 1) - log(b)
    lambda <- exp(log_lambda)
    if (any(!is.finite(lambda)) || any(lambda <= 0)) {
      return(rep(-Inf, n_panels))
    }

    panel_loglik <- renb_ll(y = y, mu = lambda, a = a, b = b,
                           panels = group)
    if (any(!is.finite(panel_loglik))) {
      return(rep(-Inf, n_panels))
    }
    as.numeric(panel_loglik)
  }

  if (any(!is.finite(reg.run.RE(full_start, y, X, group, offset_vec)))) {
    stop("RENB log-likelihood is not finite at the starting values.")
  }

  # Main model fit
  fit <- maxLik::maxLik(reg.run.RE,
                        start = full_start,
                        y = y,
                        X = X,
                        group = group,
                        offset_vec = offset_vec,
                        method = method,
                        control = list(iterlim = max.iters, 
                                       printLevel = print.level))
  
  if (is.null(fit$estimate) || length(fit$estimate) != length(full_start) ||
      any(!is.finite(fit$estimate))) {
    stop("RENB optimizer did not return valid parameter estimates.")
  }

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
  
  # Process results
  beta_est <- fit$estimate
  npars <- length(beta_est)-2
  beta_pred <- as.vector(unlist(beta_est[1:npars]))
  fit$beta_pred <- beta_pred
  fit$a <- 1 + exp(
    unlist(fit$estimate[(length(fit$estimate)-1)])
  )
  fit$b <- exp(unlist(fit$estimate[length(fit$estimate)]))
  
  mu <- exp(drop(X %*% beta_pred) + offset_vec)
  fit$predictions <- mu
  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum
  fit$modelType <- "RENB"
  fit$offset <- prep$blueprint$offset_label
  fit$offset_values <- offset_vec
  
  obj <- .createFlexCountReg(model = fit, 
                             data = data, keep_data = keep_data, 
                             call = match.call(), 
                             formula = formula)
  return(obj)
}
