#' Estimate a Random Effects Negative Binomial regression model
#'
#' @name renb
#' @param formula an R formula.
#' @param group_var the grouping variable(s) for the random effects (e.g.,
#'   individual ID or other panel ID variables).
#' @param data a dataframe that has all of the variables in the \code{formula}.
#' @param offset an optional offset term provided as a string.
#' @param method a method to use for optimization in the maximum likelihood
#'   estimation. For options, see \code{\link[maxLik]{maxLik}}. Note that "BHHH"
#'   is not available for this function due to the implementation for the random
#'   effects.
#' @param max.iters the maximum number of iterations to allow the optimization
#'   method to perform.
#' @param print.level Integer specifying the verbosity of output during
#'   optimization.
#' @param bootstraps Optional integer specifying the number of bootstrap samples
#'   to be used for estimating standard errors. If not specified, no
#'   bootstrapping is performed.
#'        
#' @import maxLik  stats modelr tibble
#' @importFrom MASS glm.nb
#' @importFrom purrr map map_df
#' @importFrom broom tidy
#' @importFrom dplyr group_by %>% across select mutate all_of reframe
#' @include renbLL.R
#' 
#' @details
#' This function estimates a random effects negative binomial (RENB) regression 
#' model. This model is based on the NB-1 model. The PDF for the RENB is:
#' \deqn{f(y_{it}|\lambda_{it}, a, b) = 
#'   \frac{\Gamma(a+b) 
#'     \Gamma(a + \sum_{t = 1}^{n_i} \\lambda_{it}) 
#'     \Gamma(b + \sum_{t=1}^{n_i}y_{it})}
#'     {\Gamma(a) \Gamma(b) \Gamma(a + b + 
#'        \sum_{t=1}^{n_i}\lambda_{it} + \sum_{t=1}^{n_i}y_{it})} \prod_{t=1}^{n_i}
#'        \frac{\Gamma(\lambda_{it}+y_{it})}{\Gamma(\lambda_{it})\Gamma(y_{it})}}
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
#' for count data with an application to the patents–R&D relationship." 
#' Econometrica: Journal of the Econometric Society (1984): 909-938.
#' 
#' @export
renb <- function(formula, group_var, data, method = 'NM', max.iters = 1000, 
                 print.level=0, bootstraps=NULL, offset=NULL) {
  
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
        data <- data %>% 
          unite("panel_id", all_of(group_var), sep = "_", remove = FALSE)
      } else {
        data <- data %>% mutate(panel_id = as.character(data[[group_var]]))
      }
    }
  }
  
  group <- data$panel_id
  x_names <- colnames(X)
  
  # Use the Negative Binomial as starting values
  p_model <- glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  a <- 2
  b <- 1
  
  full_start <- append(start, log(a - 1))
  x_names <- append(x_names, "ln(a-1)")
  full_start <- append(full_start, log(b))
  x_names <- append(x_names, "ln(b)")
  names(full_start) <- x_names


  # Extract offsets on the log scale from the same rows used by the model.
  get_offset <- function(model_frame, model_data, offset_cols) {
    formula_offset <- stats::model.offset(model_frame)
    if (!is.null(formula_offset) && !is.null(offset_cols)) {
      stop("Use either a formula offset or offset columns, not both.")
    }
    if (!is.null(formula_offset)) {
      out <- as.numeric(formula_offset)
    } else if (is.null(offset_cols)) {
      out <- rep(0, nrow(model_frame))
    } else {
      if (!is.character(offset_cols) || !length(offset_cols) ||
          anyDuplicated(offset_cols) ||
          !all(offset_cols %in% names(model_data))) {
        stop("offset must name distinct columns in data.")
      }
      columns <- as.data.frame(model_data)[, offset_cols, drop = FALSE]
      if (!all(vapply(columns, is.numeric, logical(1)))) {
        stop("Offset columns must be numeric and already on the log scale.")
      }
      out <- rowSums(columns)
    }
    if (length(out) != nrow(model_frame) || any(!is.finite(out))) {
      stop("Offsets must be finite and aligned with the model rows.")
    }
    out
  }
  offset_vec <- get_offset(mod_df, data, offset)

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

  # Bootstrap function
  plind.boot <- function(boot_data, formula, method, 
                         max.iters, print.level, offset) {
    # Prepare bootstrapped data
    boot_data <- as.data.frame(boot_data)
    mod1_frame <- stats::model.frame(formula, boot_data,
                                     na.action = stats::na.fail)
    X_boot <- stats::model.matrix(attr(mod1_frame, "terms"), mod1_frame)
    y_boot <- as.numeric(stats::model.response(mod1_frame))
    group_boot <- boot_data$panel_id
    offset_boot <- get_offset(mod1_frame, boot_data, offset)
    
    # Fit model to bootstrapped data
    int_res <- try(maxLik::maxLik(reg.run.RE,  
                                  start = fit$estimate,
                                  y = y_boot,
                                  X = X_boot,
                                  group = group_boot,
                                  offset_vec = offset_boot,
                                  method = method,
                                  control = list(iterlim = max.iters, 
                                                 printLevel = print.level)),
                   silent = TRUE)
    
    # Return NULL if the bootstrap fit failed
    if(inherits(int_res, "try-error")) return(NULL)
    
    return(int_res)
  }
  
  # Perform bootstrapping if requested - Modified bootstrap implementation
  if (!is.null(bootstraps) && is.numeric(bootstraps)) {
    # Existing row bootstrap; panel resampling requires a separate change.
    bs.data <- modelr::bootstrap(data, n = bootstraps)
    
    # Run bootstrap models with correct parameter passing
    models <- map(bs.data$strap, ~plind.boot(
      boot_data = .,
      formula = formula,
      method = method,
      max.iters = max.iters,
      print.level = print.level,
      offset = offset
    ))
    
    # Remove failed bootstrap iterations
    models <- models[!vapply(models, is.null, logical(1))]
    
    # Calculate bootstrap standard errors
    if(length(models) > 0) {
      tidied <- map_df(models, ~{
        if(!is.null(.)) {
          data.frame(
            term = names(.x$estimate),
            estimate = as.numeric(.x$estimate)
          )
        }
      }, .id = "id")
      
      SE <- tidied %>%
        group_by(term) %>%
        reframe(sd = sd(estimate))
      
      fit$bootstrapped_se <- SE
      fit$successful_bootstraps <- length(models)
    } else {
      msg <- paste(
        "All bootstrap iterations failed.", 
        "No bootstrap standard errors computed.")
      warning(msg)
      fit$bootstrapped_se <- NULL
      fit$successful_bootstraps <- 0
    }
  }
  
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
  fit$se <- sqrt(diag(vcov(fit)))
  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum
  fit$modelType <- "RENB"
  fit$offset <- offset
  fit$offset_values <- offset_vec
  
  obj <- .createFlexCountReg(model = fit, 
                             data = data, 
                             call = match.call(), 
                             formula = formula)
  return(obj)
}
