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
#' \deqn{f(y_{it}|\mu_{it}, a, b) = 
#'   \frac{\Gamma(a+b) + 
#'     \Gamma(a + \sum_{t = 1}^{n_i} \mu_{it}) + 
#'     \Gamma(b + \sum_{t=1}^{n_i}y_{it})}
#'     {\Gamma(a) \Gamma(b) \Gamma(a + b + 
#'        \sum_{t=1}^{n_i}\mu_{it} + \sum_{t=1}^{n_i}y_{it})} \prod_{t=1}^{n_i}
#'        \frac{\Gamma(\mu_{it}+y_{it})}{\Gamma(\mu_{it})\Gamma(y_{it})}}
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
#' @export
renb <- function(formula, group_var, data, method = 'NM', max.iters = 1000, 
                 print.level=0, bootstraps=NULL, offset=NULL) {
  
  # Data preparation
  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
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
  a <- 1
  b <- 1
  
  full_start <- append(start, log(a))
  x_names <- append(x_names, 'ln(a)')
  full_start <- append(full_start, log(b))
  x_names <- append(x_names, 'ln(b)')
  names(full_start) <- x_names
  
  # Log-likelihood function for the Random Effects Negative Binomial model
  reg.run.RE <- function(beta, y, X, group){
    pars <- length(beta)-2
    coefs <- as.vector(unlist(beta[1:pars]))
    a <- exp(unlist(beta[(pars+1)]))
    b <- exp(unlist(beta[(pars+2)]))
    
    mu <- exp(X %*% coefs)
    
    if(!is.null(offset)){ 
      if (length(offset)==1){
        mu <- mu * exp(data[[offset]])
      } else {
        for (i in offset){
          mu <- mu * exp(data[[i]])
        }
      }
    }
    
    mu <- as.vector(mu)
    LL <- renb_ll(y=y, mu=mu, a=a, b=b, panels=group)
    return(as.vector(LL))
  }
  
  # Main model fit
  fit <- maxLik::maxLik(reg.run.RE,
                        start = full_start,
                        y = y,
                        X = X,
                        group = group,
                        method = method,
                        control = list(iterlim = max.iters, 
                                       printLevel = print.level))
  
  # Bootstrap function - Modified to fix the error
  plind.boot <- function(boot_data, formula, method, 
                         max.iters, print.level, offset) {
    # Prepare bootstrapped data
    mod1_frame <- stats::model.frame(formula, boot_data)
    X_boot <- as.matrix(modelr::model_matrix(boot_data, formula))
    y_boot <- as.numeric(stats::model.response(mod1_frame))
    group_boot <- boot_data$panel_id
    
    # Fit model to bootstrapped data
    int_res <- try(maxLik::maxLik(reg.run.RE,  
                                  start = fit$estimate,
                                  y = y_boot,
                                  X = X_boot,
                                  group = group_boot,
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
    # Generate bootstrap samples preserving panel structure
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
  fit$a <- exp(unlist(fit$estimate[(length(fit$estimate)-1)]))
  fit$b <- exp(unlist(fit$estimate[length(fit$estimate)]))
  
  mu <- exp(X %*% beta_pred)
  fit$predictions <- mu
  fit$se <- sqrt(diag(-1/(fit$hessian)))
  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum
  fit$modelType <- "RENB"
  fit$offset <- offset
  
  obj <- .createFlexCountReg(model = fit, 
                             data = data, 
                             call = match.call(), 
                             formula = formula)
  return(obj)
}