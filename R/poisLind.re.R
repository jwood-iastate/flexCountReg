#' Function for estimating a Random Effects Poisson-Lindley regression model
#'
#' @name poisLindRE
#' @param formula an R formula.
#' @param group_var the grouping variable indicating random effects (e.g., individual ID).
#' @param data a dataframe that has all of the variables in the \code{formula}.
#' @param method a method to use for optimization in the maximum likelihood 
#' estimation. For options, see \code{\link[maxLik]{maxLik}}. Note that "BHHH"
#' is not available for this function due to the implementation for the random effects.
#' @param max.iters the maximum number of iterations to allow the optimization 
#' method to perform.
#' @param print.level Integer specifying the verbosity of output during optimization.
#' @param bootstraps Optional integer specifying the number of bootstrap samples to be used
#'        for estimating standard errors. If not specified, no bootstrapping is performed.
#' @param weights the name of the weighting variable. This is an optional 
#'   parameter for weighted regression.
#'
#' @import maxLik  stats modelr
#' @importFrom MASS glm.nb
#' @importFrom purrr map map_df
#' @importFrom broom tidy
#' @importFrom dplyr group_by %>% summarise
#' @importFrom tibble deframe
#' @include plind.R
#' 
#' @details
#' The function  poisLindRE  is similar to the  poisLind  function, but it 
#' includes an additional argument  group_var  that specifies the grouping 
#' variable for the random effects. The function estimates a Random Effects 
#' Poisson-Lindley regression model using the maximum likelihood method. 
#' The function is similar to the  poisLind  function, but it includes 
#' additional terms to account for the random effects. 
#' 
#' The PDF for the Random Effects Poisson-Lindley model is:
#' \deqn{f(y_{it}|\mu_{it},\theta)=\frac{\theta^2}{\theta+1} \prod_{t=1}^{n_i} \frac{\left(\mu_{it} \frac{\theta(\theta+1)}{\theta+2}\right)^{y_{it}}}{y_{it}!} \cdot \frac{\left(\sum_{t=1}^{n_i} y_{it}\right)! \left(\sum_{t=1}^{n_i} \mu_{it} \frac{\theta(\theta+1)}{\theta+2} + \theta + \sum_{t=1}^{n_i} y_{it} + 1\right)}{\left(\sum_{t=1}^{n_i} \mu_{it} \frac{\theta(\theta+1)}{\theta+2} + \theta\right)^{\sum_{t=1}^{n_i} y_{it} + 2}}}
#' 
#' The log-likelihood function for the Random Effects Poisson-Lindley model is:
#' \deqn{LL=2\log(\theta) - \log(\theta+1) + \sum_{t=1}^{n_i} y_{it} \log(\mu_{it}) + \sum_{t=1}^{n_i} y_{it} \log\left(\frac{\theta(\theta+1)}{\theta+2}\right) -\sum_{t=1}^{n_i} \log(y_{it}!) + \log\left(\left(\sum_{t=1}^{n_i} y_{it}\right)!\right) + \log\left(\sum_{t=1}^{n_i} \mu_{it} \frac{\theta(\theta+1)}{\theta+2} + \theta + \sum_{t=1}^{n_i} y_{it} + 1\right) - \left(\sum_{t=1}^{n_i} y_{it} + 2\right) \log\left(\sum_{t=1}^{n_i} \mu_{it} \frac{\theta(\theta+1)}{\theta+2} + \theta\right)}
#'  
#' The mean and variance are:
#' \deqn{\mu_{it}=\exp(X_{it} \beta)}
#' \deqn{V(\mu_{it})=\\mu_{it}+\left(1-\frac{2}{(\theta+2)^2}\right)\mu_{it}^2}
#' 
#' @examples
#' \donttest{
#'
#' ## Poisson-Lindley Random Effects Model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0) # create a dummy variable
#' poislind.mod <- poisLind.re(Animal ~ lnaadt + lnlength + speed50 +
#'                                 ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 group_var="ID",
#'                                 method="nm",
#'                                 max.iters = 1000)
#' summary(poislind.mod)}
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @export
poisLind.re <- function(formula, group_var, data, method = 'NM', max.iters = 1000, 
                       weights=NULL, print.level=0, bootstraps=NULL) {
  
  if (method=="BHHH" | method=="bhhh") {
    print("The `BHHH` method is not available for poisLindRE. Switching to `NM` method.")
    method="NM"
  }
  
  # Data preparation
  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  y <- as.numeric(stats::model.response(mod_df))
  group <- data[[group_var]]
  x_names <- colnames(X)

  # Preparing the weights
  if (!is.null(weights)){
    if (is.character(weights) && weights %in% colnames(data)) {
      # Use weights variable from data
      weights <- as.numeric(data[[weights]])
    } else {
      stop("Weights should be the name of a variable in the data.")
    }
  } else {
    weights <- rep(1, length(y)) # Default weights of 1
  }
  
  # Use the Negative Binomial as starting values
  p_model <- glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  t <- 116761/exp(-9.04*p_model$theta)
  theta <- ifelse(t<100,t,1) # Approximate Initial Theta
  
  full_start <- append(start, log(theta))
  x_names <- append(x_names, 'ln(theta)')
  names(full_start) <- x_names
  
  # Log-likelihood function for the Random Effects Poisson-Lindley model
  # reg.run.RE <- function(beta, y, X, group){
  #   pars <- length(beta)-1
  #   coefs <- as.vector(unlist(beta[1:pars]))
  #   theta <- exp(unlist(beta[length(beta)]))
  #   mu <- exp(X %*% coefs)
  #   
  #   # Compute the log-likelihood contributions per group
  #   unique_groups <- unique(group)
  #   
  #   probs <- sapply(unique_groups, function(g){
  #     y_group <- y[group == g]
  #     mu_group <- mu[group == g]
  #     ps <- pollind_i_group(mu_group, y_group, theta)
  #     return(ps)
  #   })
  #   
  #   group_weights <- sapply(unique_groups, function(g){
  #     sum(weights[group == g])
  #   })
  #   
  #   LL_groups <- log(probs)*group_weights
  #   
  #   
  #   # Sum log-likelihood contributions
  #   ll <- sum(LL_groups)
  #   return(ll)
  # }
  
  # Optimization using maxLik
  # Note: reg_run_RE is from Rcpp
  fit <- maxLik::maxLik(reg_run_RE,
                        start = full_start,
                        y = y,
                        X = X,
                        weights=weights,
                        group = group,
                        method = method,
                        control = list(iterlim = max.iters, 
                                       printLevel = print.level))
  
  # Optionally, compute bootstrapped standard errors
  # (Similar implementation as in poisLind function)
  plind.boot <- function(data){
    mod1_frame <- stats::model.frame(formula, data)
    X_Fixed <- stats::model.matrix(formula, data)
    y <- stats::model.response(mod1_frame)
    
    int_res <-  maxLik::maxLik(reg_run_RE,
                               start = fit$estimate,
                               y = y,
                               X = X,
                               weights=weights,
                               group = group,
                               method = method,
                               control = list(iterlim = max.iters, 
                                              printLevel = print.level))
    
    return(int_res)
  }
  
  
  if (!is.null(bootstraps) & is.numeric(bootstraps)) {
    bs.data <- modelr::bootstrap(data, n = bootstraps)
    
    
    mod1_frame <- stats::model.frame(formula, data)
    X_Fixed <- stats::model.matrix(formula, data)
    y <- stats::model.response(mod1_frame)
    
    models <- map(bs.data$strap, ~ plind.boot(data = .))
    tidied <- map_df(models, broom::tidy, .id = "id")
    
    SE <- tidied %>%
      group_by(term) %>%
      summarise(sd = sd(estimate)) %>% deframe()
    
    fit$bootstrapped_se <- SE
  }
  
  fit$bootstraps = if (!is.null(bootstraps)) bootstraps else NULL
  
  # Processing results
  beta_est <- fit$estimate
  npars <- length(beta_est)-1
  beta_pred <- as.vector(unlist(beta_est[1:npars]))
  fit$beta_pred <- beta_pred # save coefficients for predictions
  fit$theta <- exp(fit$estimate[length(fit$estimate)])

  
  mu <- exp(X %*% beta_pred)
  fit$predictions <- mu
  fit$se <- sqrt(diag(-1/(fit$hessian)))
  fit$formula <- formula
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum # The log-likelihood of the model
  fit$modelType <- "poisLindRE"
  
  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}
