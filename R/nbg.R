#' Function for estimating a variety of negative binomial models (NB-1, NB-2,
#' NB-P and generalized versions of each)
#'
#' @name nbg
#' @param formula an R formula.
#' @param data a dataframe that has all of the variables in the \code{formula}
#'   and \code{rpar_formula}. This can be the data used for estimating the model
#'   or another dataframe,
#' @param form the version of the negative binomial to estimate (\code{"nb2"}
#'   estimates the NB-2, \code{"nb1"} estimates the NB-1, \code{"nbp"} estimates
#'   the NB-P)
#' @param ln.alpha.formula an optional formula for using independent variables
#'   to estimate the natural log of the overdispersion parameter (makes the
#'   model a generalized negative binomial).
#' @param method a method to use for optimization in the maximum likelihood
#'   estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param max.iters the maximum number of iterations to allow the optimization
#'   method to perform.
#' @param weights the name of the weighting variable. This is an optional 
#'   parameter for weighted regression.
#' @param bootstraps Optional integer specifying the number of bootstrap 
#'   samples to be used for estimating bootstrapped standard errors. If this is 
#'   used, bootstrapped standard errors will be calculated.
#'
#' @details
#' The NB-1, NB-2, and NB-P versions of the negative binomial distribution are 
#' based on Greene (2008).  The details of each of these are provided below.
#' 
#' # NB-1 Distribution
#' The PMF and log-likelihood functions are:
#' \deqn{P(Y = y) = \frac{\Gamma(y + \frac{\mu}{\alpha})}{y! \, \Gamma(\frac{\mu}{\alpha})} \left( \frac{\frac{\mu}{\alpha}}{\frac{\mu}{\alpha} + \mu} \right)^{\frac{\mu}{\alpha}} \left( \frac{\mu}{\frac{\mu}{\alpha} + \mu} \right)^y}
#' \deqn{LL_{\text{NB1}}(\beta, \alpha) = \sum_{i=1}^n \left[ \ln \Gamma\left( y_i + \frac{\mu_i}{\alpha} \right) - \ln \Gamma\left( \frac{\mu_i}{\alpha} \right) - \ln y_i! + \frac{\mu_i}{\alpha} \ln \left( \frac{\frac{\mu_i}{\alpha}}{\frac{\mu_i}{\alpha} + \mu_i} \right) + y_i \ln \left( \frac{\mu_i}{\frac{\mu_i}{\alpha} + \mu_i} \right) \right]}
#' 
#' The Gradients are:
#' \deqn{\frac{\partial LL_{\text{NB1}}}{\partial \beta} = \sum_{i=1}^n \left[ \left( \psi\left(y_i + \frac{\mu_i}{\alpha}\right) - \psi\left(\frac{\mu_i}{\alpha}\right) + \ln \left( \frac{\frac{\mu_i}{\alpha}}{\frac{\mu_i}{\alpha} + \mu_i} \right) - \frac{y_i}{\frac{\mu_i}{\alpha} + \mu_i} \right) \cdot \mu_i \cdot X_i \right]}
#' \deqn{\frac{\partial LL_{\text{NB1}}}{\partial \eta} = \alpha \sum_{i=1}^n \left[ \psi\left(y_i + \frac{\mu_i}{\alpha}\right) - \psi\left(\frac{\mu_i}{\alpha}\right) + \ln \left( \frac{\frac{\mu_i}{\alpha}}{\frac{\mu_i}{\alpha} + \mu_i} \right) - \frac{\mu_i}{\frac{\mu_i}{\alpha} + \mu_i} \right]}
#' 
#' # NB-2 Distribution
#' The PMF and log-likelihood functions are:
#' \deqn{P(Y = y) = \frac{\Gamma(y + \alpha)}{y! \, \Gamma(\alpha)} \left( \frac{\alpha}{\alpha + \mu} \right)^\alpha \left( \frac{\mu}{\alpha + \mu} \right)^y}
#' \deqn{LL_{\text{NB2}} = \sum_{i=1}^n \left[ \ln \Gamma(y_i + \alpha) - \ln \Gamma(\alpha) - \ln y_i! + \alpha \ln \left( \frac{\alpha}{\alpha + \mu_i} \right) + y_i \ln \left( \frac{\mu_i}{\alpha + \mu_i} \right) \right]}
#' 
#' The Gradients are:
#' \deqn{\frac{\partial LL_{\text{NB2}}}{\partial \beta} = \sum_{i=1}^n \left[ \left( \psi(y_i + \alpha) - \psi(\alpha) + \ln \left( \frac{\alpha}{\alpha + \mu_i} \right) - \frac{y_i}{\alpha + \mu_i} \right) \cdot \mu_i \cdot X_i \right]}
#' \deqn{\frac{\partial LL_{\text{NB2}}}{\partial \eta} = \alpha \sum_{i=1}^n \left[ \psi(y_i + \alpha) - \psi(\alpha) + \ln \left( \frac{\alpha}{\alpha + \mu_i} \right) - \frac{\mu_i}{\alpha + \mu_i} \right]}
#' 
#' # NB-P distribution
#' The PMF and log-likelihood functions are:
#' \deqn{P(Y = y) = \frac{\Gamma(y + \frac{\mu^{2-p}}{\alpha})}{y! \, \Gamma(\frac{\mu^{2-p}}{\alpha})} \left( \frac{\frac{\mu^{2-p}}{\alpha}}{\frac{\mu^{2-p}}{\alpha} + \mu} \right)^{\frac{\mu^{2-p}}{\alpha}} \left( \frac{\mu}{\frac{\mu^{2-p}}{\alpha} + \mu} \right)^y}
#' \deqn{LL_{\text{NBP}}(\beta, \alpha, p) = \sum_{i=1}^n \left[ \ln \Gamma\left( y_i + \frac{\mu_i^{2-p}}{\alpha} \right) - \ln \Gamma\left( \frac{\mu_i^{2-p}}{\alpha} \right) - \ln y_i! + \frac{\mu_i^{2-p}}{\alpha} \ln \left( \frac{\frac{\mu_i^{2-p}}{\alpha}}{\frac{\mu_i^{2-p}}{\alpha} + \mu_i} \right) + y_i \ln \left( \frac{\mu_i}{\frac{\mu_i^{2-p}}{\alpha} + \mu_i} \right) \right]}
#' 
#' The Gradients are:
#' \deqn{\frac{\partial LL_{\text{NBP}}}{\partial \beta} = \sum_{i=1}^n \left[ \left( \psi\left(y_i + \frac{\mu_i^{2-p}}{\alpha}\right) - \psi\left(\frac{\mu_i^{2-p}}{\alpha}\right) + \ln \left( \frac{\frac{\mu_i^{2-p}}{\alpha}}{\frac{\mu_i^{2-p}}{\alpha} + \mu_i} \right) - \frac{y_i}{\frac{\mu_i^{2-p}}{\alpha} + \mu_i} \right) \cdot \mu_i \cdot X_i \right]}
#' \deqn{\frac{\partial LL_{\text{NBP}}}{\partial \eta} = \alpha \sum_{i=1}^n \left[ \psi\left(y_i + \frac{\mu_i^{2-p}}{\alpha}\right) - \psi\left(\frac{\mu_i^{2-p}}{\alpha}\right) + \ln \left( \frac{\frac{\mu_i^{2-p}}{\alpha}}{\frac{\mu_i^{2-p}}{\alpha} + \mu_i} \right) - \frac{\mu_i^{2-p}}{\frac{\mu_i^{2-p}}{\alpha} + \mu_i} \right]}
#' \deqn{\frac{\partial LL_{\text{NBP}}}{\partial p} = \sum_{i=1}^n \left[ \left( \psi\left(y_i + \frac{\mu_i^{2-p}}{\alpha}\right) - \psi\left(\frac{\mu_i^{2-p}}{\alpha}\right) + \ln \left( \frac{\frac{\mu_i^{2-p}}{\alpha}}{\frac{\mu_i^{2-p}}{\alpha} + \mu_i} \right) - \frac{\mu_i}{\frac{\mu_i^{2-p}}{\alpha} + \mu_i} \right) \cdot \frac{\partial}{\partial p} \left( \frac{\mu_i^{2-p}}{\alpha} \right) \right]}
#'
#' @import maxLik  stats modelr MASS
#' @importFrom purrr map map_df
#' @importFrom broom tidy
#' @importFrom dplyr group_by %>% summarise
#' @importFrom tibble deframe
#' 
#' @references 
#' Greene, W. (2008). Functional forms for the negative binomial model for count data. Economics Letters, 99(3), 585-590.
#' 
#' @examples
#'
#' ## NB-P model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
#'
#' nbp.base <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                     ShouldWidth04 + AADTover10k,
#'                     data=washington_roads, form = 'nbp', method = 'NM',
#'                     max.iters=3000)
#' summary(nbp.base)
#'
#' ## Generalized NB-P model
#'
#' nbp.overdispersion <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                                 ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 form = 'nbp',
#'                                 method = 'NM',
#'                                 max.iters=3000,
#'                                 ln.alpha.formula = ~ 1+lnlength)
#' summary(nbp.overdispersion)
#'
#' ## NB-1 Model with Bootstrapping
#' nb1.base <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                         ShouldWidth04 + AADTover10k,
#'                         data=washington_roads, form = 'nb1',
#'                         method = 'BHHH',
#'                         max.iters=3000,
#'                         bootstraps=100)
#' summary(nb1.base)
#'
#' ## Generalize NB-1 Model
#' nb1.overdispersion <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                         ShouldWidth04 + AADTover10k,
#'                         data=washington_roads, form = 'nb1',
#'                         method = 'NM',
#'                         max.iters=3000, ln.alpha.formula = ~ 1+lnlength)
#' summary(nb1.overdispersion)
#'
#' ## NB-2 Model
#' nb2.base <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                         ShouldWidth04 + AADTover10k,
#'                         data=washington_roads, form = 'nb2',
#'                         method = 'NM',
#'                         max.iters=3000)
#' summary(nb2.base)
#'
#' ## Generalize NB-2 Model
#' nb2.overdispersion <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                         ShouldWidth04 + AADTover10k,
#'                         data=washington_roads, form = 'nb2',
#'                         method = 'NM',
#'                         max.iters=3000, ln.alpha.formula = ~ 1+lnlength)
#' summary(nb2.overdispersion)
#' 
#' @export
nbg <- function(formula, data, form = 'nb2', ln.alpha.formula = NULL, 
                method = 'BHHH', max.iters=200, weights=NULL, 
                bootstraps = NULL, print.level = 0) {
  
  if (!form %in% c('nb1', 'nb2', 'nbp')) {
    stop("Invalid form specified")
  }
  
  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  y <- as.numeric(stats::model.response(mod_df))
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
  
  # Use the NB2 from MASS as starting values
  p_model <- MASS::glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  a <- log(1/p_model$theta)
  
  if (is.null(ln.alpha.formula)){
    if (!is.null(a)) {
      comb_start <- append(start, a)
      alpha_X <- NULL
      alpha_X_cols <- 1
    }
  }else{
    alpha_X <- stats::model.matrix(ln.alpha.formula, data)
    alpha_X_cols <- ncol(alpha_X)
    alpha_names <- colnames(alpha_X)
    a_coefs <- rep(0, length(alpha_names))
    a_coefs[1] <- a
    if (!is.null(a_coefs)) {
      comb_start <- c(start, a_coefs)
    }
  }
  
  if (is.null(ln.alpha.formula)) {
    x_names <- append(x_names, 'ln(alpha)')
  } else {
    x_names <- append(x_names, paste('ln(alpha):', alpha_names))
  }
  
  if (form=='nbp') {
    full_start <- append(comb_start, 1.5)
    x_names <- append(x_names, 'P')
  } else{
    full_start <- comb_start
  }
  
  names(full_start) <- x_names
  
  modparams <- as.numeric(ncol(X)) # save the number of model coefficients, not including alpha or P
  
  nbp_prob <- function(observed, predicted, log_alpha, p) {
    alpha <- exp(log_alpha)
    mu <- predicted
    y <- observed
    if (form=='nb2'){
      return(stats::dnbinom(y, size = alpha, mu = mu))
    } else if (form=='nb1'){
      return(stats::dnbinom(y, size = mu/alpha, mu = mu))
    } else{
      return(stats::dnbinom(y, size = (mu^(2-p))/alpha, mu = mu))
    }
  }
  
  reg.run <- function(params, y, X, alpha_X){
    coefs <- params[1:modparams]
    if (is.null(alpha_X)) {
      log_alpha <- params[(modparams + 1)]
    } else {
      alpha_coefs <- params[(modparams + 1):(modparams+alpha_X_cols)]
      
      log_alpha <- alpha_X %*% alpha_coefs
    }
    alpha <- exp(log_alpha)
    mu <- exp(X %*% coefs)
    
    if (form=='nbp') p <- params[length(params)] else p <- NULL
    
    predicted <- exp(X %*% coefs)
    
    probs <- nbp_prob(y, predicted, log_alpha, p)
    
    ll <- sum(weights*log(probs)) # accounting for weights (values of 1 if not provided)
    if (method == 'bhhh' | method == 'BHHH'){
      return(weights*log(probs))
    } else{return(ll)}
  }
  
  # Run the maximum likelihood estimation
  fit <- maxLik::maxLik(reg.run,
                        start = full_start,
                        y = y,
                        X = X,
                        alpha_X = alpha_X,
                        method = method,
                        control = list(iterlim = max.iters, 
                                       printLevel = print.level))
  
  beta_est <- fit$estimate
  alpha_coefs <- beta_est[(modparams + 1):(modparams + alpha_X_cols)]
  
  if (form=="nbp") p <- beta_est[length(beta_est)] else p <- NULL
  
  beta_pred <- beta_est[1:modparams]
  
  names(fit$estimate) <- x_names
  fit$beta_pred <- beta_pred # save coefficients for predictions
  fit$formula <- formula
  fit$form <- form
  fit$ln_alpha <- alpha_coefs # natural log of alpha or the parameters for it
  fit$p <- p
  fit$modelType <- "nbg"
  fit$predictions <- exp(X %*% beta_pred)
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum # The log-likelihood of the model
  
  # Optionally, compute bootstrapped standard errors
  # create function to clean data and run maxLik
  nbg.boot <- function(data){
    mod1_frame <- stats::model.frame(formula, data)
    X_Fixed <- stats::model.matrix(formula, data)
    y <- stats::model.response(mod1_frame)
    
    if (is.null(ln.alpha.formula)){
      if (!is.null(a)) {
        alpha_X <- NULL
      }
    }else{
      alpha_X <- stats::model.matrix(ln.alpha.formula, data)
    }
    
    int_res <-  maxLik::maxLik(reg.run,
                               start = fit$estimate,
                               y = y,
                               X = X,
                               alpha_X = alpha_X,
                               method = method,
                               iterlim = max.iters)
    return(int_res)
  }
  
  
  if (!is.null(bootstraps) & is.numeric(bootstraps)) {
    bs.data <- modelr::bootstrap(data, n = bootstraps)
    models <- map(bs.data$strap, ~ nbg.boot(data = .))
    tidied <- map_df(models, broom::tidy, .id = "id")
    
    SE <- tidied %>%
      group_by(term) %>%
      summarise(sd = sd(estimate)) %>% deframe()
    
    fit$bootstrapped_se <- SE[names(fit$estimate)]
  }
  
  fit$bootstraps = if (!is.null(bootstraps)) bootstraps else NULL
  
    obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}
