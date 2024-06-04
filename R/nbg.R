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
#'
#' @details
#' The NB-1, NB-2, and NB-P versions of the negative binomial distribution
#'
#'
#' @import nlme  maxLik  MASS  stats modelr
#' @export
#' @examples
#' \donttest{
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
#' ## NB-1 Model
#' nb1.base <- nbg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                         ShouldWidth04 + AADTover10k,
#'                         data=washington_roads, form = 'nb1',
#'                         method = 'NM',
#'                         max.iters=3000)
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
#' summary(nb2.overdispersion)}
nbg <- function(formula, data, form = 'nb2', ln.alpha.formula = NULL, method = 'BHHH', max.iters=200) {

  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- colnames(X)

  # Use the Poisson as starting values
  p_model <- MASS::glm.nb(formula, data = data)
  start <- unlist(p_model$coefficients)
  a <- log(1/p_model$theta)

  if (is.null(ln.alpha.formula)){
    if (!is.null(a)) {
      comb_start <- append(start, a)
    }
  }else{
    alpha_X <- stats::model.matrix(ln.alpha.formula, data)
    alpha_names <- nlme::Names(ln.alpha.formula, data)
    a_coefs <- rep(0, length(alpha_names))
    a_coefs[1] = a
    if (!is.null(a_coefs)) {
      comb_start <- c(start, a_coefs)
    }
  }

  if (!is.null(comb_start) & form=='nbp') {
    full_start <- append(comb_start, 1.5)
  } else{
    full_start <- comb_start
  }

  modparams <- as.numeric(length(start)) # save the number of model coefficients, not including alpha or P

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

  reg.run <- function(beta, y, X, est_method){

    params_split <- split(beta,ceiling(seq_along(beta) / modparams))

    coefs <- as.vector(unlist(params_split[1]))
    dist_params <- as.vector(unlist(params_split[2]))

    if (form=='nbp'){
      if (is.null(ln.alpha.formula)){
        log_alpha <- dist_params[1]
        p <- dist_params[2]
      }else{
        alpha_coefs <- dist_params[1:(length(dist_params)-1)]
        log_alpha <- alpha_X %*% alpha_coefs
        p <- dist_params[length(dist_params)]
      }
    } else{
      if (is.null(ln.alpha.formula)){
        log_alpha <- dist_params[1]
      }else{
        alpha_coefs <- dist_params
        log_alpha <- alpha_X %*% alpha_coefs
      }
      p <- NULL
    }

    predicted <- exp(X %*% coefs)

    probs <- nbp_prob(y, predicted, log_alpha, p)

    ll <- sum(log(probs))
    if (est_method == 'bhhh' | est_method == 'BHHH'){
      return(log(probs))
    } else{return(ll)}
  }

  fit <- maxLik::maxLik(reg.run,
                start = full_start,
                y = y,
                X = X,
                est_method = method,
                method = method,
                iterlim = max.iters)

  beta_est <- fit$estimate

  betas <- split(beta_est,ceiling(seq_along(beta_est) / modparams))
  beta_pred <- as.vector(unlist(betas[1]))

  if (is.null(ln.alpha.formula)){
    x_names <- append(x_names, 'ln(alpha)')
  }else{
    x_names <- append(x_names, paste('ln(alpha): ', alpha_names))
  }

  if (form=='nbp'){
    x_names = append(x_names, 'P')
  }

  names(fit$estimate) <- x_names
  fit$beta_pred <- beta_pred # save coefficients for predictions
  fit$formula <- formula
  fit$form <- form
  fit$modelType <- "nbg"
  fit$predictions <- exp(X %*% beta_pred)
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  fit$LL <- fit$maximum # The log-likelihood of the model

  # Estimate Poisson model for tests and pseudo R^2
  pois_mod <- glm(formula, data, family = poisson(link = "log"))
  base_mod <- glm(y ~ 1, family = poisson(link = "log"))

  LLpoisson <- sum(dpois(pois_mod$y, pois_mod$fitted.values, log=TRUE))
  LLbase <- sum(dpois(base_mod$y, base_mod$fitted.values, log=TRUE))

  fit$LR <- -2*(LLpoisson - fit$LL) # LR Statistic
  fit$LRdof <- length(x_names) - length(pois_mod$coefficients) # LR Degrees of Freedom
  if (fit$LR>0) {
    fit$LR_pvalue <- pchisq(fit$LR, fit$LRdof, lower.tail=FALSE)  # LR p-Value
  }else{
    fit$LR_pvalue <- 1
  }

  # Compute McFadden's Pseudo R^2, based on a Poisson intercept-only model
  fit$PseudoR2 <- 1-fit$LL/LLbase


  # Print out key model metrics
  LRpval <- ifelse(fit$LR_pvalue<0.0001, "<0.0001", round(fit$LR_pvalue,4))
  print('The Likelihood Ratio (LR) Test for H0: NB is No Better than Poisson')
  print(paste('LR = ', round(fit$LR,4)))
  print(paste('LR degrees of freedom = ', fit$LRdof))
  print(paste('LR p-value = ', LRpval))
  print(paste("Macfadden's Pseudo R^2 = ", round(fit$PseudoR2,4)))

  # Note that this has the predictions, residuals, observed outcome, LR test, and pseudo-r^2 stored with the model

  return(fit)
}
