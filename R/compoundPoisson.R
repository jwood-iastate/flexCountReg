#' Function for estimating a various compound Poisson regression models where the error term follows specific distributions using Maximum Simulated Likelihood when there is not a closed form PMF
#'
#' @name compoundPoisson
#' @param formula an R formula.
#' @param method a method to use for optimization in the maximum likelihood estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param data a dataframe that has all of the variables in the \code{formula} and \code{rpar_formula}.
#' @param distrib a distribution used with the Poisson distribution to create the compound Poisson distribution. Options include:
#' \code{"GE"} for the Generalized Exponential distribution (resulting in the Poisson-Generalized-Exponential model),
#' \code{"W"} for the weibull distribution (resulting in the Poisson-Weibull model),
#' \code{"IG"} for the Inverse Gamma distribution (resulting in the Poisson-Inverse-Gamma model)
#' \code{"Weibull"} for the weibull distribution (resulting in the Poisson-Weibull model)
#' @param ndraws the number of Halton draws to use for estimating the random parameters,
#' @param max.iters the maximum number of iterations to allow the optimization method to perform.
#'
#' @details
#' For the distributions used for compounding with the Poisson, the CDF is used to determine the quantile function. These include:
#'
#' Generalized Exponential Distribution CDF and Quantile Functions (see Gupta & Kundu (2007))
#' \deqn{F(x)=\frac{\left(e^{\beta x}-1\right)^\alpha-1}{\alpha\beta}}
#' \deqn{Q(p)=\frac{ln\left(\left(P\cdot\alpha+1\right)^{\frac{1}{\alpha}}+1\right)}{\beta}}
#' Where \eqn{\beta} is the shape parameter, \eqn{\alpha} is the scale parameter, and \eqn{\P} is the percentile. The parameter constraints are \eqn{\tau>0} and \eqn{\alpha>0}.
#'     The mean and variance of this distribution are:
#' \deqn{\mu=\frac{1}{\beta}\left(\psi(\alpha+1)-\psi(1)\right)}
#' \deqn{\sigma^2=\frac{1}{\beta}\left(\psi\prime(1)-\psi\prime(\alpha+1)\right)}
#' Where \eqn{\psi} is the digamma function and \eqn{\psi\prime} is the trigamma function
#'     The resulting Poisson-Generalized Exponential distribution has the mean and variance:
#' \deqn{E[Y]=\mu=e^{X\beta}}
#' \deqn{Var[y]=\sigma^2=\mu+\left(\frac{\psi\prime(1)-\psi\prime(\alpha+1)}{\beta\left(\psi(\alpha+1)-\psi(1)\right)^2}\right)\mu^2}
#'
#' Weibull Distribution  CDF and Quantile Functions
#' \deqn{F(x)=1-exp(-(x/\beta)^\alpha)}
#' \deqn{Q(p)=\beta(-ln(1-P))^{1/\alpha}}
#' Where \eqn{\beta} is the shape parameter, \eqn{\alpha} is the scale parameter, and \eqn{P} is the percentile. There are constraints for \eqn{\alpha>0} and \eqn{\beta>0}.
#'     The mean and variance of this distribution are:
#' \deqn{\mu=\beta\Gamma(1+1/\alpha)}
#' \deqn{\sigma^2=\beta^2\left(\Gamma(1+2/\alpha)-\left(\Gamma(1+1/\alpha)\right)^2\right)}
#'     The resulting Poisson-Weibull distribution has the mean and variance:
#' \deqn{E[Y]=\mu=e^{X\beta}}
#' \deqn{Var[y]=\sigma^2=\mu+\left(\frac{\Gamma\left(\frac{\alpha+2}{\alpha}\right)}{\Gamma\left(1+\frac{1}{\alpha}\right)^2}-1\right)\mu^2}
#'
#' Inverse-Gamma  CDF Function
#' \deqn{F(x)=\frac{\Gamma\left(\alpha,\frac{\beta}{x}\right)}{\Gamma(\alpha)}}
#' Where \eqn{\beta} is the shape parameter and \eqn{\alpha} is the scale parameter. For the variance to be defined, \eqn{\alpha>2}. Additionally, an additional constraint is \eqn{\beta>0}.
#'     The mean and variance of this distribution are:
#' \deqn{\mu=\frac{\beta}{\alpha-1}}
#' \deqn{\sigma^2=\frac{\beta^2}{(\alpha-1)^2(\alpha-2)}}
#'     The resulting Poisson-Inverse-Gamma distribution has the mean and variance:
#' \deqn{E[Y]=\mu=e^{X\beta}}
#' \deqn{Var[y]=\sigma^2=\mu+\frac{\beta^2}{\alpha-2}\mu^2}
#'
#' @references
#' Gupta, R. D., & Kundu, D. (2007). Generalized exponential distribution: Existing results and some recent developments. Journal of Statistical planning and inference, 137(11), 3537-3547.
#'
#'
#' @import nlme maxLik MASS stats randtoolbox invgamma modelr
#' @importFrom utils head  tail
#'
#' @examples
#' \donttest{
#'
#' ## Poisson-Generalized-Exponential Model
#' data("washington_roads")
#' poisgenexp.mod <- compoundPoisson(Total_crashes ~ lnaadt + lnlength,
#'                                 data=washington_roads,
#'                                 distrib="GE",
#'                                 ndraws=500, method="bfgs")
#' summary(poisgenexp.mod)
#'
#' ## Poisson-Weibull Model
#' poisweib.mod <- compoundPoisson(Total_crashes ~ lnaadt + lnlength,
#'                                 data=washington_roads,
#'                                 distrib="W",
#'                                 ndraws=500, method="bfgs")
#' summary(poisweib.mod)
#'
#' ## Poisson-Inverse-Gamma Model
#' poisginvgamma.mod <- compoundPoisson(Total_crashes ~ lnaadt + lnlength,
#'                                 data=washington_roads,
#'                                 distrib="IG",
#'                                 ndraws=500, method="bfgs")
#' summary(poisginvgamma.mod)}
#' @export
compoundPoisson <- function(formula, data, distrib="GE", ndraws = 1500, method = 'BHHH', max.iters = 1000) {

  mod_df <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  y <- as.numeric(stats::model.response(mod_df))
  x_names <- colnames(X)

  p_model <- MASS::glm.nb(formula, data)
  start <- as.vector(p_model$coefficients)

  if(distrib=="GE"){
    start <- append(start, 0) # add start value for ln(beta) - ensures beta>0
    start <- append(start, 0) # add star value for ln(alpha) - ensures alpha>0
  }
  else if (distrib=="W"){
    start <- append(start, 0) # add start value for ln(beta) - ensures beta>0
    start <- append(start, 0) # add star value for ln(alpha) - ensures alpha>0
  }
  else if (distrib=="IG"){
    start <- append(start, 0) # add start value for ln(beta) - ensures beta>0
    start <- append(start, 1) # add star value for ln(alpha-2) - ensures alpha>2
  }
  else{
    print("Please use a `distrib` value from one of the following: 'GE', 'W', or 'IG'")
    quit()
  }

  n <- as.integer(ndraws)

  poisson_prob <- function(observed, predicted) {
    log_probs <- dpois(observed, predicted, log = TRUE)

    return(log_probs)
  }

  p_poisson_lognormal <- function(p, y, X, n, est_method){

    coefs <- as.array(head(p,-2))
    distpars <- tail(p,2)
    par1 <- exp(distpars[1])
    par2 <- exp(distpars[2])
    h = randtoolbox::halton(n)

    probs <- as.vector(rep(0, length(y)))

    if(distrib=="GE"){
      i_logdist <- log((log(h*par2+1)^(1/par2)+1)/par1)
    }
    else if (distrib=="W"){
      i_logdist <- log(stats::qweibull(h, shape=par1, scale=par2))
    }
    else if (distrib=="IG"){
      i_logdist <- log(invgamma::qinvgamma(h, shape=par1, scale=(par2-2)))
    }

    for (i in i_logdist){
      coefs_i <- coefs
      coefs_i[1] <- coefs_i[1] + i
      mu <- exp(X %*% coefs_i)
      p_prob <- poisson_prob(y, mu)
      probs <- probs + exp(p_prob)/n
    }

    ll <- sum(log(probs))
    if (est_method == 'bhhh' | est_method == 'BHHH'){
      return(log(probs))
    } else{return(ll)}
  }

  fit <- maxLik::maxLik(p_poisson_lognormal,
                        start = start,
                        y = y,
                        X = X,
                        n = n,
                        est_method = method,
                        method = method,
                        control = list(iterlim = max.iters))

  beta_est <- fit$estimate
  beta_pred <- head(beta_est, -2)
  dispars <- tail(beta_est, 2)
  x_names <- append(x_names, 'ln(beta)')
  if (distrib=="IG"){
    x_names <- append(x_names, 'ln(alpha-2)')
  }
  else{
    x_names <- append(x_names, 'ln(alpha)')
  }

  names(fit$estimate) <- x_names

  fit$formula <- formula
  fit$beta_pred <- beta_pred
  fit$beta <- exp(dispars[1])
  if (distrib=="IG"){
    fit$alpha <- exp(dispars[2])+2
  }
  else{
    fit$alpha <- exp(dispars[2])
  }
  fit$distrib <- distrib
  fit$modelType <- "compoundPoisson"
  fit$predictions <- exp(X %*% beta_pred)
  fit$observed <- y
  fit$residuals <- y - fit$predictions
  # Note that this has the predictions, residuals, and observed outcome stored with the model
  
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
  print('The Likelihood Ratio (LR) Test for H0: Compound Poisson is No Better than the Poisson')
  print(paste('LR = ', round(fit$LR,4)))
  print(paste('LR degrees of freedom = ', fit$LRdof))
  print(paste('LR p-value = ', LRpval))
  print(paste("Macfadden's Pseudo R^2 = ", round(fit$PseudoR2,4)))

  return(fit)
}


