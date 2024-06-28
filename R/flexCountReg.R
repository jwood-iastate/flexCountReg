#' Estimate flexible count regression models with and without random parameters
#'
#' \code{flexCountReg} is a wrapper function for estimating a variety of count
#' regression models. This is the main function in the package. These include
#' multiple count distributions with options for random parameters.
#'
#' @name flexCountReg
#' @param formula an R formula,
#' @param rpar_formula an optional formula for the model related specifically to
#'   the random parameters (if a random parameter model is desired). If used,
#'   this should not include an outcome variable. If the intercept is random,
#'   include it in this formula. If the intercept is fixed, include it in
#'   \code{formula} but not in \code{rpar_formula}. To remove the intercept, use
#'   \code{0 + vars} or \code{-1 + vars},
#' @param data a dataframe that has all of the variables in the \code{formula}
#'   and \code{rpar_formula},
#' @param dist determines the count distribution to be used in the regression
#'   estimation and must be a string. Options include:
#' \code{"NB2"} for the standard negative binomial with a mean value of
#' \eqn{E[Y]=\mu} and variance \eqn{Var[Y]=\mu+\alpha\mu^2},
#' \code{"NB1"} for the negative binomial with a mean value of \eqn{E[Y]=\mu}
#' and variance \eqn{Var[Y]=\mu+\alpha\mu},
#' \code{"NBP"} for the negative binomial with a mean value of \eqn{E[Y]=\mu}
#' and variance \eqn{Var[Y]=\mu+\alpha\mu^{2-P}},
#' \code{"Poisson Lognormal"} for a Poisson-Lognormal model, 
#' \code{"Poisson Lindley"} for a Poisson-Lindley model,
#' \code{"GW"} for a Generalized Waring model, 
#' \code{"PGE"} for a Poisson-Generalized-Exponential model,
#' \code{"PLL"} for a Poisson-Lindley-Lognormal model, 
#' \code{"PLG"} for a Poisson-Lindley-Gamma (i.e., Negative Binomial-Lindley) model,
#' \code{"PInvGaus"} for a Poisson-Inverse-Gaussian model (Type 1), 
#' \code{"PInvGaus-2"} for a Poisson-Inverse-Gaussian model (Type 2), or
#' \code{"Sichel"} for a Sichel model.
#' @param ln.alpha.formula an optional formula for using independent variables
#'   to estimate the natural log of the overdispersion parameter (makes the
#'   model a generalized negative binomial but is not an option for random
#'   parameters models),
#' @param ln.sigma.formula an optional formula for using independent variables
#'   to estimate the natural log of the standard deviation parameter (makes the
#'   model a generalized Poisson-Lognormal).
#' @param ln.scale.formula an optional formula for using independent variables
#'   to estimate the natural log of the scale parameter in the Poisson
#'   Generalized-Exponential model.
#' @param rpardists an optional named vector whose names are the random
#'   parameters and values the distribution (for random parameter models only).
#'   The distribution options include normal ("n"), lognormal ("ln"), triangular
#'   ("t"), and uniform ("u"). If this is not provided, normal distributions are
#'   used for all random coefficients,
#' @param ndraws the number of Halton draws to use for estimating the random
#'   parameters,
#' @param scrambled if the Halton draws should be scrambled or not.
#'   \code{scrambled = FALSE} results in standard Halton draws while
#'   \code{scrambled = FALSE} results in scrambled Halton draws,
#' @param correlated if the random parameters should be correlated
#'   (\code{correlated = FALSE} results in uncorrelated random coefficients,
#'   \code{correlated = TRUE} results in correlated random coefficients),
#' @param method a method to use for optimization in the maximum likelihood
#'   estimation. For options, see \code{\link[maxLik]{maxLik}},
#' @param max.iters the maximum number of iterations to allow the optimization
#'   method to perform,
#' @param start.vals an optional vector of starting values for the regression
#'   coefficients
#' @param print.level determines the level of verbosity for printing details of
#'   the optimization as it is computed. A value of 0 does not print out any
#'   information, a value of 1 prints minimal information, and a value of 2
#'   prints the most information.
#'
#' @include poisLogn.R genWaring.R rpnb.R nbg.R sichel.R poisLind.R poisLindLnorm.R poisLindGamma.R  
#'
#' @examples
#' \donttest{
#'
#' ## Poisson-Lognormal Model
#' data("washington_roads")
#' washington_roads$AADTover10k <- ifelse(washington_roads$AADT>10000,1,0)
#' poislog.mod <- flexCountReg(Total_crashes ~ lnaadt + lnlength + speed50 +
#'                                 ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 dist="Poisson Lognormal",
#'                                 ndraws=500)
#' summary(poislog.mod)
#'
#' ## Poisson-Lindley Model
#' poislind.mod <- flexCountReg(Animal ~ lnaadt + lnlength + speed50 +
#'                             ShouldWidth04 + AADTover10k,
#'                             data=washington_roads, dist="Poisson Lindley",
#'                             method="sann")
#' summary(poislind.mod)
#'
#' ## Poisson-Generalized-Exponential Model
#' poisgenexp.mod <- flexCountReg(Total_crashes ~ lnaadt + lnlength +
#'                                 speed50 + ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 dist="PGE",
#'                                 ndraws = 100, 
#'                                 method = 'nm')
#' summary(poisgenexp.mod)
#'
#'
#' ## Generalized Waring Model
#' genwaring.mod <- flexCountReg(Total_crashes ~ lnaadt + lnlength +
#'                                 speed50 + ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 dist="GW", method="nm")
#' summary(genwaring.mod)
#'
#' ## Poisson-Lindley-Lognormal Model
#' poislindlnorm.mod <- flexCountReg(Animal ~ lnaadt + lnlength + speed50 +
#'                             ShouldWidth04 + AADTover10k,
#'                             data=washington_roads, dist="PLL",
#'                             ndraws=50, method="nm")
#' summary(poislindlnorm.mod)
#'
#' ## Poisson-Lindley-Gamma Model
#' poislindgamma.mod <- flexCountReg(Animal ~ lnaadt + lnlength + speed50 +
#'                             ShouldWidth04 + AADTover10k,
#'                             data=washington_roads, dist="PLG",
#'                             ndraws=50, method="nm")
#' summary(poislindgamma.mod)
#' 
#' # Poisson-Inverse-Guassian Model (Type 1)
#' poisinvgaus.mod <- flexCountReg(Total_crashes ~ lnaadt + lnlength +
#'                                 speed50 + ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 dist="PInvGaus", 
#'                                 method="nm")
#' summary(poisinvgaus.mod)
#'
#' # Poisson-Inverse-Guassian Model (Type 2)
#' poisinvgaus.mod2 <- flexCountReg(Total_crashes ~ lnaadt + lnlength +
#'                                 speed50 + ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 dist="PInvGaus-2", 
#'                                 method="nm")
#' summary(poisinvgaus.mod2)
#'
#' # Sichel Model
#' sichel.mod <- flexCountReg(Total_crashes ~ lnaadt + lnlength +
#'                                 speed50 + ShouldWidth04 + AADTover10k,
#'                                 data=washington_roads,
#'                                 dist="Sichel", 
#'                                 method="bfgs")
#' summary(sichel.mod)
#' }
#'
#' @export
flexCountReg <- function(formula, rpar_formula=NULL, data, dist = "NB2",
                         ln.alpha.formula = NULL,
                         ln.sigma.formula = NULL,
                         ln.scale.formula = NULL,
                         rpardists = NULL,
                         ndraws = 1500, scrambled = FALSE,
                         correlated = FALSE, method = 'BHHH', max.iters = 200,
                         start.vals = NULL, print.level = 0){

  if (dist=="NB1" && !is.null(rpar_formula)){ # random parameter NB-1 model
    model <- rpnb(formula=formula, rpar_formula=rpar_formula, data=data,
                  ndraws=ndraws, scrambled=scrambled, form = 'nb1',
                  correlated=correlated, method=method, max.iters=max.iters,
                  start.vals=start.vals, print.level=print.level)
    model$dist <- dist
    names(model$estimate) <- model$x_names
    model$random <- TRUE
  }
  else if (dist=="NB2" && !is.null(rpar_formula)){ # random parameter NB-2 model
    model <- rpnb(formula=formula, rpar_formula=rpar_formula, data=data,
                  ndraws=ndraws, scrambled=scrambled, form = 'nb2',
                  correlated=correlated, method=method, max.iters=max.iters,
                  start.vals=start.vals, print.level=print.level)
    model$dist <- dist
    names(model$estimate) <- model$x_names
    model$random <- TRUE
  }
  else if (dist=="NBP" && !is.null(rpar_formula)){ # random parameter NB-P model
    model <- rpnb(formula=formula, rpar_formula=rpar_formula, data=data,
                  ndraws=ndraws, scrambled=scrambled, form = 'nbp',
                  correlated=correlated, method=method, max.iters=max.iters,
                  start.vals=start.vals, print.level=print.level)
    model$dist <- dist
    names(model$estimate) <- model$x_names
    model$random <- TRUE
  }
  else if(dist=="NB2"){ # NB-2 model
    model <- nbg(formula=formula, data=data, ln.alpha.formula=ln.alpha.formula,
                 method=method, max.iters=max.iters, form = 'nb2')
    model$dist <- dist
    model$random <- FALSE
  }
  else if(dist=="NB1"){ # NB-1 model
    model <- nbg(formula=formula, data=data, ln.alpha.formula=ln.alpha.formula,
                 method=method, max.iters=max.iters, form = 'nb1')
    model$dist <- dist
    model$random <- FALSE
  }
  else if(dist=="NBP"){ # NB-P model
    model <- nbg(formula=formula, data=data, ln.alpha.formula=ln.alpha.formula,
                 method=method, max.iters=max.iters, form = 'nbp')
    model$dist <- dist
    model$random <- FALSE
  }
  else if(dist=="Poisson Lognormal"){ # Poisson-Lognormal Model
    model <- poisLogn(formula=formula, data=data, ndraws=ndraws,
                      ln.sigma.formula = ln.sigma.formula,
                      method=method, max.iters=max.iters)
    model$dist <- dist
    model$random <- FALSE
  }
  else if(dist=="Poisson Lindley"){ # Poisson-Lindley Model
    model <- poisLind(formula=formula, data=data, method=method,
                      max.iters=max.iters)
    model$dist <- dist
    model$random <- FALSE
  }
  else if(dist=="GW"){ # Generalized Waring Model
    model <- genWaring(formula=formula, data=data, method=method,
                      max.iters=max.iters)
    model$dist <- dist
    model$random <- FALSE
  }
  else if(dist=="PGE"){ # Poisson-Generalized Exponential Model
    model <- poisGE(formula=formula, data=data, 
                    ln.scale.formula = ln.scale.formula,
                             ndraws=ndraws, method=method, max.iters=max.iters)
    model$dist <- dist
    model$random <- FALSE
  }
  else if(dist=="PLL"){ # Poisson-Lindley-Lognormal Model
    model <- poisLindLnorm(formula=formula, data=data,
                             ndraws=ndraws, method=method, max.iters=max.iters)
    model$dist <- dist
    model$random <- FALSE
  }
  else if(dist=="PLG"){ # Poisson-Lindley-Gamma Model
    model <- poisLindGamma(formula=formula, data=data,
                           ndraws=ndraws, method=method, max.iters=max.iters)
    model$dist <- dist
    model$random <- FALSE
  }
  else if(dist=="PInvGaus"){ # Poisson-Inverse-Gaussian Model - Type 1
    model <- poisInvGaus(formula=formula, data=data,
                           method=method, 
                           max.iters=max.iters)
    model$dist <- dist
    model$random <- FALSE
  }
  else if(dist=="PInvGaus-2"){ # Poisson-Inverse-Gaussian Model - Type 2
    model <- poisInvGaus(formula=formula, data=data,
                         form="Type 2", method=method, 
                         max.iters=max.iters)
    model$dist <- dist
    model$random <- FALSE
  }
  else if (dist=="Sichel"){ # Sichel Model
    model <- sichel(formula=formula, data=data,
                         method=method,
                         max.iters=max.iters)
    model$dist <- dist
    model$random <- FALSE
  }
  else{
    print("Please use one of the following options for `dist`: 'NB1', 'NB2', 'NBP', 'Poisson Lognormal', 
          'Poisson Lindley', 'GW', 'PGE', 'PW', 'PIG', 'PLL', 'PLG', 'PInvGaus', 'PInvGaus-2', or 'Sichel'")
    stop()
    }

  return(model)
}



