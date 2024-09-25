#' Count regression models
#' 
#' @name countreg
#' @aliases countreg
#' 
#' @param formula a symbolic description of the model to be fitted.
#' @param data a data frame containing the variables in the model.
#' @param family the name of the distribution/model type to estimate. The 
#'        default "NB2" is the standard negative binomial distribution with a 
#'        log link. other options are listed below. 
#' @param offset the name of a variable in the data frame that should be used 
#'        as an offset (i.e., included but forced to have a coefficient of 1).
#' @param weights the name of a variable in the data frame that should be used
#'        as a frequency weight.
#' @param verbose an optional parameter. If `TRUE`, the function will print out 
#'        the progress of the model fitting. Default is `FALSE`.
#' @param dis_param_formula_1 a symbolic description of the model for the 
#'        natural log of the dispersion parameter or first parameter of the 
#'        count distribution used. Further details are provided below.
#' @param dis_param_formula_2 a symbolic description of the model for the second
#'        parameter of the count distribution used. Further details are provided 
#'        below.
#' @param underreport_formula an optional formula to estimate the underreporting 
#'        for any of the count model options. The underreporting is estimated as 
#'        a function (logit or probit) of the predictors in the model. For the 
#'        model to be tractable, the independent variables cannot be the exact 
#'        same as the count model. The default is `NULL`.
#' @param underreport_family the name of the distribution/model type to estimate 
#'        the underreporting portion of the model when `underreport_formula` is 
#'        specified. The default is "logit" for a binary logistic regression 
#'        model. The other option is "probit" for a probit model.
#' @param rpar_obs_formula an R formula with the independent variables that are 
#'        observation-specific random parameters. This should not include an 
#'        outcome variable. If the intercept is random (and varies over all 
#'        observations), include it in this formula. If the intercept is fixed, 
#'        include it in `formula` and not in `rpar_obs_formula`. To remove the 
#'        intercept, use \code{0 + vars} or \code{-1 + vars}.
#' @param rpar_panel_formula an R formula with the independent variables that  
#'        are panel-specific random parameters. This should not include an 
#'        outcome variable. If the intercept is random (and varies across the 
#'        panels), include it in this formula. Only use this if the `panelID` is 
#'        also specified.
#' @param rpar_underreport_formula an R formula with the independent variables   
#'        that are random parameters in the underreporting model. This should  
#'        not include an outcome variable. If the intercept for the 
#'        underreporting model is random include it in this formula and remove 
#'        it from `underreport_formula`. Only use this if the 
#'        `underreport_formula` is also specified.
#' @param rpar_obs_dists an optional named vector whose names are the random 
#'        parameters and distribution to be used for estimating the random 
#'        parameters for the variables in `rpar_obs_formula`. The distribution 
#'        options include normal ("n"), lognormal ("ln"), triangular ("t"), 
#'        uniform ("u"), and gamma ("g"). If this is not provided, normal 
#'        distributions are used for all random coefficients in 
#'        `rpar_obs_formula`. If `correlated`=`TRUE`, the distributions are  
#'        set to the normal distribution, regardless of the distributions  
#'        specified with `rpar_obs_dists`.
#' @param rpar_panel_dists an optional named vector whose names are the random 
#'        parameters and distribution to be used for estimating the random 
#'        parameters for the variables in `rpar_panel_formula`. The distribution 
#'        options include normal ("n"), lognormal ("ln"), triangular ("t"), 
#'        uniform ("u"), and gamma ("g"). If this is not provided, normal 
#'        distributions are used for all random coefficients in 
#'        `rpar_panel_formula`. If `correlated`=`TRUE`, the distributions are 
#'        set to the normal distribution, regardless of the distributions  
#'        specified with `rpar_panel_dists`.
#' @param rpar_underreport_dists an optional named vector whose names are the  
#'        random parameters and distribution to be used for estimating the  
#'        random parameters for the variables in `rpar_underreport_formula`. The  
#'        distribution options include normal ("n"), lognormal ("ln"), 
#'        triangular ("t"),  uniform ("u"), and gamma ("g"). If this is not 
#'        provided, normal distributions are used for all random coefficients in 
#'        `rpar_underreport_formula`. If `correlated`=`TRUE`, the distributions 
#'        are set to the normal distribution, regardless of the distributions  
#'        specified with `rpar_underreport_dists`.
#' @param correlated an optional parameter. When a random parameter model is 
#'        specified, a `TRUE` value indicates that the random parameters should
#'         be correlated. This uses correlated normal distributions with the 
#'         correlation captured using a Cholesky decomposition matrix. Default 
#'         is `FALSE`.
#' @param ndraws The number of Halton draws for integrating the distribution 
#'        being compounded with the Poisson distribution when there is not a 
#'        closed-form solution. This is also used in estimating the random 
#'        parameters for random parameter models. Default is 1500. It is 
#'        recommended to test different numbers of draws to determine if the 
#'        model is stable (i.e., doesn't change or has minimal change as the 
#'        number of draws changes within a reasonable range).
#' @param scrambled TRUE or FALSE. Indicates if the Halton draws should be 
#'        scrambled. This removes correlation between Halton draws with larger 
#'        base prime numbers. Default is TRUE.
#' @param panelID an optional string or vector of strings with the name(s) of 
#'        variables in the data that identify the panels. This is only 
#'        applicable to random effects and random parameter models. Default is 
#'        `NULL`. 
#' @param method Optimization method to be used for maximum likelihood 
#'        estimation. See `maxLik` documentation for options. The default is 
#'        "NM" for the Nelder-Mead method.
#' @param max.iters Maximum number of iterations for the optimization method.
#' @param start.vals Optional vector of starting values for the optimization.
#' @param stderr Type of standard errors to use. The default is "normal". Other 
#'        options include "boot" for bootstrapped standard errors, or "robust" 
#'        for robust standard errors.
#' @param bootstraps Optional integer specifying the number of bootstrap samples 
#'        to be used for estimating standard errors when `stderr`= "boot". Note
#'        that this currently does not work when an offset variable is used.
#' 
#' @description
#' The purpose of this function is to estimate count regression models using 
#' maximum likelihood estimation (MLE) or Maximum Simulated Likelihood 
#' Estimation (MSLE). The function can estimate the following models:
#' \itemize{
#'  \item Poisson (Poisson)
#'  \item Negative Binomial 1 (NB1)
#'  \item Negative Binomial 2 (NB2)
#'  \item Negative Binomial P (NBP)
#'  \item Poisson-Lognormal (PLN)
#'  \item Poisson-Generalized-Exponential (PGE)
#'  \item Poisson-Inverse-Gaussian Type 1 (PIG1)
#'  \item Poisson-Inverse-Gaussian Type 2 (PIG2)
#'  \item Poisson-Lindley (PL)
#'  \item Poisson-Lindley-Gamma (PLG), also known as the Negative 
#'        Binomial-Lindley (NBL)
#'  \item Poisson-Lindley-Lognormal (PLL)
#'  \item Poisson-Weibull (PW)
#'  \item Sichel (SI)
#'  \item Generalized Waring (GW)
#'  \item Conway-Maxwell-Poisson (COM)
#' }
#' 
#' @details
#' For the `family` argument, the following options are available:
#' \itemize{
#'  \item "Poisson" for the Poisson distribution with a log link.
#'  \item "NB1" for Negative Binomial 1 distribution with a log link.
#'  \item "NB2" for Negative Binomial 2 distribution with a log link (i.e., the 
#'        standard negative binomial model).
#'  \item "NBP" for Negative Binomial P distribution with a log link.
#'  \item "PLN" for Poisson-Lognormal distribution with a log link.
#'  \item "PGE" for Poisson-Generalized-Exponential distribution with a log 
#'        link.
#'  \item "PIG1" for Poisson-Inverse-Gaussian Type-1 distribution with a log link.
#'  \item "PIG2" for Poisson-Inverse-Gaussian Type-2 distribution with a log link.
#'  \item "PL" for Poisson-Lindley distribution with a log link.
#'  \item "PLG" for Poisson-Lindley-Gamma distribution with a log link.
#'  \item "PLL" for Poisson-Lindley-Lognormal distribution with a log link.
#'  \item "PW" for Poisson-Weibull distribution with a log link.
#'  \item "SI" for Sichel distribution with a log link.
#'  \item "GW" for Generalized Waring distribution with a log link.
#'  \item "COM" for Conway-Maxwell-Poisson (COM) distribution with a log link.
#'  }
#'  
#'  The `dis_param_formula_1` and `dis_param_formula_2` parameters are used to 
#'  estimate the dispersion parameter or other parameters of the count 
#'  distribution used. This leads to the distributions parameters being 
#'  functions rather than constants in the model. For example, 
#'  if the user wants to estimate the overdispersion parameter of the Negative 
#'  Binomial 2 distribution as a function of the variable `x1` and `x2`, 
#'  the user would specify `dis_param_formula_1 = ~ x1 + x2`. In the case of the 
#'  Negative Binomial distributions, the model is known as a Generalized 
#'  Negative Minomial model when the overdispersion parmater is specified as a 
#'  function.
#'  
#'  The function linking the distribution parameters to the predictors is:
#'  \deqn{Param = \exp((Intercept) + \sum \beta X)}
#'  
#'  The parameters for the different models are as follows:
#'  
#'  For `dis_param_formula_1`, the models are for the parameters:
#'  \itemize{
#'  \item Not Applicable for the Poisson model.
#'  \item \eqn{ln(\alpha)} for the Negative Binomial 1 model.
#'  \item \eqn{ln(\alpha)} for the Negative Binomial 2 model.
#'  \item \eqn{ln(\alpha)} for the Negative Binomial P model.
#'  \item \eqn{ln(\sigma)} for the Poisson-Lognormal model.
#'  \item shape parameter for the Poisson-Generalized-Exponential model.
#'  \item \eqn{ln(\eta)} for the Poisson-Inverse-Gaussian model.
#'  \item \eqn{ln(\theta)} for the Poisson-Lindley model.
#'  \item \eqn{ln(\theta)} for the Poisson-Lindley-Gamma model.
#'  \item \eqn{ln(\theta)} for the Poisson-Lindley-Lognormal model.
#'  \item \eqn{ln(\alpha)} for the Poisson-Weibull model.
#'  \item \eqn{\gamma} for the Sichel model.
#'  \item \eqn{k} for the Generalized Waring model.
#'  \item \eqn{ln(\nu)} for the Conway-Maxwell-Poisson model.
#'  }
#'  
#'  For `dis_param_formula_2`, the models are for the parameters:
#'  \itemize{
#'  \item Not Applicable for the Poisson model.
#'  \item Not Applicable for the Negative Binomial 1 model.
#'  \item Not Applicable for the Negative Binomial 2 model.
#'  \item p for the Negative Binomial P model.
#'  \item Not Applicable for the Poisson-Lognormal model.
#'  \item scale parameter for the Poisson-Generalized-Exponential model.
#'  \item Not Applicable for the Poisson-Inverse-Gaussian model.
#'  \item Not Applicable for the Poisson-Lindley model.
#'  \item \eqn{ln(\alpha)} for the Poisson-Lindley-Gamma model.
#'  \item \eqn{ln(\sigma)} for the Poisson-Lindley-Lognormal model.
#'  \item \eqn{ln(\sigma)} for the Poisson-Weibull model.
#'  \item \eqn{ln(\sigma)} for the Sichel model.
#'  \item \eqn{ln(\rho)} for the Generalized Waring model.
#'  \item Not Applicable for the Conway-Maxwell-Poisson model.
#'  }
#'  
#'  The `ndraws` parameter is used to estimate the distribution when there is 
#'  not a closed-form solution. This uses Halton draws to integrate the 
#'  distribution being compounded with the Poisson distribution. The default is 
#'  1500. The models this is applicable for include:
#'  \itemize{
#'  \item Poisson-Lognormal
#'  \item Poisson-Generalized-Exponential
#'  \item Poisson-Lindley-Gamma (more efficient and stable than using 
#'        hypergeometric functions)
#'  \item Poisson-Lindley-Lognormal
#'  \item Poisson-Weibull
#'  }
#'  The `ndraws` parameter is also used for the estimation of random parameters. 
#'  When there are 4 or more random parameters, it is recommended to use 
#'  `scrambeled` = TRUE to use scrambled Halton draws (which removes correlation 
#'  that is inherent between Halton draws with higher base prime numbers)
#' 
#' @section Random Parameters Model Details
#' For random parameters, the model captures heterogeneity of effects and other 
#' unobserved heterogeneity correlated with specified parameters. These can 
#' accommodate panel (i.e., longitudinal) model specifications. When the model 
#' is based on panel data, random parameters can be allowed to vary across 
#' individuals (i.e., each individual has a constant coefficient/effect related 
#' to that variable but the coefficient varies across individuals) or within 
#' individuals (where every observation, regardless of the individual, has 
#' different coefficient values). The differences for this are captured using 
#' the model parameters `rpar_panel_formula` for the random parameters that vary 
#' across individuals (only applicable if the panel ID is also specified) and 
#' `rpar_obs_formula` for the random parameters that vary across all 
#' observations (i.e., within individuals). 
#' 
#' Additionally, random parameters for the underreporting model (if specified) 
#' can also be estimated using `rpar_underreport_formula`. When this is used 
#' with a panel model specification, the random parameters for underreporting 
#' are treated as across-individual random parameters. This increases the 
#' tractability of the model. 
#' 
#' For the random parameters, they can be correlated or uncorrelated. If they 
#' correlated (specified with `correlated`=`TRUE`), the random parameters all 
#' follow correlated normal distributions. If they are not correlated, the 
#' distributions can be specified for the variables in each of the random 
#' parameters formulas using `rpar_obs_dists` for `rpar_obs_formula`, 
#' `rpar_panel_dists` for `rpar_panel_formula`, and `rpar_underreport_dists` 
#' for `rpar_underreport_formula`. Options for the random parameter 
#' distributions include:
#'  \itemize{
#'  \item "n" for the normal distribution
#'  \item "ln" for the lognormal distribution
#'  \item "t" for the symmetrical triangular distribution 
#'  \item "u" for the uniform distribution
#'  \item "g" for the gamma distribution
#'  }
#' The normal distribution is most commonly used. The lognormal and gamma 
#' distributions are useful for coefficients that have a theoretical reason that 
#' the coefficients should always be positive or always be negative. The uniform 
#' distribution often works well for indicator variables.
#' 
#' To specify the distributions, provide a named vector of coefficient names 
#' with the associated distributions. For instance, if:
#' 
#' `rpar_obs_formula` ~ 1 + x1 + x2 
#' 
#' And the distributions desired are uniform for the intercept, gamma for x1, 
#' and normal for x2, the distributions would be specified as 
#' 
#' `rpar_obs_dists` = c(intercept="u", x1="g", x2="n")
#' 
#' **Special Case Warning**
#' When the intercept is included in `rpar_obs_formula` and is normally 
#' distributed and the `family`="PL", the results will be intractable. The 
#' Poisson-Lognormal model is equivalent to a random parameters Poisson model 
#' that has a normally distributed random intercept that varies across all 
#' observations. 
#' 
#' @return
#' An object of class `countreg` which is a list with the following components:
#' \itemize{
#'  \item model: the fitted model object.
#'  \item data: the data frame used to fit the model.
#'  \item call: the matched call.
#'  \item formula: the formula used to fit the model.
#' }
#' 
#' @import modelr randtoolbox  plm
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom purrr map map_df compact
#' @importFrom broom tidy
#' @importFrom dplyr group_by %>% summarise
#' @importFrom tibble deframe
#' @importFrom maxLik maxLik
#' @importFrom stringr str_replace_all
#' @importFrom sandwich sandwich
#' @include pinvgaus.R ppoislogn.R plindLnorm.R plindGamma.R psichel.R Generalized-Waring.R ppoisGE.R psichel.R plind.R helpers.R
#' 
#' @examples
#' # Load the Washington data
#' data("washington_roads")
#' washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
#' # Estimate an NB2 model with a dispersion parameter as a function of the 
#' # variable `speed50` (i.e., generalized NB2), verbose output, and use the 
#' # BFGS optimization method
#' nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'                data = washington_roads, family = "NB2",
#'                dis_param_formula_1 = ~ speed50, verbose = TRUE, method='BFGS')
#' summary(nb2)
#' 
#' # Estimate a Poisson-Lognormal model (a low number of draws is used to speed 
#' # up the estimation for examples - not recommended in practice)
#' pln <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'               data = washington_roads, family = "PLN", ndraws=10)
#' summary(pln)  
#' 
#' # Estimate an Poisson-Lognormal with underreporting (probit)
#' plogn_underreport <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'               data = washington_roads, family = "NB2",
#'               underreport_formula = ~ speed50 + AADT10kplus, underreport_family = "probit")
#' summary(plogn_underreport)
#' 
#' # Estimate a Conwa-Maxwell-Poisson model
#' com_model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'               data = washington_roads, family = "COM")
#' summary(com_model)
#' 
#' @references
#' Greene, W. (2008). Functional forms for the negative binomial model for count data. Economics Letters, 99(3), 585-590.
#' 
#' Pararai, M., Famoye, F., & Lee, C. (2006). Generalized Poisson regression model for underreported counts. Advances and applications in Statistics, 6(3), 305-322.
#' 
#' Pararai, M., Famoye, F., & Lee, C. (2010). Generalized poisson-poisson mixture model for misreported counts with an application to smoking data. Journal of Data Science, 8(4), 607-617.
#' 
#' Rigby, R. A., Stasinopoulos, D. M., & Akantziliotou, C. (2008). A framework for modelling overdispersed count data, including the Poisson-shifted generalized inverse Gaussian distribution. Computational Statistics & Data Analysis, 53(2), 381-393.
#' 
#' Wood, J.S., Eric T. Donnell, & Christopher J. Fariss. "A method to account for and estimate underreporting in crash frequency research." Accident Analysis & Prevention 95 (2016): 57-66.
#' 
#' Zou, Y., Lord, D., & Zhang, Y. (2012). Analyzing highly dispersed crash data using the Sichel generalized additive models for location, scale and shape. 
#' 
#' 
#' @importFrom Rcpp sourceCpp
#' @useDynLib flexCountReg
#' @export
countreg <- function(formula, data, family = "NB2", offset = NULL, weights = NULL, 
                      verbose = FALSE, dis_param_formula_1 = NULL, 
                      dis_param_formula_2 = NULL, underreport_formula = NULL,
                     underreport_family = "logit",
                     ndraws = 1500, method = "NM", 
                      max.iters = 1000, start.vals = NULL, 
                     stderr = "normal", bootstraps = NULL) {
  
  if (verbose){
    print.level=2
  }else{
    print.level=0
  }

  family <- toupper(str_replace_all(family, "[^[:alnum:]]", "")) # remove non-alphanumeric characters from the family name and ensure all upper case
  method <- toupper(str_replace_all(method, "[^abcfghmnrsABCFGHMNRS]", "")) # clean the method name
  if (!(method %in% c("SANN", "NM", "BFGS", "BFGSR", "CG", "NR", "BHHH"))){
    print('Method must be one of: "SANN", "NM", "BFGS", "BFGSR", "CG", "NR", or "BHHH". Switching to "NM.')
    method = "NM"
  }
  
  # Get the parameters and probability function
  params <- get_params(family)
  probFunc <- get_probFunc(family)
  
  # Prepare model matrices
  mod1_frame <- stats::model.frame(formula, data)
  X <- as.matrix(modelr::model_matrix(data, formula))
  x_names <- colnames(X)
  
  # If underreporting model, create a data matrix for the underreporting model
  if (!is.null(underreport_formula)){
    X_underreport <- as.matrix(modelr::model_matrix(data, underreport_formula))
    N_underreport <- ncol(X_underreport)
    underreport_names <- paste0("Underreporting:", colnames(X_underreport))
  }
  else{
    N_underreport <- 0
    underreport_names <- NULL
  }
  
  # If an offset is specified, create a vector for the offset
  if (!is.null(offset)){
    X_offset <- data[, offset]
  }
  
  if (!is.null(dis_param_formula_1)) {
    mod_alpha_frame <- as.matrix(modelr::model_matrix(data, dis_param_formula_1))
    x_names <- c(x_names, paste0(params[1], ":",colnames(mod_alpha_frame)))
  } else {
    x_names <- append(x_names, params[1])
  }
  
  if (!is.null(dis_param_formula_2)) {
    mod_sigma_frame <- as.matrix(modelr::model_matrix(data, dis_param_formula_2))
    x_names <- c(x_names, paste0(params[2], ":",colnames(mod_sigma_frame)))
  } else {
    x_names <- append(x_names, params[2])
  }
  
  if(!is.null(underreport_formula)){
    x_names <- c(x_names, underreport_names)
  }
  
  y <- stats::model.response(mod1_frame)
  
  # Generate Halton draws to use as quantile values
  haltons <- randtoolbox::halton(ndraws)
  normed_haltons <- dnorm(haltons)
  
  # Some numbers that will be useful
  N_predictors <- ncol(X)
  # N_obs <- length(y)
  N_params <- length(Filter(Negate(is.null), params)) # No. of params for the distribution
  
  if (is.null(dis_param_formula_1)){
    N_alpha = 1
  }
  else{
    N_alpha <- ncol(mod_alpha_frame)
  }
  
  if (N_params==2){
    if(is.null(dis_param_formula_2)){
      if (is.null(params[2])){
        N_sigma = 0
      }else{
        N_sigma = 1
      }
    }
    else{
      N_sigma <- ncol(mod_sigma_frame)
    }
  }else{
    N_sigma = 0
  }
  
  # Initialize starting values if not provided
  start <- if (is.null(start.vals)) {
    
    # Use the NB2 from MASS as starting values
    p_model <- glm.nb(formula, data = data)
    start <- unlist(p_model$coefficients)
    
    if (family=="GW" | family=="SI"){
      start <- c(start, rep(0.1, N_alpha + N_sigma + N_underreport))
    }else{
      start <- c(start, rep(0, N_alpha + N_sigma + N_underreport))
    }
    
  } 
  else {
    start <- start.vals
  }
  
  x_names <- x_names  %>% 
    compact() #%>%  # Remove NULL values
  
  names(start) <-x_names[1:length(start)] # Shouldn't need to do this, but this ensures it runs without having errors from names

  # Handling Weights
  if (is.null(weights)){
    weights.df <- rep(1, length(y))
  }else{
    weights.df <- data[,weights]
  }
  
  # Define the main function for computing log-likelihood
  logLikFunc <- function(p) {
    coefs <- as.array(p)
    fixed_coefs <- as.vector(head(coefs, N_predictors))
    
    if (!is.null(dis_param_formula_1)){
      alpha_coefs <- as.vector(coefs[(N_predictors + 1):(N_predictors + N_alpha)])
      alpha <- exp(mod_alpha_frame %*% alpha_coefs)
    } else {
      alpha <- exp(coefs[(N_predictors + 1)])
    }
    
    if (!is.null(dis_param_formula_2)){
      sigma_coefs <- as.vector(coefs[(N_predictors + N_alpha + 1):(N_predictors + N_alpha + N_sigma)])
      sigma <- exp(mod_sigma_frame %*% sigma_coefs)
    } else {
      sigma <- exp(coefs[(N_predictors + N_alpha + 1)])
    }
    
    if (N_underreport>0){
      beta_underreport <- coefs[(length(coefs)-N_underreport+1):(length(coefs))]
      lin_underreport <- X_underreport %*% beta_underreport
      
      if (underreport_family == "logit"){
        underreport_prob <- 1/(1+exp(-lin_underreport))
      }else{
        underreport_prob <- pnorm(lin_underreport, lower.tail = FALSE)
        
      }
    }else{underreport_prob=1} # If no underreporting model, set the probability to 1
    
    if (!is.null(offset)){
      predicted <- exp(X %*% fixed_coefs + X_offset)*underreport_prob
    } else {
      predicted <- exp(X %*% fixed_coefs)*underreport_prob
    }
    
    # Call the probability function
    probs <- probFunc(y, predicted, alpha, sigma, haltons, normed_haltons)
    
    return(log(probs)*weights.df)
  }
  # Run the maximum likelihood estimation
  fit <- maxLik::maxLik(logLikFunc, 
                        start = start,
                        method = method, 
                        control = list(iterlim = max.iters, 
                                       printLevel = print.level))
  
  # Optionally, compute bootstrapped standard errors
  
  if (stderr=="boot" & is.numeric(bootstraps)) {
    bs.data <- modelr::bootstrap(data, n = bootstraps)
    models <- map(bs.data$strap, ~ mod.boot(data = .))
    tidied <- map_df(models, broom::tidy, .id = "id")
    
    SE <- tidied %>%
      group_by(term) %>%
      summarise(sd = sd(estimate)) %>% deframe()
    
    fit$bootstrapped_se <- SE
  }
  if (stderr == "normal") {
    fit$bootstrapped_se <- sqrt(diag(-1/(fit$hessian)))
  }
  else if (stderr == "bootstrapped") {
    fit$bootstrapped_se <- SE
  }
  else if (stderr == "Robust"){
    fit$bootstrapped_se <- sqrt(diag(sandwich::sandwich(fit)))
  }
  
  fit$coefficients = fit$estimate
  fit$se = if (!is.null(bootstraps) & is.numeric(bootstraps)) fit$bootstrapped_se else sqrt(diag(-1/(fit$hessian)))
  fit$logLik = fit$maximum
  fit$converged = fit$convergence
  fit$model = family
  fit$family = family
  fit$method = method
  fit$data = data
  fit$formula = formula
  fit$dis_param_formula_1 = dis_param_formula_1
  fit$dis_param_formula_2 = dis_param_formula_2
  fit$ndraws = ndraws
  fit$bootstraps = if (!is.null(bootstraps)) bootstraps else NULL
  fit$offset <- offset
  fit.stderr <- stderr
  fit$modelType <- "countreg"
  fit$underreport_formula <- underreport_formula
  fit$underreport_family <- underreport_family
  fit$params <- params
  
  obj = .createFlexCountReg(model = fit, data = data, call = match.call(), formula = formula)
  return(obj)
}

