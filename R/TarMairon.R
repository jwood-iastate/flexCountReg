#' Count regression models With or Without Random Parameters
#' 
#' @name isildursBane
#' @aliases isildursBane
#' @description
#' This function estimates count regression models for a wide range of count 
#' distributions. It allows for estimation of underreporting, random parameters, 
#' random parameters with heterogeneity in means and variances, 
#' panel-specifications, random parameters in the underreporting model, and 
#' modeling the distribution parameters as functions of predictors. It also has 
#' options for handling sampling weights, specifying the optimization method to 
#' use, options for correlated random parameters, and using scrambled Halton 
#' draws for estimating random parameters (which removes correlation in the 
#' pseudo-random Halton draws). The count distributions cover a wide range of 
#' behaviors that capture overdispersion, underdispersion, equidispersion, heavy 
#' tails, heavy densities at 0 (i.e., "excess zeros"), etc.
#' 
#' Three distributions for counts with excess zeros,
#' Five that are complex - requiring numerical approximation,
#' Nine possible formulas to capture complex models,
#' One function for the R package focused on count regression;
#' One function to rule them all,
#' One function to find them,
#' One function to bring them all,
#' And in the estimates bind them;
#' In the flexCountReg package, where count models are found.
#' 
#' 
#' @param formula a symbolic description of the model to be fitted.
#' @param data a data frame containing the variables in the model.
#' @param family the name of the distribution/model type to estimate. The 
#'        default "NB2" is the standard negative binomial distribution with a 
#'        log link. other options are listed below. 
#' @param offsets the name of a variable or a vector of names of variables in 
#'        the dataframe that should be used as an offset (i.e., included but 
#'        forced to have a coefficient of 1).
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
#' @param rpar_obs_formula an optional R formula with the independent variables 
#'        that are observation-specific random parameters. This should not 
#'        include an outcome variable. If the intercept is random (and varies 
#'        over all observations), include it in this formula. If the intercept 
#'        is fixed, include it in `formula` and not in `rpar_obs_formula`. To 
#'        remove the intercept, use \code{0 + vars} or \code{-1 + vars}.
#' @param rpar_panel_formula an optional R formula with the independent 
#'        variables that are panel-specific random parameters. This should not 
#'        include an outcome variable. If the intercept is random (and varies 
#'        across the panels), include it in this formula. Only use this if the 
#'        `panelID` is also specified.
#' @param rpar_underreport_formula an optional R formula with the independent    
#'        variables that are random parameters in the underreporting model. This   
#'        should not include an outcome variable. If the intercept for the 
#'        underreporting model is random include it in this formula and remove 
#'        it from `underreport_formula`. Only use this if the 
#'        `underreport_formula` is also specified.
#' @param het_in_means_formula an optional R formula with independent variables 
#'        that are used to model the heterogeneity in means for the random 
#'        parameters specified in `rpar_obs_formula` and `rpar_panel_formula`. 
#' @param het_in_var_formula an optional R formula with independent variables 
#'        that are used to model the heterogeneity in variance for the random 
#'        parameters specified in `rpar_obs_formula` and `rpar_panel_formula`.
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
#' @details
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
#' @section Random Parameters Model Details:
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
#' **Heterogeiety in Means and Variances**
#' When the heterogeneity in means and/or heterogeneity in variances are 
#' specified using the `het_in_means_formula` and `het_in_var_formula` options, 
#' the mean and standard deviations for each observation are scaled based on 
#' these formulas. The mean for random parameter \eqn{j} for observation \eqn{i} 
#' becomes:
#' \deqn{\beta_{random, \mu, i, j}=e^{X_{het \ in \ means}\beta_{het \ in \ means}}\beta_{random, \mu, j}}
#' 
#' Similarly, the standard deviation for for random parameter \eqn{j} for 
#' observation \eqn{i} becomes:
#' \deqn{\beta_{random, \sigma, i, j}=e^{X_{het \ in \ var}\beta_{het \ in \ var}}\beta_{random, \sigma, j}}
#' 
#' @section Count Distributions and Underreporting Models:  
#' The NB-1, NB-2, and NB-P versions of the negative binomial distribution are 
#' based on Greene (2008).  The details of each of these are provided below.
#' 
#' **NB-1 Model**
#' The PMF and log-likelihood functions are:
#' \deqn{P(Y = y) = \frac{\Gamma(y + \frac{\mu}{\alpha})}{y! \, \Gamma(\frac{\mu}{\alpha})} \left( \frac{\frac{\mu}{\alpha}}{\frac{\mu}{\alpha} + \mu} \right)^{\frac{\mu}{\alpha}} \left( \frac{\mu}{\frac{\mu}{\alpha} + \mu} \right)^y}
#' \deqn{LL_{\text{NB1}}(\beta, \alpha) = \sum_{i=1}^n \left[ \ln \Gamma\left( y_i + \frac{\mu_i}{\alpha} \right) - \ln \Gamma\left( \frac{\mu_i}{\alpha} \right) - \ln y_i! + \frac{\mu_i}{\alpha} \ln \left( \frac{\frac{\mu_i}{\alpha}}{\frac{\mu_i}{\alpha} + \mu_i} \right) + y_i \ln \left( \frac{\mu_i}{\frac{\mu_i}{\alpha} + \mu_i} \right) \right]}
#' 
#' The mean is:
#' \deqn{\mu = exp(X\beta)}
#' 
#' The variance is:
#' \deqn{\text{Var}(Y) = \mu + \alpha\mu}
#' 
#' **NB-2 Model**
#' The PMF and log-likelihood functions are:
#' \deqn{P(Y = y) = \frac{\Gamma(y + \alpha)}{y! \, \Gamma(\alpha)} \left( \frac{\alpha}{\alpha + \mu} \right)^\alpha \left( \frac{\mu}{\alpha + \mu} \right)^y}
#' \deqn{LL_{\text{NB2}} = \sum_{i=1}^n \left[ \ln \Gamma(y_i + \alpha) - \ln \Gamma(\alpha) - \ln y_i! + \alpha \ln \left( \frac{\alpha}{\alpha + \mu_i} \right) + y_i \ln \left( \frac{\mu_i}{\alpha + \mu_i} \right) \right]}
#' 
#' The mean is:
#' \deqn{\mu = exp(X\beta)}
#' 
#' The variance is:
#' \deqn{\text{Var}(Y) = \mu + \alpha\mu^2}
#' 
#' **NB-P Model**
#' The PMF and log-likelihood functions are:
#' \deqn{P(Y = y) = \frac{\Gamma(y + \frac{\mu^{2-p}}{\alpha})}{y! \, \Gamma(\frac{\mu^{2-p}}{\alpha})} \left( \frac{\frac{\mu^{2-p}}{\alpha}}{\frac{\mu^{2-p}}{\alpha} + \mu} \right)^{\frac{\mu^{2-p}}{\alpha}} \left( \frac{\mu}{\frac{\mu^{2-p}}{\alpha} + \mu} \right)^y}
#' \deqn{LL_{\text{NBP}}(\beta, \alpha, p) = \sum_{i=1}^n \left[ \ln \Gamma\left( y_i + \frac{\mu_i^{2-p}}{\alpha} \right) - \ln \Gamma\left( \frac{\mu_i^{2-p}}{\alpha} \right) - \ln y_i! + \frac{\mu_i^{2-p}}{\alpha} \ln \left( \frac{\frac{\mu_i^{2-p}}{\alpha}}{\frac{\mu_i^{2-p}}{\alpha} + \mu_i} \right) + y_i \ln \left( \frac{\mu_i}{\frac{\mu_i^{2-p}}{\alpha} + \mu_i} \right) \right]}
#' 
#' The mean is:
#' \deqn{\mu = exp(X\beta)}
#' 
#' The variance is:
#' \deqn{\text{Var}(Y) = \mu + \alpha\mu^P}
#' 
#' **Poisson-Lognormal (PLN) Model**
#' The compound Probability Mass Function(PMF) for the Poisson-Lognormal distribution is:
#' \deqn{f(y|\lambda,\sigma)=\int_0^\infty \frac{\lambda^y x^y e^{-\lambda x}}{y!}\frac{exp\left(-\frac{ln^2(x)}{2\sigma^2} \right)}{x\sigma\sqrt{2\pi}}dx}
#'
#' Where \eqn{\sigma} is a parameter for the lognormal distribution with the restriction \eqn{\sigma>0}, and \eqn{y} is a non-negative integer.
#'
#' The expected value of the distribution is:
#' \deqn{E[y]=e^{X\beta+\sigma^2/2} = \mu e^{\sigma^2/2}}
#' 
#' When `ln.sigma.formula` is used, the parameter \eqn{\sigma} is modeled as:
#' \deqn{ln(\sigma)=\beta_0+\beta_1 x_1 + \cdots + \beta_n x_n}
#' 
#' Thus, the resulting value for the parameter \eqn{\sigma} is:
#' \deqn{\sigma=e^{\beta_0+\beta_1 x_1 + \cdots + \beta_n x_n}}
#' 
#' The t-statistics and p-values for the coefficients related to ln(sigma) are,
#' by default, testing if the coefficients are different from a value of 0. This
#' has little practical meaning given that they are coefficients for ln(sigma).
#' They are not testing if the coefficients have statistical significance in
#' terms of improvement over a Poisson model. The Likelihood-Ratio test results
#' provided in the output provide a test comparing if the Poisson-Lognormal
#' model provides a statistically significant improvement in model fit over the
#' Poisson model.
#' 
#' **Poisson Generalized-Exponential (PGE) Model**
#' #' The Generalized Exponential distribution can be written as a function with a
#' shape parameter \eqn{\alpha>0} and scale parameter \eqn{\gamma>0}. The
#' distribution has strictly positive continuous values. The PDF of the
#' distribution is:
#' \deqn{f(x|\alpha,\gamma)=\frac{\alpha}{\gamma}\left(1-e^{-\frac{x}{\gamma}}\right)^{\alpha-1}e^{-\frac{x}{\gamma}}} 
#' 
#' Thus, the compound Probability Mass Function(PMF) for the PGE distribution
#' is:
#' \deqn{f(y|\lambda,\alpha,\beta)=\int_0^\infty \frac{\lambda^y x^y e^{-\lambda x}}{y!}\frac{\alpha}{\gamma}\left(1-e^{-\frac{x}{\gamma}}\right)^{\alpha-1}e^{-\frac{x}{\gamma}} dx}
#' 
#' The expected value of the distribution is:
#' \deqn{E[y]=\mu=\lambda \left(\frac{\psi(\alpha+1)-\psi(1)}{\gamma}\right)}
#' 
#' Where \eqn{\psi(\cdot)} is the digamma function.
#' 
#' The variance is:
#' \deqn{\sigma^2=\lambda \left(\frac{\psi(\alpha+1)-\psi(1)}{\gamma}\right) + \left(\frac{-\psi'(\alpha+1)+\psi'(1)}{\gamma^2}\right)\lambda^2}
#' 
#' Where \eqn{\psi'(\cdot)} is the trigamma function.
#' 
#' To ensure that \eqn{\mu=e^{X\beta}}, \eqn{\lambda} is replaced with:
#' \deqn{\lambda=\frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)}}
#' 
#' This results in:
#' \deqn{f(y|\mu,\alpha,\beta)=\int_0^\infty \frac{\left(\frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)}\right)^y x^y e^{-\left(\frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)}\right) x}}{y!}\frac{\alpha}{\gamma}\left(1-e^{-\frac{x}{\gamma}}\right)^{\alpha-1}e^{-\frac{x}{\gamma}} dx}
#' 
#' Halton draws are used to perform simulation over the lognormal distribution to solve the integral.
#' 
#' **Poisson-Inverse-Gaussian Type 1 (PIG1)  and Type 2 (PIG2) Models**
#' The Poisson-Inverse-Gaussian regression model is based on the Poisson-Inverse-Gaussian Distribution. 
#' 
#' The expected value of the distribution in the regression utilizes a log-link function. Thus, the mean is:
#' \deqn{\mu=e^{X\beta}}
#'
#' The variance function for the Type 1 distribution (which is the default) is:
#' \deqn{\sigma^2=\mu+\eta\mu}
#' 
#' While the variance for the Type 2 distribution is:
#' \deqn{\sigma^2=\mu+\eta\mu^2}
#' 
#' The parameter \eqn{\eta} is estimated as the natural logarithm transformed value, \eqn{\ln(\eta)}, to ensure that \eqn{\eta>0}.
#' 
#' **Poisson-Lindley (PL) Model**
#' The Poisson-Lindley regression is based on a compound Poisson-Lindley 
#' distribution. It handles count outcomes with high levels of zero 
#' observations (or other high densities at low outcome values) that standard 
#' count regression methods, including the negative binomial, may struggle to 
#' adequately capture or model.
#'
#' The compound Probability Mass Function(PMF) for the Poisson-Lindley (PL) 
#' distribution is:
#' \deqn{f(y|\theta,\lambda)=\frac{\theta^2\lambda^y(\theta+\lambda+ 
#' y+1)}{(\theta+1)(\theta+\lambda)^{y+2}}}
#'
#' Where \eqn{\theta} and \eqn{\lambda} are distribution parameters with the 
#' restrictions that \eqn{\theta>0} and \eqn{\lambda>0}, and \eqn{y} is a 
#' non-negative integer.
#'
#' The expected value of the distribution is:
#' \deqn{\mu=\frac{\lambda(\theta+2)}{\theta(\theta+1)}}
#'
#' If a log-link function is used, the mean is:
#' \deqn{\mu=e^{X\beta}=\frac{\lambda(\theta+2)}{\theta(\theta+1)}}
#'
#' Thus, the parameter \eqn{\lambda} in the PL distribution when applied to 
#' regression analysis is:
#' \deqn{\lambda=\frac{\mu\theta(\theta+1)}{\theta+2}}
#' 
#' Using the replacement and simplifying results in:
#' \deqn{f(y \mid \theta, \mu) = \\ \frac{\theta^2 (\mu \theta (\theta+1))^y 
#' (\theta^2 (1+\mu) + \theta (2+\mu) + (\theta+2) (y+1))}{(\theta+1) 
#' (\theta+2)^{y+1} (\theta^2 (1+\mu) + \theta (2+\mu))^{y+2}}}
#' And
#' \deqn{LL=2 \log(\theta) + y (\log(\mu) + \log(\theta) + \log(\theta+1)) + 
#' \log(\theta^2 (1+\mu) + \theta (2+\mu) + (\theta+2) (y+1)) - \log(\theta+1) 
#' - (y+1) \log(\theta+2) - (y+2) \log(\theta^2 (1+\mu) + \theta (2+\mu))}
#'
#' The variance function is defined as:
#' \deqn{\sigma^2=\mu+\left(1-\frac{2}{(\theta+2)^2}\right)\mu^2}
#' 
#' It should be noted that the p-value for the parameter `ln(theta)` in the 
#' model summary is testing if the parameter `theta` is equal to a value of 1. 
#' This has no practical meaning. The Likelihood-Ratio (LR) test compares the 
#' Poisson-Lindley regression with a Poisson regression with the same 
#' independent variables. Thus, the PR test result indicates the statistical 
#' significance for the improvement in how well the model fits the data over a 
#' Poisson regression. This indicates the statistical significance of the 
#' `theta` parameter.
#' 
#' **Poisson-Lindley-Gamma (PLG) Model**
#' The Poisson-Lindley-Gamma regression is based on a compound Poisson-Lindley-
#' Gamma distribution. Details of the distribution can be seen at 
#' \code{\link[flexCountReg]{dplindGamma}}.
#'
#' The mean for the regression model is:
#' \deqn{\mu=e^{X\beta}}
#'
#' The variance function is defined as:
#' \deqn{\sigma^2=\mu+\left(\alpha+1-\frac{2}{(\theta+2)^2}\right)\mu^2}
#'
#' It should be noted that the p-value for the parameters `ln(theta)` and 
#' `ln(alpha)` in the model summary are testing if the parameter `theta` and 
#' `alpha` are equal to a value of 1.
#' 
#' **Poisson-Lindley-Lognormal (PLL) Model**
#' The Poisson-Lindley-Lognormal regression is based on a compound Poisson-
#' Lindley-Lognormal distribution. Details of the distribution can be seen at 
#' \code{\link[flexCountReg]{dplindLnorm}}.
#'
#' The mean for the regression model is:
#' \deqn{\mu=e^{X\beta}}
#'
#' The variance function is defined as:
#' \deqn{\sigma^2=\mu+\left(\frac{1-\frac{2}{(\theta+2)^2}}{e^{\frac{\sigma^2}{2}}}+e^{\sigma^2}-1\right)\mu^2}
#'
#' It should be noted that the p-value for the parameters `ln(theta)` and 
#' `ln(sigma)` in the model summary are testing if the parameter `theta` and 
#' `sigma` are equal to a value of 1.
#'
#' **Poisson-Weibull (PW) Model**
#' The Poisson-Weibull distribution uses the Weibull distribution as a mixing 
#' distribution for a Poisson process. It is useful for modeling overdispersed 
#' count data. The density function (probability mass function) for the 
#' Poisson-Weibull distribution is given by:
#' \deqn{P(y|\lambda,\alpha,\sigma) = \int_0^\infty \frac{e^{-\lambda x} \lambda^y x^y }{y!} \left(\frac{\alpha}{\sigma}\right)\left(\frac{x}{\sigma}\right)^{\alpha-1}e^{-\left(\frac{x}{\sigma}\right)^\alpha} dx}
#' where \eqn{f(x| \alpha, \sigma)} is the PDF of the Weibull distribution and 
#' \eqn{\lambda} is the mean of the Poisson distribution.
#' 
#' For the Poisson-Weibull Regression model, the expected values is:
#' \deqn{E[Y]=\lambda\sigma\Gamma\left(1+\frac{1}{\alpha}\right)}
#' Where \eqn{\lambda} is the mean of the Poisson distribution, \eqn{\alpha} is 
#' the shape parameter, and \eqn{\sigma} is the scale parameter.
#' 
#' To ensure that the regression model predicts the mean value, the regression 
#' utilizes:
#' \deqn{\mu=\exp{X\gamma}=\lambda\sigma\Gamma\left(1+\frac{1}{\alpha}\right)}
#' Where \eqn{X} is a matrix of independent variables and \eqn{\gamma} is a 
#' vector of coefficients. 
#' 
#' This leads to:
#' \deqn{\lambda=\frac{\mu}{\sigma\Gamma\left(1+\frac{1}{\alpha}\right)}}
#' 
#' The variance for the Poisson-Weibull regression is:
#' \deqn{V[Y]=\mu+\left(\frac{\Gamma\left(1+\frac{2}{\alpha}\right)}{\Gamma\left(1+\frac{1}{\alpha}\right)^2}-1\right)\mu^2}
#' 
#' **Sichel (SI) Model**
#' The compound Probability Mass Function (PMF) for the Sichel distribution uses 
#' the formulation from Zhou et al. (2011) and Rigby et al. (2008):
#' \deqn{f(y|\mu, \sigma, \gamma)=\frac{\left(\frac{\mu}{c}\right)^y K_{y+\gamma}(\alpha)}{K_\gamma(1/\sigma)y!(\alpha\sigma)^{y+\gamma}}}
#' 
#' Where \eqn{\sigma} and \eqn{\gamma} are distribution parameters with \eqn{-\infty < \gamma < \infty} and \eqn{\sigma>0}, 
#' \eqn{c=\frac{K_{\gamma+1}(1/\sigma)}{K_\gamma(1/\sigma)}}, \eqn{\alpha^2=\sigma^{-2}+2\mu(c\sigma)^{-1}}, 
#' a mean value of \eqn{\mu}, \eqn{y} is a non-negative integer, and 
#' \eqn{K_j(x)} is a modified Bessel function of the third kind with order 
#' \eqn{j} and argument \eqn{x}.
#'
#' The variance of the distribution is:
#' \deqn{\sigma^2=\mu+\left(\frac{2\sigma(\gamma+1)}{c}+\frac{1}{c^2}-1\right)\mu^2}
#' 
#' **Generalized Waring (GW) Model**
#' The following are the versions of the PMF, mean, and variance used for the 
#' Generaalized Waring model. This is adjusted from the typical formulation by 
#' replacing parameter \code{k} with \eqn{\mu}
#' \deqn{PMF=\frac{\Gamma(\alpha+y)\Gamma(k+y)\Gamma(\rho+k)\Gamma(\alpha+\rho)}{y!\Gamma(\alpha)\Gamma(k)\Gamma(\rho)\Gamma(\alpha+k+\rho+y)}}
#' \deqn{\mu=e^{X\beta}=\frac{\alpha k}{\rho-1}}
#' \deqn{\sigma^2=\frac{\alpha k(\alpha+k+\rho-1)}{(\rho-1)^2(\rho-2)}}
#' 
#' The distribution parameters are often considered to capture the randomness 
#' (parameter \deqn{\alpha}), proneness (parameter \deqn{}k), and liability 
#' (parameter \deqn{\rho}) of the data.
#'
#' #' If we use:
#' \deqn{\alpha=\frac{\mu k}{\rho-1}}
#' 
#' The PMF becomes:
#' 
#' \deqn{PMF=\frac{\Gamma\left(\frac{\mu k}{\rho-1}+y\right)\Gamma(k+y)\Gamma(\rho+k)\Gamma\left(\frac{\mu k}{\rho-1}+\rho\right)}{y!\Gamma\left(\frac{\mu k}{\rho-1}\right)\Gamma(k)\Gamma(\rho)\Gamma\left(\frac{\mu k}{\rho-1}+k+\rho+y\right)}}
#' 
#' #' This results in a regression model where:
#' \deqn{\mu=e^{X\beta}}
#' \deqn{\sigma^2=\frac{\mu k^2\left(\frac{\mu k}{\rho-1}+k+\rho-1\right)}{(\rho-1)^3(\rho-2)}=\left(\frac{k^3+\rho k^2- k^2}{(\rho-1)^3(\rho-2)}\right)\mu+\left(\frac{k^3}{(\rho-1)^4(\rho-2)}\right)\mu^2}
#'
#' Note that when \deqn{p=1} or \deqn{p=2}, the distribution is undefined.
#' 
#' **Conway-Maxwell-Poisson (COM) Model**
#' The following is the the PMF for the COM model.
#' \deqn{f(x|\lambda, \nu)=\frac{\lambda^x}{(x!)^\nu Z(\lambda,\nu)}}
#' 
#' Where \eqn{\lambda} and \eqn{\nu} are distribution parameters with \eqn{\lambda>0} and \eqn{\nu>0}, and \eqn{Z(\lambda,\nu)} is the normalizing constant.
#' 
#' The normalizing constant is given by:
#' \deqn{Z(\lambda,\nu)=\sum_{n=0}^{\infty}\frac{\lambda^n}{(n!)^\nu}}
#' 
#' The mean and variance are:
#' \deqn{\mu=e^{X\beta}=\lambda \frac{\delta}{\delta \lambda} \log(Z(\lambda,\nu))}
#' \deqn{\sigma^2=\lambda \frac{\delta}{\delta \lambda} \mu}
#' 
#' Note that the COM distribution parameter \eqn{\lambda} is solved for using \eqn{\mu} and \eqn{\nu}, so the regression model provides direct predictions for the mean.
#' 
#' **Underreporting**
#' Models for underreporting combine a binary probability model (logit or 
#' probit) with a count model. This is accomplished using a model for the 
#' probability of crashes being reported multiplied by the estimated mean for 
#' the count model, based on the observed data. This is discussed in Wood et. 
#' al. (2016), Pararai et. al., (2006), and Pararai et. al., (2010). The 
#' underreporting model is based on:
#' \deqn{\mu_{true}=\mu_{observed}\cdot P(\text{event is reported})}
#' 
#' This allows the inference of both the true event count and the probability of
#' the event being reported as a function of independent variables. 
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
#' @import modelr randtoolbox
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom purrr map map_df compact
#' @importFrom broom tidy
#' @importFrom dplyr group_by %>% reframe
#' @importFrom tibble deframe
#' @importFrom maxLik maxLik
#' @importFrom stringr str_replace_all
#' @importFrom sandwich sandwich
#' @include pinvgaus.R ppoislogn.R plindLnorm.R plindGamma.R psichel.R Generalized-Waring.R ppoisGE.R psichel.R plind.R helpers.R RP_Model_Helper_Functions.R
#' 
#' @examples
#' # Load the Washington data
#' data("washington_roads")
#' washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
#' 
#' # Estimate an NB2 model with a dispersion parameter as a function of the 
#' # variable `speed50` (i.e., generalized NB2), verbose output, and use the 
#' # BFGS optimization method. No random parameters.
#' nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'                data = washington_roads, family = "NB2",
#'                dis_param_formula_1 = ~ speed50, verbose = TRUE, method='BFGS')
#' summary(nb2)
#' 
#' # Estimate a Poisson-Lognormal model (a low number of draws is used to speed 
#' # up the estimation for examples - not recommended in practice). No random 
#' # parameters.
#' pln <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'               data = washington_roads, family = "PLN", ndraws=10)
#' summary(pln)  
#' 
#' # Estimate an Poisson-Lognormal with underreporting (probit)
#' # No random parameters.
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
isildursBane <- function(formula, 
                     data, 
                     family = "NB2", 
                     offsets = NULL, 
                     weights = NULL, 
                     verbose = FALSE, 
                     dis_param_formula_1 = NULL, 
                     dis_param_formula_2 = NULL, 
                     underreport_formula = NULL,
                     underreport_family = "logit",
                     rpar_obs_formula = NULL,
                     rpar_panel_formula = NULL,
                     rpar_underreport_formula = NULL,
                     het_in_means_formula = NULL,
                     het_in_var_formula = NULL,
                     rpar_obs_dists = NULL,
                     rpar_panel_dists = NULL,
                     rpar_underreport_dists = NULL,
                     correlated = FALSE,
                     ndraws = 1500,
                     scrambled = TRUE,
                     panelID = NULL,  
                     method = "NM", 
                     max.iters = 1000, 
                     start.vals = NULL, 
                     stderr = "normal", 
                     bootstraps = NULL) {
  
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
  
  # Prep Data for Modeling
  prepped_data <- data_prep(formula, data, rpar_obs_formula,  rpar_panel_formula,  
                            underreport_formula, dis_param_formula_1, dis_param_formula_2,
                            rpar_underreport_formula, het_in_means_formula, 
                            het_in_var_formula, offsets, panelID)
  
  y                       <- prepped_data$y
  X_fixed                 <-prepped_data$X_fixed 
  X_rand_crosssectional   <-prepped_data$X_rand_crosssectional
  X_rand_panel            <-prepped_data$X_rand_panel 
  X_underreport           <-prepped_data$X_underreport
  X_param1                <-prepped_data$X_param1
  X_param2                <-prepped_data$X_param2
  X_underreport_rpar      <-prepped_data$X_underreport_rpar
  X_het_means             <-prepped_data$X_het_means
  X_het_variance          <-prepped_data$X_het_variance
  X_offsets               <-prepped_data$X_offsets
  
  # Get the parameters and probability function
  params <- get_params(family)
  probFunc <- get_probFunc(family)
  
  
  # Get number of parameters for model elements
  N_list <- get_Nparams <- function(prepped_data, family)

  N_fixed               <- N_list$N_fixed 
  N_rand_crosssectional <- N_list$N_rand_crosssectional 
  N_rand_panel          <- N_list$N_rand_panel
  N_underreport         <- N_list$N_underreport 
  N_param1              <- N_list$N_param1 
  N_param2              <- N_list$N_param2 
  N_underreport_rpar    <- N_list$N_underreport_rpar 
  N_het_means           <- N_list$N_het_means 
  N_het_variance        <- N_list$N_het_variance 
  N_family_params       <- N_list$family_params
  N_rand_total          <- N_list$N_rand_total
  N_Offsets             <- N_list$N_Offsets
  
  # If an offset is specified, create a vector for the offset
  if (!is.null(offset)){
    X_offset <- data[, offset]
  }
  
  if (is.null(start.vals)){
    start <- get_start_vals(formula, prepped_data, data, family,
                            rpar_obs_formula,  rpar_panel_formula, correlated)
  }
  
  # check random parameters. For random parameter equations, if the distributions are not specified, set them to "n". Otherwise, set them to the specified distribution
  if (N_rand_crosssectional>0 & names(X_rand_crosssectional) != names(rpar_obs_dists)){
    rpar_obs_dists <- rep("n", N_rand_crosssectional)
    names(rpar_obs_dists) <- names(X_rand_crosssectional)
  }
  
  if (N_rand_panel>0 & names(X_rand_panel) != names(rpar_panel_dists)){
    rpar_panel_dists <- rep("n", N_rand_panel)
    names(rpar_panel_dists) <- names(X_rand_panel)
  }
  if (N_underreport>0 & names(X_underreport) != names(rpar_underreport_dists)){
    rpar_underreport_dists <- rep("n", N_underreport)
    names(rpar_underreport_dists) <- names(X_underreport)
  }
  
  if(!is.null(rpar_obs_dists) & !is.null(rpar_panel_dists) & !is.null(rpar_underreport_dists)){
    rpar_dists <- c(rpar_obs_dists, rpar_panel_dists, rpar_underreport_dists)
  }else if(!is.null(rpar_obs_dists) & !is.null(rpar_panel_dists)){
    rpar_dists <- c(rpar_obs_dists, rpar_panel_dists)
  }else if(!is.null(rpar_obs_dists) & !is.null(rpar_underreport_dists)){
    rpar_dists <- c(rpar_obs_dists, rpar_underreport_dists)
  }else if(!is.null(rpar_panel_dists) & !is.null(rpar_underreport_dists)){
    rpar_dists <- c(rpar_panel_dists, rpar_underreport_dists)
  }else if(!is.null(rpar_obs_dists)){
    rpar_dists <- rpar_obs_dists
  }else if(!is.null(rpar_panel_dists)){
    rpar_dists <- rpar_panel_dists
  }else if(!is.null(rpar_underreport_dists)){
    rpar_dists <- rpar_underreport_dists
  }else{
    stop("No random parameters were specified. Use the rpar_obs_formula, rpar_panel_formula, or rpar_underreport_formula arguments to specify random parameters.")
  }
  
  # Generate halton draws for random parameters
  rpar_haltons <- get_halton_draws(ndraws, dim=(N_rand_total+1), scrambled=scrambled)
  
  # Seperate draws for use later
  rpardraws <- rpar_haltons$rpardraws
  distdraws <- rpar_haltons$distdraws
  normed_distdraws <- stats::pnorm(distdraws)
  
  ##############################################################################
  ## The remainder needs to be revised to fit all of the models ################
  ##############################################################################
  
  
  # Handling Weights
  if (is.null(weights)){
    weights.df <- rep(1, length(y))
  }else{
    weights.df <- data[,weights]
  }

  # Define the main function for computing log-likelihood
  logLikFunc <- function(p) {
    coefs <- as.array(p)
    coef_vals <- coef_vals(coefs, formula, prepped_data, data,
                          family=family, correlated=correlated)
    
    beta_fixed = coef_vals$beta_fixed
    beta_underreport = coef_vals$beta_underreport
    beta_crss_mean = coef_vals$beta_crss_mean
    beta_panel_mean = coef_vals$beta_panel_mean
    beta_underreport_rpar = coef_vals$beta_underreport_rpar
    beta_crss_std = coef_vals$beta_crss_std
    beta_panel_std = coef_vals$beta_panel_std
    beta_underreport_rpar_std =  coef_vals$beta_underreport_rpar_std
    beta_het_means =  coef_vals$beta_het_means
    beta_het_variance = coef_vals$beta_het_variance
    beta_param1 = coef_vals$beta_param1
    beta_param2 = coef_vals$beta_param2
    beta_chol = coef_vals$beta_chol
    
    
    mu_fixed  <- exp(X_fixed %*% beta_fixed)
    
    if (!is.null(dis_param_formula_1)){
      alpha <- exp(X_param1 %*% beta_param1)
    } else {
      alpha <- exp(beta_param1)
    }
    
    if (!is.null(dis_param_formula_2)){
      sigma <- exp(X_param2 %*% beta_param2)
    } else {
      sigma <- exp(beta_param2)
    }
    
    # set up random parameters
    beta_rand_means = c()
    beta_rand_sd = c()
    if(N_rand_crosssectional>0){
      beta_rand_means = c(beta_rand_means, beta_crss_mean)
      beta_rand_sd = c(beta_rand_sd, beta_crss_std)
    }
    if(N_rand_panel>0){
      beta_rand_means = c(beta_rand_means, beta_panel_mean)
      beta_rand_sd = c(beta_rand_sd, beta_panel_std)
    }
    if(N_underreport>0){
      beta_rand_means = c(beta_rand_means, beta_underreport)
      beta_rand_sd = c(beta_rand_sd, beta_underreport_rpar_std)
    }
    
    # Computations for heterogeneity in means and variances
    if(N_het_means>0){ # If there is heterogeneity in means, the observation specific means are the mean coefficients multiplied by the X\beta values
      XB_het_means = exp(X_het_means %*% beta_het_means)
      # Create a matrix of the means for each observation for each random parameter by multiplying the mean coefficients by the X\beta values
      rp_means_it = X_het_means %outer% beta_rand_means
    }else{
      rp_means_it = beta_rand_means
    }
    
    if(N_het_variance>0){
      beta_rand_sd = exp(X_het_variance  %*% beta_het_variance)
      rp_sd_it = X_het_variance %outer% beta_rand_sd
    }else{
      rp_sd_it = beta_rand_sd
    }
    
    # If no heterogeneity in means or variances, get the random parameter draws
    if(N_het_means==0 & N_het_variance==0){
      draws_converted = generate_draws(rpardraws, rp_means_it, rp_sd_it, rpar_dists)
      
      if(N_rand_crosssectional>0){
        draws_crss = draws_converted[,1:N_rand_crosssectional]
      }
      if(N_rand_panel>0){
        draws_panel = draws_converted[,N_rand_crosssectional+1:(N_rand_crosssectional+N_rand_panel)]
      }
      if(N_underreport>0){
        draws_underreport = draws_converted[,(N_rand_crosssectional+N_rand_panel+1):(N_rand_crosssectional+N_rand_panel+N_underreport)]
      }
      
      if (!is.null(underreport_formula)){
        lin_underreport <- X_underreport %*% beta_underreport
        
        if(N_underreport>0){
          lin_underreport <- lin_underreport + X_underreport_rpar %*% beta_underreport_rpar
          
          if (underreport_family == "logit"){
            underreport_prob_i <- 1/(1+exp(-lin_underreport))
          }else{
            underreport_prob_i <- pnorm(lin_underreport, lower.tail = FALSE)
          }
          underreport_prob <- rowMeans(underreport_prob_i)
          
        }else{
          underreport_prob_i <- 1/(1+exp(-lin_underreport))
        }
        mu_fixed <- mu_fixed * underreport_prob
      }
    
    # if there is a panel specification, group by the panel ID to do the computations
      
      # If there is also individual random parameters, treat those random draws seperately
    
    
    } # end of no heterogeneity of means or variances
    

    
    if (N_underreport>0){
      lin_underreport <- X_underreport %*% beta_underreport
      
      if (underreport_family == "logit"){
        # if there are random parameters for underreporting, add them to the linear predictor then compute the probabilities, then take the average for each row
        if (N_underreport_rpar>0){
          
          lin_underreport <- lin_underreport + X_underreport_rpar %*% beta_underreport_rpar
          underreport_prob <- 1/(1+exp(-lin_underreport))
        
        
        underreport_prob <- 1/(1+exp(-lin_underreport))
      }else{
        underreport_prob <- pnorm(lin_underreport, lower.tail = FALSE)
        
      }
    }else{underreport_prob=1} # If no underreporting model, set the probability to 1
    
    # multiply the mu_fixed by the probability
    mu_fixed <- mu_fixed * underreport_prob
    
    if(!is.null(offsets)){ # Correct mu for all offsets
      if (length(offsets)==1){
        mu_fixed <- mu_fixed * exp(data[[offsets]])
      }
      else{
        for (i in offsets){
          mu_fixed <- mu_fixed * exp(rowSums(X_offsets))
        }
      }
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
      reframe(sd = sd(estimate))
    
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

