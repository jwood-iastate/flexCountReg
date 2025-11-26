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
#' @param offset the name of a variable, or vector of variable names, in the 
#'        data frame that should be used as an offset (i.e., included but forced 
#'        to have a coefficient of 1). The normal method of setting an offset in 
#'        the equation can also be used (overrides the offset option).
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
#' @param ndraws The number of Halton draws for integrating the distribution 
#'        being compounded with the Poisson distribution when there is not a 
#'        closed-form solution. Default is 1500. It is recommended to test 
#'        different numbers of draws to determine if the model is stable (i.e., 
#'        doesn't change or has minimal change as the number of draws changes 
#'        within a reasonable range).
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
#'  \item Poisson-Inverse-Gamma (PIG)
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
#' \item "POISSON" for Poisson distribution with a log link.
#'  \item "NB1" for Negative Binomial 1 distribution with a log link.
#'  \item "NB2" for Negative Binomial 2 distribution with a log link (i.e., the 
#'        standard negative binomial model).
#'  \item "NBP" for Negative Binomial P distribution with a log link.
#'  \item "PLN" for Poisson-Lognormal distribution with a log link.
#'  \item "PGE" for Poisson-Generalized-Exponential distribution with a log 
#'        link.
#'  \item "PIG1" for Poisson-Inverse-Gaussian Type-1 distribution with a log link.
#'  \item "PIG2" for Poisson-Inverse-Gaussian Type-2 distribution with a log link.
#'  \item "PIG" for Poisson-Inverse-Gamma distribution with a log link.
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
#'  \item \eqn{ln(\alpha)} for the Negative Binomial 1 model.
#'  \item \eqn{ln(\alpha)} for the Negative Binomial 2 model.
#'  \item \eqn{ln(\alpha)} for the Negative Binomial P model.
#'  \item \eqn{ln(\sigma)} for the Poisson-Lognormal model.
#'  \item shape parameter for the Poisson-Generalized-Exponential model.
#'  \item \eqn{ln(\eta)} for the Poisson-Inverse-Gaussian model.
#'  \item \eqn{ln(\eta)} for the Poisson-Inverse-Gamma model.
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
#'  \item Not Applicable for the Negative Binomial 1 model.
#'  \item Not Applicable for the Negative Binomial 2 model.
#'  \item p for the Negative Binomial P model.
#'  \item Not Applicable for the Poisson-Lognormal model.
#'  \item scale parameter for the Poisson-Generalized-Exponential model.
#'  \item Not Applicable for the Poisson-Inverse-Gaussian model.
#'  \item Not Applicable for the Poisson-Inverse-Gamma model.
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
#'  \item Poisson-Lindley-Gamma (more efficient than using hypergeometric 
#'  functions)
#'  \item Poisson-Lindley-Lognormal
#'  \item Poisson-Weibull
#'  }
#' 
#' @section Model Details: 
#' ## Poisson Model
#' This implements the Poisson regression model using Maximum Liklelihood 
#' Estimation, as opposed to the Iteratively Reweighted Least Squares (IRLS) 
#' method used in the `glm` function.
#' 
#' The PMF and log-likelihood functions are:
#' \deqn{P(Y = y) = \frac{e^{-\mu} \mu^y}{y!}}
#' \deqn{LL_{\text{Poisson}}(\beta) = \sum_{i=1}^n \left[ -\mu_i + y_i \ln(\mu_i) - \ln(y_i!) \right]}
#' 
#' The mean is:
#' \deqn{\mu = exp(X\beta)}
#' 
#' The variance is:
#' \deqn{\text{Var}(Y) = \mu}
#' 
#' ## Negative Binomial Models**
#' 
#' The NB-1, NB-2, and NB-P versions of the negative binomial distribution are 
#' based on Greene (2008).  The details of each of these are provided below.
#' 
#' ### NB-1 Model
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
#' ### NB-2 Model
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
#' ### NB-P Model
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
#' ## Poisson-Lognormal (PLN) Model
#' The compound Probability Mass Function(PMF) for the Poisson-Lognormal 
#' distribution is:
#' \deqn{f(y|\lambda,\sigma)=\int_0^\infty \frac{\lambda^y x^y e^{-\lambda x}}{y!}\frac{exp\left(-\frac{ln^2(x)}{2\sigma^2} \right)}{x\sigma\sqrt{2\pi}}dx}
#'
#' Where \eqn{\sigma} is a parameter for the lognormal distribution with the 
#' restriction \eqn{\sigma>0}, and \eqn{y} is a non-negative integer.
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
#' ## Poisson Generalized-Exponential (PGE) Model
#' The Generalized Exponential distribution can be written as a function with 
#' a shape parameter \eqn{\alpha>0} and scale parameter \eqn{\gamma>0}. The
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
#' Halton draws are used to perform simulation over the lognormal distribution 
#' to solve the integral.
#' 
#' ## Poisson-Inverse-Gaussian Type 1 (PIG1)  and Type 2 (PIG2) Models
#' The Poisson-Inverse-Gaussian regression model is based on the 
#' Poisson-Inverse-Gaussian Distribution. 
#' 
#' The expected value of the distribution in the regression utilizes a log-link 
#' function. Thus, the mean is:
#' \deqn{\mu=e^{X\beta}}
#'
#' The variance function for the Type 1 distribution (which is the default) is:
#' \deqn{\sigma^2=\mu+\eta\mu}
#' 
#' While the variance for the Type 2 distribution is:
#' \deqn{\sigma^2=\mu+\eta\mu^2}
#' 
#' The parameter \eqn{\eta} is estimated as the natural logarithm transformed 
#' value, \eqn{\ln(\eta)}, to ensure that \eqn{\eta>0}.
#' 
#' ## Poisson-Inverse-Gamma (PIG) Model
#' The PDF of the distribution is:
#' \deqn{f(x|\eta,\mu)=\frac{2\left(\mu\left(\frac{1}{\eta}+1\right)\right)^{\frac{x+\frac{1}{eta}+2}{2}}}{x!\Gamma\left(\frac{1}{\eta}+2\right)}K_{x-\frac{1}{\eta}-2}\left(2\sqrt{\mu\left(\frac{1}{\eta}+1\right)}\right)}
#' 
#' Where \eqn{\eta} is a shape parameter with the restriction that \eqn{\eta>0}, 
#' \eqn{\mu>0} is the mean value,  \eqn{y} is a non-negative integer, and 
#' \eqn{K_i(z)} is the modified Bessel function of the second kind. This 
#' formulation uses the mean directly.
#'
#' The variance of the distribution is:
#' \deqn{\sigma^2=\mu+\eta\mu^2}
#' 
#' ## Poisson-Lindley (PL) Model
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
#' ## Poisson-Lindley-Gamma (PLG) Model
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
#' ## Poisson-Lindley-Lognormal (PLL) Model
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
#' ## Poisson-Weibull (PW) Model
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
#' ## Sichel (SI) Model
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
#' ## Generalized Waring (GW) Model
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
#' If we use:
#' \deqn{\alpha=\frac{\mu k}{\rho-1}}
#' 
#' The PMF becomes:
#' 
#' \deqn{PMF=\frac{\Gamma\left(\frac{\mu k}{\rho-1}+y\right)\Gamma(k+y)\Gamma(\rho+k)\Gamma\left(\frac{\mu k}{\rho-1}+\rho\right)}{y!\Gamma\left(\frac{\mu k}{\rho-1}\right)\Gamma(k)\Gamma(\rho)\Gamma\left(\frac{\mu k}{\rho-1}+k+\rho+y\right)}}
#' 
#' This results in a regression model where:
#' \deqn{\mu=e^{X\beta}}
#' \deqn{\sigma^2=\frac{\mu k^2\left(\frac{\mu k}{\rho-1}+k+\rho-1\right)}{(\rho-1)^3(\rho-2)}=\left(\frac{k^3+\rho k^2- k^2}{(\rho-1)^3(\rho-2)}\right)\mu+\left(\frac{k^3}{(\rho-1)^4(\rho-2)}\right)\mu^2}
#'
#' Note that when \deqn{p=1} or \deqn{p=2}, the distribution is undefined.
#' 
#' ## Conway-Maxwell-Poisson (COM) Model
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
#' Note that the COM distribution parameter \eqn{\lambda} is solved for using 
#' \eqn{\mu} and \eqn{\nu}, so the regression model provides direct predictions 
#' for the mean.
#' 
#' ## Underreporting
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
#' @importFrom dplyr mutate %>% row_number group_by across all_of summarize ungroup reframe pull
#' @importFrom tibble deframe
#' @importFrom maxLik maxLik
#' @importFrom stringr str_replace_all str_extract
#' @importFrom sandwich sandwich
#' @include pinvgaus.R pinvgamma.R ppoislogn.R plindLnorm.R plindGamma.R psichel.R Generalized-Waring.R ppoisGE.R psichel.R plind.R helpers.R
#' 
#' @examples
#' \donttest{
#' # Load the Washington data
#' data("washington_roads")
#' washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
#' # Estimate an NB2 model with a dispersion parameter as a function of the 
#' # variable `speed50` (i.e., generalized NB2), verbose output, and use the 
#' # BFGS optimization method
#' nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'                data = washington_roads, family = "NB2",
#'                dis_param_formula_1 = ~ speed50, verbose = TRUE, 
#'                method='BFGS')
#' summary(nb2)
#' 
#' 
#' # Estimate a Poisson-Lognormal model (a low number of draws is used to speed 
#' # up the estimation for examples - not recommended in practice)
#' pln <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'               data = washington_roads, family = "PLN", ndraws=10)
#' summary(pln)  
#' 
#' # Estimate an Poisson-Lognormal with underreporting (probit)
#' plogn_underreport <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + 
#'                AADT10kplus,
#'               data = washington_roads, family = "NB2",
#'               underreport_formula = ~ speed50 + AADT10kplus, 
#'               underreport_family = "probit")
#' summary(plogn_underreport)
#' 
#' # Estimate a Conway-Maxwell-Poisson model
#' com_model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + 
#'                       AADT10kplus,
#'                       data = washington_roads, family = "COM", method="BHHH")
#' summary(com_model)
#' #}
#' 
#' @references
#' Greene, W. (2008). Functional forms for the negative binomial model for count 
#' data. Economics Letters, 99(3), 585-590.
#' 
#' Pararai, M., Famoye, F., & Lee, C. (2006). Generalized Poisson regression 
#' model for underreported counts. Advances and applications in Statistics, 
#' 6(3), 305-322.
#' 
#' Pararai, M., Famoye, F., & Lee, C. (2010). Generalized poisson-poisson 
#' mixture model for misreported counts with an application to smoking data. 
#' Journal of Data Science, 8(4), 607-617.
#' 
#' Rigby, R. A., Stasinopoulos, D. M., & Akantziliotou, C. (2008). A framework 
#' for modelling overdispersed count data, including the Poisson-shifted 
#' generalized inverse Gaussian distribution. Computational Statistics & Data 
#' Analysis, 53(2), 381-393.
#' 
#' Wood, J.S., Eric T. Donnell, & Christopher J. Fariss. "A method to account 
#' for and estimate underreporting in crash frequency research." Accident 
#' Analysis & Prevention 95 (2016): 57-66.
#' 
#' Zou, Y., Lord, D., & Zhang, Y. (2012). Analyzing highly dispersed crash data 
#' using the Sichel generalized additive models for location, scale and shape. 
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
  
  # Get the parameters and probability function
  family <- toupper(str_replace_all(family, "[^[:alnum:]]", "")) # remove non-alphanumeric characters from the family name and ensure all upper case
  params <- get_params(family)
  
  method <- toupper(str_replace_all(method, "[^abcfghmnrsABCFGHMNRS]", "")) # clean the method name
  if (!(method %in% c("SANN", "NM", "BFGS", "BFGSR", "CG", "NR", "BHHH"))){
    print('Method must be one of: "SANN", "NM", "BFGS", "BFGSR", "CG", "NR", or "BHHH". Switching to "NM.')
    method = "NM"
  }
  
  
  
  # Prepare model matrices
  mod1_frame <- model.frame(formula, data)
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
    X_offset <- as.matrix(data[, offset, drop = FALSE])
  }
  
  if(any(grepl("offset", deparse(formula)))) { # If offset() is used in the formula, use that as the offset
    offset_variable <- str_extract(deparse(formula), "(?<=offset\\().*?(?=\\))")
    X_offset <- as.matrix(data[, offset, drop = FALSE])
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
  
  if (is.null(params[[1]])) {
    N_alpha = 0
  } else if (is.null(dis_param_formula_1)) {
    N_alpha = 1
  } else {
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
    
    # Use glm (Poisson) as starting values for betas
    p_model <- glm(formula, data = data, family = poisson(link = "log"))
    start <- p_model$coefficients
    
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
    weights.df <- data %>% pull(weights)
  }
  
  # Define the main function for computing log-likelihood
  logLikFunc <- function(p) {
    # Ensure local_probFunc is available in this scope
    local_probFunc <- get_probFunc(family)
    
    # Critical Check: If get_probFunc returned NULL (e.g. invalid family matching), we cannot proceed.
    # This prevents the obscure "could not find function 'local_probFunc'" or "attempt to apply non-function" error later.
    if(is.null(local_probFunc)) {
      stop(paste0("Probability function not found for family: ", family, ". Please check family name and helpers.R definition."))
    }
    
    coefs <- as.array(p)
    fixed_coefs <- as.vector(head(coefs, N_predictors))
    
    # Initialize alpha and sigma
    alpha <- NULL
    sigma <- NULL
    
    if (N_alpha > 0) {
      if (!is.null(dis_param_formula_1)){
        alpha_coefs <- as.vector(coefs[(N_predictors + 1):(N_predictors + N_alpha)])
        alpha <- exp(mod_alpha_frame %*% alpha_coefs)
      } else {
        alpha <- exp(coefs[(N_predictors + 1)])
      }
    }
    
    if (N_sigma > 0) {
      if (!is.null(dis_param_formula_2)){
        sigma_coefs <- as.vector(coefs[(N_predictors + N_alpha + 1):(N_predictors + N_alpha + N_sigma)])
        sigma <- exp(mod_sigma_frame %*% sigma_coefs)
      } else {
        sigma <- exp(coefs[(N_predictors + N_alpha + 1)])
      }
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
      if (length(offset)>1){
        X_offset_i = rowSums(X_offset)
        predicted <- exp(X %*% fixed_coefs + X_offset_i)*underreport_prob
      }else{
        predicted <- exp(X %*% fixed_coefs + X_offset)*underreport_prob
      }
    } else {
      predicted <- exp(X %*% fixed_coefs)*underreport_prob
    }
    
    # Call the probability function using the local variable
    probs <- local_probFunc(y=y, predicted=predicted, alpha=alpha, sigma=sigma, 
                            haltons=haltons, normed_haltons=normed_haltons)
    
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
    # Reuse current estimates as starting values for bootstrap
    start_vals_boot <- fit$estimate
    
    models <- map(bs.data$strap, ~ mod.boot(data = ., 
                                            formula = formula, 
                                            family = family, 
                                            offset = offset, 
                                            weights = weights, 
                                            dis_param_formula_1 = dis_param_formula_1, 
                                            dis_param_formula_2 = dis_param_formula_2, 
                                            underreport_formula = underreport_formula, 
                                            underreport_family = underreport_family, 
                                            ndraws = ndraws, 
                                            method = method, 
                                            max.iters = max.iters, 
                                            start.vals = start_vals_boot))
    tidied <- map_df(models, broom::tidy, .id = "id")
    
    SE <- tidied %>%
      group_by(term) %>%
      reframe(sd = sd(estimate, na.rm=TRUE))
    
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