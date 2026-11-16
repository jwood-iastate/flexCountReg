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
#'   data frame that should be used as an offset (i.e., included but forced to
#'   have a coefficient of 1). When several names are given, the offsets are
#'   summed. Alternatively, an offset can be given in the formula with
#'   `offset()`. Use only one of the two approaches: supplying both is an
#'   error.
#' @param weights the name of a variable in the data frame, or a numeric vector
#'   with one value per row of `data`, to be used as frequency weights.
#' @param verbose an optional parameter. If `TRUE`, the function will print out
#'   the progress of the model fitting. Default is `FALSE`.
#' @param dis_param_formula_1 a symbolic description of the model for the
#'   natural log of the dispersion parameter or first parameter of the count
#'   distribution used. Further details are provided below.
#' @param dis_param_formula_2 a symbolic description of the model for the second
#'   parameter of the count distribution used. Further details are provided
#'   below.
#' @param underreport_formula an optional formula to estimate the underreporting
#'   for any of the count model options. The underreporting is estimated as a
#'   function (logit or probit) of the predictors in the model. For the model to
#'   be tractable, the independent variables cannot be the exact same as the
#'   count model. The default is `NULL`.
#' @param underreport_family the name of the distribution/model type to estimate
#'   the underreporting portion of the model when `underreport_formula` is
#'   specified. The default is "logit" for a binary logistic regression model.
#'   The other option is "probit" for a probit model.
#' @param ndraws The number of Halton draws for integrating the distribution
#'   being compounded with the Poisson distribution when there is not a
#'   closed-form solution. Default is 1500. It is recommended to test different
#'   numbers of draws to determine if the model is stable (i.e., doesn't change
#'   or has minimal change as the number of draws changes within a reasonable
#'   range).
#' @param method Optimization method to be used for maximum likelihood
#'   estimation. One of "NM" (Nelder-Mead, the default), "BFGS", "BFGSR", "CG",
#'   "NR", "BHHH", or "SANN". See `maxLik` documentation for details.
#' @param max.iters Maximum number of iterations for the optimization method.
#' @param start.vals Optional vector of starting values for the optimization.
#'   It must have one value per parameter. If the vector is named, the values
#'   are matched to the parameters by name.
#' @param stderr Type of standard errors to use. The default is "normal"
#'   (based on the Hessian). Other options include "boot" for bootstrapped
#'   standard errors, or "robust" for robust (sandwich) standard errors.
#' @param bootstraps Optional integer specifying the number of bootstrap samples
#'   to be used for estimating standard errors when `stderr`= "boot". Defaults
#'   to 200 when `stderr = "boot"` and `bootstraps` is not given.
#' @param keep_data logical; if \code{FALSE}, the training data are not stored
#'   in the fitted object (saving memory). \code{predict()} then needs
#'   \code{newdata}, and \code{fitted()}/\code{residuals()} are unavailable.
#' @param zero_formula an optional one-sided formula for a zero model
#'   (\code{~ 1} for a constant), turning the model into a zero-inflated or
#'   hurdle model (see \code{zero_type}) for any \code{family}. The zero model
#'   is a logit model. Its coefficients are named \code{Zero:<variable>}.
#'   Cannot be combined with \code{underreport_formula}.
#' @param zero_type \code{"inflated"} (default) or \code{"hurdle"}. In a
#'   zero-inflated model, a structural zero occurs with probability
#'   \eqn{\pi = logit^{-1}(z\gamma)} and otherwise the count follows the
#'   family: \eqn{P(0) = \pi + (1 - \pi) f(0)} and
#'   \eqn{P(y) = (1 - \pi) f(y)} for \eqn{y > 0}. In a hurdle model,
#'   \eqn{q = logit^{-1}(z\gamma)} is the probability of a positive count and
#'   positive counts follow the zero-truncated family:
#'   \eqn{P(0) = 1 - q} and \eqn{P(y) = q f(y) / (1 - f(0))}.
#' @param engine the engine to use for the integration if the Poisson Lognormal 
#'   (PLN) distribution is used. Options include "halton" which uses Halton 
#'   draws or "poilog" which uses the poilog package. The default is "halton".
#' 
#' 
#' @description
#' The purpose of this function is to estimate count regression models using
#' maximum likelihood estimation (MLE) or Maximum Simulated Likelihood
#' Estimation (MSLE). The function can estimate the following models:
#' \itemize{
#'  \item Poisson (Poisson)
#'  \item Generalized Poisson Version 1 (GP1)
#'  \item Generalized Poisson Version 2 (GP2)
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
#' \item "GP1" for Generalized Poisson Version 1 Distribution with a log link
#' \item "GP2" for Generalized Poisson Version 2 Distribution with a log link
#'  \item "NB1" for Negative Binomial 1 distribution with a log link.
#'  \item "NB2" for Negative Binomial 2 distribution with a log link (i.e., the 
#'        standard negative binomial model).
#'  \item "NBP" for Negative Binomial P distribution with a log link.
#'  \item "PLN" for Poisson-Lognormal distribution with a log link.
#'  \item "PGE" for Poisson-Generalized-Exponential distribution with a log 
#'        link.
#'  \item "PIG1" for Poisson-Inverse-Gaussian Type-1 distribution with a log
#'  link.
#'  \item "PIG2" for Poisson-Inverse-Gaussian Type-2 distribution with a log
#'  link.
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
#'  functions rather than constants in the model. For example, if the user wants
#'  to estimate the overdispersion parameter of the Negative Binomial 2
#'  distribution as a function of the variable `x1` and `x2`, the user would
#'  specify `dis_param_formula_1 = ~ x1 + x2`. In the case of the Negative
#'  Binomial distributions, the model is known as a Generalized Negative
#'  Binomial model when the overdispersion parameter is specified as a function.
#'  
#'  The function linking the distribution parameters to the predictors is:
#'  \deqn{Param = \exp((Intercept) + \sum \beta X)}
#'  
#'  The parameters for the different models are as follows:
#'  
#'  For `dis_param_formula_1`, the models are for the parameters:
#'  \itemize{
#'  \item \eqn{ln(\phi+1)} for the Generalized Poisson Version 1.
#'  \item \eqn{ln(\alpha+1)} for the Generalized Poisson Version 2.
#'  \item \eqn{ln(\alpha)} for the Negative Binomial 1 model.
#'  \item \eqn{ln(\alpha)} for the Negative Binomial 2 model.
#'  \item \eqn{ln(\alpha)} for the Negative Binomial P model.
#'  \item \eqn{ln(\sigma)} for the Poisson-Lognormal model.
#'  \item \eqn{ln(shape)} for the Poisson-Generalized-Exponential model.
#'  \item \eqn{ln(\eta)} for the Poisson-Inverse-Gaussian model.
#'  \item \eqn{ln(\eta)} for the Poisson-Inverse-Gamma model.
#'  \item \eqn{ln(\theta)} for the Poisson-Lindley model.
#'  \item \eqn{ln(\theta)} for the Poisson-Lindley-Gamma model.
#'  \item \eqn{ln(\theta)} for the Poisson-Lindley-Lognormal model.
#'  \item \eqn{ln(\alpha)} for the Poisson-Weibull model.
#'  \item \eqn{\gamma} for the Sichel model. This parameter can be negative,
#'  so it is modeled directly (identity link), \eqn{\gamma = W\delta}, rather
#'  than through \eqn{\exp()}.
#'  \item \eqn{ln(k)} for the Generalized Waring model.
#'  \item \eqn{ln(\nu)} for the Conway-Maxwell-Poisson model.
#'  }
#'  
#'  For `dis_param_formula_2`, the models are for the parameters:
#'  \itemize{
#'  \item Not Applicable for the Generalized Poisson Version 1 model.
#'  \item Not Applicable for the Generalized Poisson Version 2 model.
#'  \item Not Applicable for the Negative Binomial 1 model.
#'  \item Not Applicable for the Negative Binomial 2 model.
#'  \item \eqn{ln(p)} for the Negative Binomial P model.
#'  \item Not Applicable for the Poisson-Lognormal model.
#'  \item \eqn{ln(scale)} for the Poisson-Generalized-Exponential model.
#'  \item Not Applicable for the Poisson-Inverse-Gaussian model.
#'  \item Not Applicable for the Poisson-Inverse-Gamma model.
#'  \item Not Applicable for the Poisson-Lindley model.
#'  \item \eqn{ln(\alpha)} for the Poisson-Lindley-Gamma model.
#'  \item \eqn{ln(\sigma)} for the Poisson-Lindley-Lognormal model.
#'  \item \eqn{ln(\sigma)} for the Poisson-Weibull model.
#'  \item \eqn{ln(\sigma)} for the Sichel model.
#'  \item \eqn{ln(\rho-1)} for the Generalized Waring model.
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
#' This implements the Poisson regression model using Maximum Likelihood 
#' Estimation (MLE), as opposed to the Iteratively Reweighted Least Squares (IRLS) 
#' method used in the `glm` function.
#' 
#' The PMF and log-likelihood functions are:
#' \deqn{P(Y = y) = \frac{e^{-\mu} \mu^y}{y!}}
#' \deqn{LL_{\text{Poisson}}(\beta) = 
#'        \sum_{i=1}^n \left[ -\mu_i + y_i \ln(\mu_i) - \ln(y_i!) \right]}
#' 
#' The mean is:
#' \deqn{\mu = exp(X\beta)}
#' 
#' The variance is:
#' \deqn{\text{Var}(Y) = \mu}
#' 
#' ## Generalized Poisson Version 1 Model
#' This implements the Generalized Poisson Version 1 model using MLE.
#' 
#' The PMF is:
#' \deqn{
#' f(y|\phi,\mu)=\frac{\mu(\mu+\phi y)^{y-1} exp\left(-\frac{\mu+\phi y}
#' {1+\phi}\right)}{(1+\phi)^y y!}
#' }
#'
#' Where \eqn{\phi} is a dispersion parameter with the restriction that
#' \eqn{\phi>-1}, \eqn{\mu>0} is the mean value, and \eqn{y} is a non-negative
#' integer. The mean is modeled as:
#' \deqn{\mu = exp(X\beta)} 
#'
#' The variance of the GP-1 distribution is:
#' \deqn{\sigma^2=(1+\phi)^2 \mu}
#' 
#' If \eqn{\phi>0}, the distribution is overdispersed. If \eqn{\phi=0}, the 
#' distribution is equidispersed. If \eqn{\phi<0}, the distribution is 
#' underdispersed.
#' 
#' Furthermore, \eqn{\phi>-1} is required for this distribution. When 
#' \eqn{\phi<0}, there is also a maximum value of support for the integer value 
#' \eqn{y}. This is \eqn{y_max=\left\lfloor -\frac{\mu}{\phi} \right\rfloor}.
#' To account for the range of possible values, this is modeled as 
#' \deqn{ln(\phi+1)}.
#' 
#' ## Generalized Poisson Version 2 Model
#' This implements the Generalized Poisson Version 2 model using MLE.
#' 
#' The PMF is:
#' \deqn{
#' f(y|\mu,\alpha)=
#' \frac{\mu(\mu+\alpha\mu y)^{y-1}
#' \exp\left(-\frac{\mu+\alpha\mu y}{1+\alpha\mu}\right)}
#' {(1+\alpha\mu)^y y!}
#' }
#'
#' Where \eqn{\alpha} is a dispersion parameter with the restriction that
#' \eqn{1+\alpha\mu>0}, \eqn{\mu>0} is the mean value, and \eqn{y} is a
#' non-negative integer. The mean is modeled as:
#' \deqn{\mu = exp(X\beta)} 
#'
#' The variance of the GP-2 distribution is:
#' \deqn{\sigma^2=\mu(1+\alpha\mu)^2}
#'
#' If \eqn{\alpha > 0}, the distribution is overdispersed. If
#' \eqn{\alpha = 0}, the distribution is equidispersed and reduces to the
#' ordinary Poisson distribution. If \eqn{\alpha < 0} and
#' \eqn{1 + \alpha\mu > 0}, the distribution is underdispersed and has
#' finite support.
#'
#' Under the underdispersed case, the largest possible value of \eqn{y}
#' is bounded above by \eqn{\left\lceil -1/\alpha \right\rceil - 1},
#' i.e., the largest integer satisfying \eqn{1 + \alpha y > 0}.
#' To account for the range of possible values, this is modeled as 
#' \deqn{ln(\alpha+1)}.
#' 
#' ## Negative Binomial Models
#' 
#' The NB-1, NB-2, and NB-P versions of the negative binomial distribution are 
#' based on Greene (2008).  The details of each of these are provided below.
#' 
#' ### NB-1 Model
#' The PMF and log-likelihood functions are:
#' \deqn{P(Y = y) = 
#'     \frac{\Gamma(y + \frac{\mu}{\alpha})}{y! \, \Gamma(\frac{\mu}{\alpha})}
#'     \left( \frac{\frac{\mu}{\alpha}}{\frac{\mu}{\alpha} + \mu} \right)
#'     ^{\frac{\mu}{\alpha}} \left( \frac{\mu}{\frac{\mu}{\alpha} + \mu} 
#'     \right)^y}
#' \deqn{LL_{\text{NB1}}(\beta, \alpha) = \sum_{i=1}^n \left[ \ln \Gamma\left( 
#' y_i + \frac{\mu_i}{\alpha} \right) - \ln \Gamma\left( \frac{\mu_i}{\alpha} 
#' \right) - \ln y_i! + \frac{\mu_i}{\alpha} \ln \left( 
#' \frac{\frac{\mu_i}{\alpha}}{\frac{\mu_i}{\alpha} + \mu_i} \right) + y_i 
#' \ln \left( \frac{\mu_i}{\frac{\mu_i}{\alpha} + \mu_i} \right) \right]}
#' 
#' The mean is:
#' \deqn{\mu = exp(X\beta)}
#' 
#' The variance is:
#' \deqn{\text{Var}(Y) = \mu + \alpha\mu}
#' 
#' ### NB-2 Model
#' The PMF and log-likelihood functions are:
#' \deqn{P(Y = y) = \frac{\Gamma(y + \alpha)}{y! \, \Gamma(\alpha)} 
#' \left( \frac{\alpha}{\alpha + \mu} \right)^\alpha 
#' \left( \frac{\mu}{\alpha + \mu} \right)^y}
#' \deqn{LL_{\text{NB2}} = \sum_{i=1}^n \left[ \ln \Gamma(y_i + \alpha) - 
#' \ln \Gamma(\alpha) - \ln y_i! + \alpha \ln \left( \frac{\alpha}{\alpha + 
#' \mu_i} \right) + y_i \ln \left( \frac{\mu_i}{\alpha + \mu_i} \right) \right]}
#' 
#' The mean is:
#' \deqn{\mu = exp(X\beta)}
#' 
#' The variance is:
#' \deqn{\text{Var}(Y) = \mu + \alpha\mu^2}
#' 
#' ### NB-P Model
#' The PMF and log-likelihood functions are:
#' \deqn{P(Y = y) = \frac{\Gamma(y + \frac{\mu^{2-p}}{\alpha})}{y! \, 
#' \Gamma(\frac{\mu^{2-p}}{\alpha})} \left( \frac{\frac{\mu^{2-p}}{\alpha}}
#' {\frac{\mu^{2-p}}{\alpha} + \mu} \right)^{\frac{\mu^{2-p}}{\alpha}} 
#' \left( \frac{\mu}{\frac{\mu^{2-p}}{\alpha} + \mu} \right)^y}
#' \deqn{LL_{\text{NBP}}(\beta, \alpha, p) = \sum_{i=1}^n \left[ \ln 
#' \Gamma\left( y_i + \frac{\mu_i^{2-p}}{\alpha} \right) - \ln \Gamma\left( 
#' \frac{\mu_i^{2-p}}{\alpha} \right) - \ln y_i! + \frac{\mu_i^{2-p}}{\alpha} 
#' \ln \left( \frac{\frac{\mu_i^{2-p}}{\alpha}}{\frac{\mu_i^{2-p}}{\alpha} + 
#' \mu_i} \right) + y_i \ln \left( \frac{\mu_i}{\frac{\mu_i^{2-p}}{\alpha} + 
#' \mu_i} \right) \right]}
#' 
#' The mean is:
#' \deqn{\mu = exp(X\beta)}
#' 
#' The variance is:
#' \deqn{\text{Var}(Y) = \mu + \alpha\mu^p}
#' 
#' ## Poisson-Lognormal (PLN) Model
#' The compound Probability Mass Function(PMF) for the Poisson-Lognormal 
#' distribution is:
#' \deqn{f(y|\lambda,\sigma)=\int_0^\infty \frac{\lambda^y x^y e^{-\lambda x}}
#' {y!}\frac{exp\left(-\frac{ln^2(x)}{2\sigma^2} \right)}{x\sigma\sqrt{2\pi}}dx}
#'
#' Where \eqn{\sigma} is a parameter for the lognormal distribution with the 
#' restriction \eqn{\sigma>0}, and \eqn{y} is a non-negative integer.
#'
#' The expected value of the distribution is:
#' \deqn{E[y]=e^{X\beta+\sigma^2/2} = \mu e^{\sigma^2/2}}
#' 
#' The t-statistics and p-values for the coefficients related to ln(sigma) are,
#' by default, testing if the coefficients are different from a value of 0. This
#' has little practical meaning given that they are coefficients for ln(sigma).
#' They are not testing if the coefficients have statistical significance in
#' terms of improvement over a Poisson model. To test whether the
#' Poisson-Lognormal model fits significantly better than the Poisson model,
#' use a likelihood-ratio test with \code{\link{regCompTest}}.
#' 
#' ## Poisson Generalized-Exponential (PGE) Model
#' The Generalized Exponential distribution can be written as a function with 
#' a shape parameter \eqn{\alpha>0} and scale parameter \eqn{\gamma>0}. The
#' distribution has strictly positive continuous values. The PDF of the
#' distribution is:
#' \deqn{f(x|\alpha,\gamma)=\frac{\alpha}{\gamma}\left(1-e^{-\frac{x}{\gamma}}
#' \right)^{\alpha-1}e^{-\frac{x}{\gamma}}} 
#' 
#' Thus, the compound Probability Mass Function(PMF) for the PGE distribution
#' is:
#' \deqn{f(y|\lambda,\alpha,\beta)=\int_0^\infty \frac{\lambda^y x^y 
#' e^{-\lambda x}}{y!}\frac{\alpha}{\gamma}\left(1-e^{-\frac{x}{\gamma}}
#' \right)^{\alpha-1}e^{-\frac{x}{\gamma}} dx}
#' 
#' The expected value of the distribution is:
#' \deqn{E[y]=\mu=\lambda \left(\frac{\psi(\alpha+1)-\psi(1)}{\gamma}\right)}
#' 
#' Where \eqn{\psi(\cdot)} is the digamma function.
#' 
#' The variance is:
#' \deqn{\sigma^2=\lambda \left(\frac{\psi(\alpha+1)-\psi(1)}{\gamma}\right) + 
#' \left(\frac{-\psi'(\alpha+1)+\psi'(1)}{\gamma^2}\right)\lambda^2}
#' 
#' Where \eqn{\psi'(\cdot)} is the trigamma function.
#' 
#' To ensure that \eqn{\mu=e^{X\beta}}, \eqn{\lambda} is replaced with:
#' \deqn{\lambda=\frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)}}
#' 
#' This results in:
#' \deqn{f(y|\mu,\alpha,\beta)=\int_0^\infty \frac{\left(\frac{\gamma 
#' e^{X\beta}}{\psi(\alpha+1)-\psi(1)}\right)^y x^y e^{-\left(\frac{\gamma 
#' e^{X\beta}}{\psi(\alpha+1)-\psi(1)}\right) x}}{y!}\frac{\alpha}{\gamma}
#' \left(1-e^{-\frac{x}{\gamma}}\right)^{\alpha-1}e^{-\frac{x}{\gamma}} dx}
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
#' \deqn{f(x|\eta,\mu)=\frac{2\left(\mu\left(\frac{1}{\eta}+1\right)\right)
#' ^{\frac{x+\frac{1}{\eta}+2}{2}}}{x!\Gamma\left(\frac{1}{\eta}+2\right)}
#' K_{x-\frac{1}{\eta}-2}\left(2\sqrt{\mu\left(\frac{1}{\eta}+1\right)}\right)}
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
#' independent variables. Thus, the LR test result indicates the statistical 
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
#' \deqn{\sigma^2=\mu+\left(2\alpha+1-\frac{2(1+\alpha)}
#' {(\theta+2)^2}\right)\mu^2}
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
#' \deqn{\sigma^2=\mu+\left(\frac{1-\frac{2}{(\theta+2)^2}}{e^{\frac{\sigma^2}
#' {2}}}+e^{\sigma^2}-1\right)\mu^2}
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
#' \deqn{P(y|\lambda,\alpha,\sigma) = \int_0^\infty \frac{e^{-\lambda x} 
#' \lambda^y x^y }{y!} \left(\frac{\alpha}{\sigma}\right)\left(\frac{x}{\sigma}
#' \right)^{\alpha-1}e^{-\left(\frac{x}{\sigma}\right)^\alpha} dx}
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
#' \deqn{V[Y]=\mu+\left(\frac{\Gamma\left(1+\frac{2}{\alpha}\right)}
#' {\Gamma\left(1+\frac{1}{\alpha}\right)^2}-1\right)\mu^2}
#' 
#' ## Sichel (SI) Model
#' The compound Probability Mass Function (PMF) for the Sichel distribution uses
#' the formulation from Zou, Lord and Zhang (2012) and Rigby et al. (2008):
#' \deqn{f(y|\mu, \sigma, \gamma)=\frac{\left(\frac{\mu}{c}\right)^y 
#' K_{y+\gamma}(\alpha)}{K_\gamma(1/\sigma)y!(\alpha\sigma)^{y+\gamma}}}
#' 
#' Where \eqn{\sigma} and \eqn{\gamma} are distribution parameters with 
#' \eqn{\sigma>0}, 
#' \deqn{c=\frac{K_{\gamma+1}(1/\sigma)}{K_\gamma(1/\sigma)}}, 
#' \deqn{\alpha^2=\sigma^{-2}+2\mu(c\sigma)^{-1}}, 
#' a mean value of \eqn{\mu}, \eqn{y} is a non-negative integer, and 
#' \deqn{K_j(x)} is a modified Bessel function of the third kind with order 
#' \deqn{j} and argument \deqn{x}.
#'
#' The variance of the distribution is:
#' \deqn{\sigma^2=\mu+\left(\frac{2\sigma(\gamma+1)}{c}+\frac{1}{c^2}-1\right)
#' \mu^2}
#' 
#' ## Generalized Waring (GW) Model
#' The following are the PMF, mean, and variance used for the Generalized
#' Waring model. The regression formulation replaces parameter \eqn{\alpha}
#' with the mean \eqn{\mu}, while retaining parameters \eqn{k} and \eqn{\rho}.
#'
#' In the original parameterization, the PMF is:
#' \deqn{
#' P(Y=y)=
#' \frac{
#' \Gamma(\alpha+y)\Gamma(k+y)\Gamma(\rho+k)\Gamma(\alpha+\rho)
#' }{
#' y!\Gamma(\alpha)\Gamma(k)\Gamma(\rho)
#' \Gamma(\alpha+k+\rho+y)
#' },
#' }
#' where \eqn{y=0,1,2,\ldots} and \eqn{\alpha,k,\rho>0}.
#'
#' The mean and variance are:
#' \deqn{
#' E(Y)=\frac{\alpha k}{\rho-1}, \qquad \rho>1,
#' }
#' \deqn{
#' \operatorname{Var}(Y)=
#' \frac{\alpha k(\alpha+\rho-1)(k+\rho-1)}
#' {(\rho-1)^2(\rho-2)}, \qquad \rho>2.
#' }
#'
#' In accident applications, the variance can be decomposed into components
#' attributed to randomness, proneness, and liability. These are variance
#' components, rather than individual interpretations of \eqn{\alpha},
#' \eqn{k}, and \eqn{\rho}.
#'
#' To parameterize the model in terms of its mean, set:
#' \deqn{
#' \mu=e^{X\beta}, \qquad
#' \alpha=\frac{\mu(\rho-1)}{k}.
#' }
#'
#' Substituting into the PMF gives:
#' \deqn{
#' P(Y=y \mid X)=
#' \frac{
#' \Gamma\left(\frac{\mu(\rho-1)}{k}+y\right)
#' \Gamma(k+y)\Gamma(\rho+k)
#' \Gamma\left(\frac{\mu(\rho-1)}{k}+\rho\right)
#' }{
#' y!\Gamma\left(\frac{\mu(\rho-1)}{k}\right)
#' \Gamma(k)\Gamma(\rho)
#' \Gamma\left(\frac{\mu(\rho-1)}{k}+k+\rho+y\right)
#' }.
#' }
#'
#' This results in a regression model with:
#' \deqn{
#' E(Y \mid X)=\mu=e^{X\beta},
#' }
#' \deqn{
#' \operatorname{Var}(Y \mid X)=
#' \frac{\mu(k+\mu)(k+\rho-1)}{k(\rho-2)}
#' =
#' \frac{k+\rho-1}{\rho-2}\mu+
#' \frac{k+\rho-1}{k(\rho-2)}\mu^2,
#' \qquad \rho>2.
#' }
#'
#' The original distribution is valid for all \eqn{\alpha,k,\rho>0},
#' including \eqn{\rho=1} and \eqn{\rho=2}. Its mean is infinite when
#' \eqn{0<\rho\leq1}. When \eqn{1<\rho\leq2}, the mean is finite but
#' the variance is infinite. A finite variance requires \eqn{\rho>2}.
#' The mean-based regression parameterization requires \eqn{\mu>0},
#' \eqn{k>0}, and \eqn{\rho>1}.
#' 
#' ## Conway-Maxwell-Poisson (COM) Model
#' The following is the the PMF for the COM model.
#' \deqn{f(x|\lambda, \nu)=\frac{\lambda^x}{(x!)^\nu Z(\lambda,\nu)}}
#' 
#' Where \eqn{\lambda} and \eqn{\nu} are distribution parameters with 
#' \eqn{\lambda>0} and \eqn{\nu>0}, and \eqn{Z(\lambda,\nu)} is the normalizing 
#' constant.
#' 
#' The normalizing constant is given by:
#' \deqn{Z(\lambda,\nu)=\sum_{n=0}^{\infty}\frac{\lambda^n}{(n!)^\nu}}
#' 
#' The mean and variance are:
#' \deqn{\mu=e^{X\beta}=\lambda \frac{\delta}{\delta \lambda} \log(Z(\lambda,
#' \nu))}
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
#' \deqn{\mu_{observed}=\mu_{true}\cdot P(\text{event is reported})}
#' 
#' This allows the inference of both the true event count and the probability of
#' the event being reported as a function of independent variables.
#'  
#' @returns
#' An object of class `flexCountReg` which is a list with the following components:
#' \itemize{
#'  \item model: the fitted model object.
#'  \item data: the data frame used to fit the model.
#'  \item call: the matched call.
#'  \item formula: the formula used to fit the model.
#' }
#' 
#' @importFrom randtoolbox halton
#' @importFrom stats model.frame model.matrix model.response
#' @importFrom tibble deframe
#' @importFrom maxLik maxLik
#' @importFrom sandwich sandwich
#' @include pinvgaus.R pinvgamma.R ppoislogn.R plindLnorm.R plindGamma.R 
#' @include psichel.R Generalized-Waring.R ppoisGE.R psichel.R plind.R helpers.R
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
#' # Estimate a Poisson-Lognormal model (integration with the poilog package)
#' if (requireNamespace("poilog", quietly = TRUE)) {
#'   pln <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
#'                   data = washington_roads, family = "PLN", engine="poilog")
#'   summary(pln)
#' }
#' 
#' # Estimate a Negative Binomial with underreporting (probit)
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
#' 
#' # Estimate a Generalized-Poisson Version 1 model
#' gp1_model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + 
#'                       AADT10kplus,
#'                       data = washington_roads, family = "GP1", method="BHHH")
#' summary(gp1_model)
#' 
#' # Estimate a Generalized-Poisson Version 2 model
#' gp2_model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + 
#'                       AADT10kplus,
#'                       data = washington_roads, family = "GP2", method="BHHH")
#' summary(gp2_model)
#' }
#' 
#' @references
#' Greene, W. (2008). Functional forms for the negative binomial model for count
#' data. Economics Letters, 99(3), 585-590.
#' 
#' Pararai, M., Famoye, F., & Lee, C. (2006). Generalized Poisson regression 
#' model for underreported counts. Advances and applications in Statistics, 
#' 6(3), 305-322.
#' 
#' Pararai, M., Famoye, F., & Lee, C. (2010). Generalized Poisson-Poisson 
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
#' Accident Analysis & Prevention, 45.
#' 
#' 
#' @importFrom Rcpp sourceCpp
#' 
#' @useDynLib flexCountReg
#' @export
countreg <- function(formula, data, family = "NB2", offset = NULL, 
                     weights = NULL, verbose = FALSE, 
                     dis_param_formula_1 = NULL, 
                     dis_param_formula_2 = NULL, underreport_formula = NULL,
                     underreport_family = "logit",
                     ndraws = 1500, method = "NM", 
                     max.iters = 1000, start.vals = NULL, 
                     stderr = "normal", bootstraps = NULL,
                     engine = "halton", keep_data = TRUE,
                     zero_formula = NULL,
                     zero_type = c("inflated", "hurdle")) {
  
  print.level <- if (verbose) 2 else 0
  
  # ---- Validate options ------------------------------------------------------
  # Remove non-alphanumeric characters from the family name, upper case
  family <- toupper(gsub("[^[:alnum:]]", "", family))
  if (!family %in% .fcr_families) {
    stop("Unknown `family` \"", family, "\". Must be one of: ",
         toString(.fcr_families), ".", call. = FALSE)
  }
  params <- get_params(family)
  
  engine <- match.arg(engine, c("halton", "poilog"))
  if (engine == "poilog" && family == "PLN") {
    rlang::check_installed("poilog", reason = "for `engine = \"poilog\"`.")
  }
  method <- .fcr_check_method(method)
  underreport_family <- match.arg(tolower(underreport_family),
                                  c("logit", "probit"))
  stderr <- match.arg(tolower(stderr), c("normal", "robust", "boot"))
  zero_type <- match.arg(zero_type)
  if (!is.null(zero_formula) && !is.null(underreport_formula)) {
    stop("`zero_formula` and `underreport_formula` cannot be combined: both ",
         "act on the zeros, so they are not separately identified.",
         call. = FALSE)
  }
  if (stderr == "boot") {
    if (is.null(bootstraps)) {
      bootstraps <- 200L
      message("`bootstraps` not given; using 200 bootstrap samples.")
    }
    if (!is.numeric(bootstraps) || length(bootstraps) != 1L ||
        bootstraps < 2) {
      stop("`bootstraps` must be a single number of at least 2.",
           call. = FALSE)
    }
    bootstraps <- as.integer(bootstraps)
  }
  
  if (is.null(params[[1]]) && !is.null(dis_param_formula_1)) {
    stop("`dis_param_formula_1` was given, but family \"", family,
         "\" has no distribution parameter to model.", call. = FALSE)
  }
  if (is.null(params[[2]]) && !is.null(dis_param_formula_2)) {
    stop("`dis_param_formula_2` was given, but family \"", family,
         "\" has no second distribution parameter to model.", call. = FALSE)
  }
  
  # ---- Data: one complete-case model frame for everything --------------------
  prep <- prepare_model_data(formula, data,
                             aux_formulas = list(
                               dis1 = dis_param_formula_1,
                               dis2 = dis_param_formula_2,
                               underreport = underreport_formula,
                               zero = zero_formula),
                             offset = offset, weights = weights)
  # Keep the formulas with `.` expanded, for prediction, methods and refits
  formula <- prep$formula
  dis_param_formula_1 <- prep$aux_formulas$dis1
  dis_param_formula_2 <- prep$aux_formulas$dis2
  underreport_formula <- prep$aux_formulas$underreport
  zero_formula <- prep$aux_formulas$zero
  if (prep$n_dropped > 0) {
    message(prep$n_dropped, " row(s) with missing values were removed.")
  }
  
  X <- prep$X
  y <- prep$y
  X_offset <- prep$offset
  weights.df <- prep$weights
  mod_alpha_frame <- prep$aux$dis1
  mod_sigma_frame <- prep$aux$dis2
  X_underreport <- prep$aux$underreport
  X_zero <- prep$aux$zero
  
  # ---- Parameter bookkeeping -------------------------------------------------
  N_predictors <- ncol(X)
  x_names <- colnames(X)
  
  if (is.null(params[[1]])) {
    N_alpha <- 0
  } else if (is.null(mod_alpha_frame)) {
    N_alpha <- 1
    x_names <- c(x_names, params[[1]])
  } else {
    N_alpha <- ncol(mod_alpha_frame)
    x_names <- c(x_names, paste0(params[[1]], ":", colnames(mod_alpha_frame)))
  }
  
  if (is.null(params[[2]])) {
    N_sigma <- 0
  } else if (is.null(mod_sigma_frame)) {
    N_sigma <- 1
    x_names <- c(x_names, params[[2]])
  } else {
    N_sigma <- ncol(mod_sigma_frame)
    x_names <- c(x_names, paste0(params[[2]], ":", colnames(mod_sigma_frame)))
  }
  
  if (!is.null(X_underreport)) {
    N_underreport <- ncol(X_underreport)
    x_names <- c(x_names, paste0("Underreporting:", colnames(X_underreport)))
  } else {
    N_underreport <- 0
  }
  if (!is.null(X_zero)) {
    N_zero <- ncol(X_zero)
    x_names <- c(x_names, paste0("Zero:", colnames(X_zero)))
  } else {
    N_zero <- 0
  }
  N_par <- length(x_names)
  idx_under <- N_predictors + N_alpha + N_sigma + seq_len(N_underreport)
  idx_zero <- N_predictors + N_alpha + N_sigma + N_underreport +
    seq_len(N_zero)
  
  # Generate Halton draws to use as quantile values
  haltons <- randtoolbox::halton(ndraws)
  normed_haltons <- qnorm(haltons)
  
  # ---- Starting values -------------------------------------------------------
  if (is.null(start.vals)) {
    # Poisson GLM with the same rows, offset and weights for the betas
    p_model <- stats::glm.fit(X, y, weights = weights.df, offset = X_offset,
                              family = stats::poisson(link = "log"))
    start <- p_model$coefficients
    start[is.na(start)] <- 0
    extra <- N_alpha + N_sigma + N_underreport
    start <- c(start, rep_len(.fcr_get_family(family)$start, extra))
    if (N_zero > 0) {
      # Zero model: intercept at a structural-zero share of about 0.2
      # (inflated) or at the observed share of positive counts (hurdle)
      z0 <- rep(0, N_zero)
      int <- which(colnames(X_zero) == "(Intercept)")
      if (length(int)) {
        z0[int] <- if (zero_type == "hurdle") {
          stats::qlogis(min(max(mean(y > 0), 0.01), 0.99))
        } else {
          -1.5
        }
      }
      start <- c(start, z0)
    }
    names(start) <- x_names
  } else {
    start <- .fcr_check_start(start.vals, x_names)
  }
  
  # ---- Log-likelihood ------------------------------------------------------
  local_probFunc <- get_probFunc(family)
  
  logLikFunc <- function(p) {
    coefs <- as.vector(p)
    fixed_coefs <- coefs[seq_len(N_predictors)]
    
    alpha <- NULL
    sigma <- NULL
    
    if (N_alpha > 0) {
      idx <- N_predictors + seq_len(N_alpha)
      alpha <- if (!is.null(mod_alpha_frame)) {
        exp(as.vector(mod_alpha_frame %*% coefs[idx]))
      } else {
        exp(coefs[idx])
      }
    }
    
    if (N_sigma > 0) {
      idx <- N_predictors + N_alpha + seq_len(N_sigma)
      sigma <- if (!is.null(mod_sigma_frame)) {
        exp(as.vector(mod_sigma_frame %*% coefs[idx]))
      } else {
        exp(coefs[idx])
      }
    }
    
    if (N_underreport > 0) {
      underreport_prob <- .report_prob(
        as.vector(X_underreport %*% coefs[idx_under]), underreport_family)
    } else {
      underreport_prob <- 1
    }
    
    predicted <- exp(as.vector(X %*% fixed_coefs) + X_offset) *
      underreport_prob
    
    # Parameter values outside a distribution's support make the C++ code
    # stop(); return NA instead so the optimizer can step back.
    probs <- tryCatch(
      local_probFunc(y = y, predicted = predicted, alpha = alpha,
                     sigma = sigma, haltons = haltons,
                     normed_haltons = normed_haltons, engine = engine),
      error = function(e) rep(NA_real_, length(y)))
    
    if (N_zero > 0) {
      f0 <- tryCatch(
        local_probFunc(y = rep(0, length(y)), predicted = predicted,
                       alpha = alpha, sigma = sigma, haltons = haltons,
                       normed_haltons = normed_haltons, engine = engine),
        error = function(e) rep(NA_real_, length(y)))
      zp <- stats::plogis(as.vector(X_zero %*% coefs[idx_zero]))
      probs <- .fcr_zero_adjust(y, as.vector(probs), as.vector(f0), zp,
                                zero_type)
    }
    
    ll <- log(as.vector(probs))
    ll[weights.df == 0] <- 0
    ll * weights.df
  }
  
  # Analytic gradient for the closed-form Poisson, NB1 and NB2 families
  # (without underreporting or a zero model); numeric otherwise.
  gradFunc <- NULL
  if (family %in% c("POISSON", "NB1", "NB2") && N_underreport == 0 &&
      N_zero == 0) {
    W_alpha <- if (N_alpha == 0) NULL else if (is.null(mod_alpha_frame))
      matrix(1, length(y), 1) else mod_alpha_frame
    gradFunc <- function(p) {
      .fcr_analytic_grad(family, as.vector(p), y, X, X_offset, W_alpha,
                         weights.df)
    }
  }
  
  # ---- Estimation ------------------------------------------------------------
  fit <- maxLik::maxLik(logLikFunc, grad = gradFunc,
                        start = start,
                        method = method, 
                        control = list(iterlim = max.iters, 
                                       printLevel = print.level))
  names(fit$estimate) <- x_names
  
  fit$converged <- .fcr_converged(fit, method)
  if (!isTRUE(fit$converged)) {
    warning("The optimizer did not report convergence (", trimws(fit$message),
            "). Consider a different `method` or more `max.iters`.",
            call. = FALSE)
  }
  
  # ---- Standard errors: one covariance matrix per fit ------------------------
  if (stderr == "boot") {
    boot_data <- data[prep$rows_used, , drop = FALSE]
    boot_w <- if (is.numeric(weights)) prep$weights else weights
    draws <- .fcr_bootstrap(
      n = nrow(boot_data), B = bootstraps, coef_names = x_names,
      refit = function(idx) {
        w_i <- if (is.numeric(boot_w)) boot_w[idx] else boot_w
        m <- countreg(formula = formula, data = boot_data[idx, , drop = FALSE],
                      family = family, offset = offset, weights = w_i,
                      verbose = FALSE,
                      dis_param_formula_1 = dis_param_formula_1,
                      dis_param_formula_2 = dis_param_formula_2,
                      underreport_formula = underreport_formula,
                      underreport_family = underreport_family,
                      zero_formula = zero_formula, zero_type = zero_type,
                      ndraws = ndraws, method = method,
                      max.iters = max.iters, start.vals = fit$estimate,
                      stderr = "normal", engine = engine,
                      keep_data = FALSE)
        m$model$estimate
      })
    V <- stats::cov(draws, use = "pairwise.complete.obs")
    fit$n_boot_ok <- sum(stats::complete.cases(draws))
    vcov_type <- "bootstrap"
  } else if (stderr == "robust") {
    V <- sandwich::sandwich(fit)
    vcov_type <- "robust"
  } else {
    V <- .fcr_hessian_vcov(fit)
    vcov_type <- "hessian"
  }
  dimnames(V) <- list(x_names, x_names)
  
  fit$vcov <- V
  fit$vcov_type <- vcov_type
  fit$stderr <- stderr
  fit$se <- .fcr_se(V)
  names(fit$se) <- x_names
  fit$bootstrapped_se <- if (stderr == "boot") fit$se else NULL
  
  # ---- Store results -------------------------------------------------------
  fit$coefficients <- fit$estimate
  fit$logLik <- fit$maximum
  fit$model <- family
  fit$family <- family
  fit$method <- method
  fit$formula <- formula
  fit$dis_param_formula_1 <- dis_param_formula_1
  fit$dis_param_formula_2 <- dis_param_formula_2
  fit$ndraws <- ndraws
  fit$engine <- engine
  fit$bootstraps <- bootstraps
  fit$offset <- prep$blueprint$offset_label
  fit$modelType <- "countreg"
  fit$underreport_formula <- underreport_formula
  fit$zero_formula <- zero_formula
  fit$zero_type <- if (!is.null(zero_formula)) zero_type
  fit$underreport_family <- underreport_family
  # Sign convention for probit underreporting: 1 = P(reported) = pnorm(z'g)
  # (flexCountReg >= 0.1.3). Fits without this field used 1 - pnorm(z'g).
  fit$underreport_sign <- 1L
  fit$params <- params
  fit$blueprint <- prep$blueprint
  fit$y <- y
  fit$weights <- weights.df
  fit$nobs <- length(y)
  fit$rows_used <- prep$rows_used
  fit$n_dropped <- prep$n_dropped
  
  obj <- .createFlexCountReg(model = fit, data = data, call = match.call(),
                             keep_data = keep_data,
                             formula = formula)
  return(obj)
}
