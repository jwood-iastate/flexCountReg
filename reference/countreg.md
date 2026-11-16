# Count regression models

The purpose of this function is to estimate count regression models
using maximum likelihood estimation (MLE) or Maximum Simulated
Likelihood Estimation (MSLE). The function can estimate the following
models:

- Poisson (Poisson)

- Generalized Poisson Version 1 (GP1)

- Generalized Poisson Version 2 (GP2)

- Negative Binomial 1 (NB1)

- Negative Binomial 2 (NB2)

- Negative Binomial P (NBP)

- Poisson-Lognormal (PLN)

- Poisson-Generalized-Exponential (PGE)

- Poisson-Inverse-Gaussian Type 1 (PIG1)

- Poisson-Inverse-Gaussian Type 2 (PIG2)

- Poisson-Inverse-Gamma (PIG)

- Poisson-Lindley (PL)

- Poisson-Lindley-Gamma (PLG), also known as the Negative
  Binomial-Lindley (NBL)

- Poisson-Lindley-Lognormal (PLL)

- Poisson-Weibull (PW)

- Sichel (SI)

- Generalized Waring (GW)

- Conway-Maxwell-Poisson (COM)

## Usage

``` r
countreg(
  formula,
  data,
  family = "NB2",
  offset = NULL,
  weights = NULL,
  verbose = FALSE,
  dis_param_formula_1 = NULL,
  dis_param_formula_2 = NULL,
  underreport_formula = NULL,
  underreport_family = "logit",
  ndraws = 1500,
  method = "NM",
  max.iters = 1000,
  start.vals = NULL,
  stderr = "normal",
  bootstraps = NULL,
  engine = "halton",
  keep_data = TRUE,
  zero_formula = NULL,
  zero_type = c("inflated", "hurdle")
)
```

## Arguments

- formula:

  a symbolic description of the model to be fitted.

- data:

  a data frame containing the variables in the model.

- family:

  the name of the distribution/model type to estimate. The default "NB2"
  is the standard negative binomial distribution with a log link. other
  options are listed below.

- offset:

  the name of a variable, or vector of variable names, in the data frame
  that should be used as an offset (i.e., included but forced to have a
  coefficient of 1). When several names are given, the offsets are
  summed. Alternatively, an offset can be given in the formula with
  [`offset()`](https://rdrr.io/r/stats/offset.html). Use only one of the
  two approaches: supplying both is an error.

- weights:

  the name of a variable in the data frame, or a numeric vector with one
  value per row of `data`, to be used as frequency weights.

- verbose:

  an optional parameter. If `TRUE`, the function will print out the
  progress of the model fitting. Default is `FALSE`.

- dis_param_formula_1:

  a symbolic description of the model for the natural log of the
  dispersion parameter or first parameter of the count distribution
  used. Further details are provided below.

- dis_param_formula_2:

  a symbolic description of the model for the second parameter of the
  count distribution used. Further details are provided below.

- underreport_formula:

  an optional formula to estimate the underreporting for any of the
  count model options. The underreporting is estimated as a function
  (logit or probit) of the predictors in the model. For the model to be
  tractable, the independent variables cannot be the exact same as the
  count model. The default is `NULL`.

- underreport_family:

  the name of the distribution/model type to estimate the underreporting
  portion of the model when `underreport_formula` is specified. The
  default is "logit" for a binary logistic regression model. The other
  option is "probit" for a probit model.

- ndraws:

  The number of Halton draws for integrating the distribution being
  compounded with the Poisson distribution when there is not a
  closed-form solution. Default is 1500. It is recommended to test
  different numbers of draws to determine if the model is stable (i.e.,
  doesn't change or has minimal change as the number of draws changes
  within a reasonable range).

- method:

  Optimization method to be used for maximum likelihood estimation. One
  of "NM" (Nelder-Mead, the default), "BFGS", "BFGSR", "CG", "NR",
  "BHHH", or "SANN". See `maxLik` documentation for details.

- max.iters:

  Maximum number of iterations for the optimization method.

- start.vals:

  Optional vector of starting values for the optimization. It must have
  one value per parameter. If the vector is named, the values are
  matched to the parameters by name.

- stderr:

  Type of standard errors to use. The default is "normal" (based on the
  Hessian). Other options include "boot" for bootstrapped standard
  errors, or "robust" for robust (sandwich) standard errors.

- bootstraps:

  Optional integer specifying the number of bootstrap samples to be used
  for estimating standard errors when `stderr`= "boot". Defaults to 200
  when `stderr = "boot"` and `bootstraps` is not given.

- engine:

  the engine to use for the integration if the Poisson Lognormal (PLN)
  distribution is used. Options include "halton" which uses Halton draws
  or "poilog" which uses the poilog package. The default is "halton".

- keep_data:

  logical; if `FALSE`, the training data are not stored in the fitted
  object (saving memory).
  [`predict()`](https://rdrr.io/r/stats/predict.html) then needs
  `newdata`, and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)/[`residuals()`](https://rdrr.io/r/stats/residuals.html)
  are unavailable.

- zero_formula:

  an optional one-sided formula for a zero model (`~ 1` for a constant),
  turning the model into a zero-inflated or hurdle model (see
  `zero_type`) for any `family`. The zero model is a logit model. Its
  coefficients are named `Zero:<variable>`. Cannot be combined with
  `underreport_formula`.

- zero_type:

  `"inflated"` (default) or `"hurdle"`. In a zero-inflated model, a
  structural zero occurs with probability \\\pi = logit^{-1}(z\gamma)\\
  and otherwise the count follows the family: \\P(0) = \pi + (1 - \pi)
  f(0)\\ and \\P(y) = (1 - \pi) f(y)\\ for \\y \> 0\\. In a hurdle
  model, \\q = logit^{-1}(z\gamma)\\ is the probability of a positive
  count and positive counts follow the zero-truncated family: \\P(0) =
  1 - q\\ and \\P(y) = q f(y) / (1 - f(0))\\.

## Value

An object of class `flexCountReg` which is a list with the following
components:

- model: the fitted model object.

- data: the data frame used to fit the model.

- call: the matched call.

- formula: the formula used to fit the model.

## Details

For the `family` argument, the following options are available:

- "POISSON" for Poisson distribution with a log link.

- "GP1" for Generalized Poisson Version 1 Distribution with a log link

- "GP2" for Generalized Poisson Version 2 Distribution with a log link

- "NB1" for Negative Binomial 1 distribution with a log link.

- "NB2" for Negative Binomial 2 distribution with a log link (i.e., the
  standard negative binomial model).

- "NBP" for Negative Binomial P distribution with a log link.

- "PLN" for Poisson-Lognormal distribution with a log link.

- "PGE" for Poisson-Generalized-Exponential distribution with a log
  link.

- "PIG1" for Poisson-Inverse-Gaussian Type-1 distribution with a log
  link.

- "PIG2" for Poisson-Inverse-Gaussian Type-2 distribution with a log
  link.

- "PIG" for Poisson-Inverse-Gamma distribution with a log link.

- "PL" for Poisson-Lindley distribution with a log link.

- "PLG" for Poisson-Lindley-Gamma distribution with a log link.

- "PLL" for Poisson-Lindley-Lognormal distribution with a log link.

- "PW" for Poisson-Weibull distribution with a log link.

- "SI" for Sichel distribution with a log link.

- "GW" for Generalized Waring distribution with a log link.

- "COM" for Conway-Maxwell-Poisson (COM) distribution with a log link.

The `dis_param_formula_1` and `dis_param_formula_2` parameters are used
to estimate the dispersion parameter or other parameters of the count
distribution used. This leads to the distributions parameters being
functions rather than constants in the model. For example, if the user
wants to estimate the overdispersion parameter of the Negative Binomial
2 distribution as a function of the variable `x1` and `x2`, the user
would specify `dis_param_formula_1 = ~ x1 + x2`. In the case of the
Negative Binomial distributions, the model is known as a Generalized
Negative Binomial model when the overdispersion parameter is specified
as a function.

The function linking the distribution parameters to the predictors is:
\$\$Param = \exp((Intercept) + \sum \beta X)\$\$

The parameters for the different models are as follows:

For `dis_param_formula_1`, the models are for the parameters:

- \\ln(\phi+1)\\ for the Generalized Poisson Version 1.

- \\ln(\alpha+1)\\ for the Generalized Poisson Version 2.

- \\ln(\alpha)\\ for the Negative Binomial 1 model.

- \\ln(\alpha)\\ for the Negative Binomial 2 model.

- \\ln(\alpha)\\ for the Negative Binomial P model.

- \\ln(\sigma)\\ for the Poisson-Lognormal model.

- \\ln(shape)\\ for the Poisson-Generalized-Exponential model.

- \\ln(\eta)\\ for the Poisson-Inverse-Gaussian model.

- \\ln(\eta)\\ for the Poisson-Inverse-Gamma model.

- \\ln(\theta)\\ for the Poisson-Lindley model.

- \\ln(\theta)\\ for the Poisson-Lindley-Gamma model.

- \\ln(\theta)\\ for the Poisson-Lindley-Lognormal model.

- \\ln(\alpha)\\ for the Poisson-Weibull model.

- \\\gamma\\ for the Sichel model. This parameter can be negative, so it
  is modeled directly (identity link), \\\gamma = W\delta\\, rather than
  through \\\exp()\\.

- \\ln(k)\\ for the Generalized Waring model.

- \\ln(\nu)\\ for the Conway-Maxwell-Poisson model.

For `dis_param_formula_2`, the models are for the parameters:

- Not Applicable for the Generalized Poisson Version 1 model.

- Not Applicable for the Generalized Poisson Version 2 model.

- Not Applicable for the Negative Binomial 1 model.

- Not Applicable for the Negative Binomial 2 model.

- \\ln(p)\\ for the Negative Binomial P model.

- Not Applicable for the Poisson-Lognormal model.

- \\ln(scale)\\ for the Poisson-Generalized-Exponential model.

- Not Applicable for the Poisson-Inverse-Gaussian model.

- Not Applicable for the Poisson-Inverse-Gamma model.

- Not Applicable for the Poisson-Lindley model.

- \\ln(\alpha)\\ for the Poisson-Lindley-Gamma model.

- \\ln(\sigma)\\ for the Poisson-Lindley-Lognormal model.

- \\ln(\sigma)\\ for the Poisson-Weibull model.

- \\ln(\sigma)\\ for the Sichel model.

- \\ln(\rho-1)\\ for the Generalized Waring model.

- Not Applicable for the Conway-Maxwell-Poisson model.

The `ndraws` parameter is used to estimate the distribution when there
is not a closed-form solution. This uses Halton draws to integrate the
distribution being compounded with the Poisson distribution. The default
is 1500. The models this is applicable for include:

- Poisson-Lognormal

- Poisson-Generalized-Exponential

- Poisson-Lindley-Gamma (more efficient than using hypergeometric
  functions)

- Poisson-Lindley-Lognormal

- Poisson-Weibull

## Model Details

### Poisson Model

This implements the Poisson regression model using Maximum Likelihood
Estimation (MLE), as opposed to the Iteratively Reweighted Least Squares
(IRLS) method used in the `glm` function.

The PMF and log-likelihood functions are: \$\$P(Y = y) = \frac{e^{-\mu}
\mu^y}{y!}\$\$ \$\$LL\_{\text{Poisson}}(\beta) = \sum\_{i=1}^n \left\[
-\mu_i + y_i \ln(\mu_i) - \ln(y_i!) \right\]\$\$

The mean is: \$\$\mu = exp(X\beta)\$\$

The variance is: \$\$\text{Var}(Y) = \mu\$\$

### Generalized Poisson Version 1 Model

This implements the Generalized Poisson Version 1 model using MLE.

The PMF is: \$\$ f(y\|\phi,\mu)=\frac{\mu(\mu+\phi y)^{y-1}
exp\left(-\frac{\mu+\phi y} {1+\phi}\right)}{(1+\phi)^y y!} \$\$

Where \\\phi\\ is a dispersion parameter with the restriction that
\\\phi\>-1\\, \\\mu\>0\\ is the mean value, and \\y\\ is a non-negative
integer. The mean is modeled as: \$\$\mu = exp(X\beta)\$\$

The variance of the GP-1 distribution is: \$\$\sigma^2=(1+\phi)^2
\mu\$\$

If \\\phi\>0\\, the distribution is overdispersed. If \\\phi=0\\, the
distribution is equidispersed. If \\\phi\<0\\, the distribution is
underdispersed.

Furthermore, \\\phi\>-1\\ is required for this distribution. When
\\\phi\<0\\, there is also a maximum value of support for the integer
value \\y\\. This is \\y_max=\left\lfloor -\frac{\mu}{\phi}
\right\rfloor\\. To account for the range of possible values, this is
modeled as \$\$ln(\phi+1)\$\$.

### Generalized Poisson Version 2 Model

This implements the Generalized Poisson Version 2 model using MLE.

The PMF is: \$\$ f(y\|\mu,\alpha)= \frac{\mu(\mu+\alpha\mu y)^{y-1}
\exp\left(-\frac{\mu+\alpha\mu y}{1+\alpha\mu}\right)} {(1+\alpha\mu)^y
y!} \$\$

Where \\\alpha\\ is a dispersion parameter with the restriction that
\\1+\alpha\mu\>0\\, \\\mu\>0\\ is the mean value, and \\y\\ is a
non-negative integer. The mean is modeled as: \$\$\mu = exp(X\beta)\$\$

The variance of the GP-2 distribution is:
\$\$\sigma^2=\mu(1+\alpha\mu)^2\$\$

If \\\alpha \> 0\\, the distribution is overdispersed. If \\\alpha =
0\\, the distribution is equidispersed and reduces to the ordinary
Poisson distribution. If \\\alpha \< 0\\ and \\1 + \alpha\mu \> 0\\, the
distribution is underdispersed and has finite support.

Under the underdispersed case, the largest possible value of \\y\\ is
bounded above by \\\left\lceil -1/\alpha \right\rceil - 1\\, i.e., the
largest integer satisfying \\1 + \alpha y \> 0\\. To account for the
range of possible values, this is modeled as \$\$ln(\alpha+1)\$\$.

### Negative Binomial Models

The NB-1, NB-2, and NB-P versions of the negative binomial distribution
are based on Greene (2008). The details of each of these are provided
below.

#### NB-1 Model

The PMF and log-likelihood functions are: \$\$P(Y = y) =
\frac{\Gamma(y + \frac{\mu}{\alpha})}{y! \\ \Gamma(\frac{\mu}{\alpha})}
\left( \frac{\frac{\mu}{\alpha}}{\frac{\mu}{\alpha} + \mu} \right)
^{\frac{\mu}{\alpha}} \left( \frac{\mu}{\frac{\mu}{\alpha} + \mu}
\right)^y\$\$ \$\$LL\_{\text{NB1}}(\beta, \alpha) = \sum\_{i=1}^n
\left\[ \ln \Gamma\left( y_i + \frac{\mu_i}{\alpha} \right) - \ln
\Gamma\left( \frac{\mu_i}{\alpha} \right) - \ln y_i! +
\frac{\mu_i}{\alpha} \ln \left(
\frac{\frac{\mu_i}{\alpha}}{\frac{\mu_i}{\alpha} + \mu_i} \right) + y_i
\ln \left( \frac{\mu_i}{\frac{\mu_i}{\alpha} + \mu_i} \right)
\right\]\$\$

The mean is: \$\$\mu = exp(X\beta)\$\$

The variance is: \$\$\text{Var}(Y) = \mu + \alpha\mu\$\$

#### NB-2 Model

The PMF and log-likelihood functions are: \$\$P(Y = y) =
\frac{\Gamma(y + \alpha)}{y! \\ \Gamma(\alpha)} \left(
\frac{\alpha}{\alpha + \mu} \right)^\alpha \left( \frac{\mu}{\alpha +
\mu} \right)^y\$\$ \$\$LL\_{\text{NB2}} = \sum\_{i=1}^n \left\[ \ln
\Gamma(y_i + \alpha) - \ln \Gamma(\alpha) - \ln y_i! + \alpha \ln \left(
\frac{\alpha}{\alpha + \mu_i} \right) + y_i \ln \left(
\frac{\mu_i}{\alpha + \mu_i} \right) \right\]\$\$

The mean is: \$\$\mu = exp(X\beta)\$\$

The variance is: \$\$\text{Var}(Y) = \mu + \alpha\mu^2\$\$

#### NB-P Model

The PMF and log-likelihood functions are: \$\$P(Y = y) =
\frac{\Gamma(y + \frac{\mu^{2-p}}{\alpha})}{y! \\
\Gamma(\frac{\mu^{2-p}}{\alpha})} \left( \frac{\frac{\mu^{2-p}}{\alpha}}
{\frac{\mu^{2-p}}{\alpha} + \mu} \right)^{\frac{\mu^{2-p}}{\alpha}}
\left( \frac{\mu}{\frac{\mu^{2-p}}{\alpha} + \mu} \right)^y\$\$
\$\$LL\_{\text{NBP}}(\beta, \alpha, p) = \sum\_{i=1}^n \left\[ \ln
\Gamma\left( y_i + \frac{\mu_i^{2-p}}{\alpha} \right) - \ln \Gamma\left(
\frac{\mu_i^{2-p}}{\alpha} \right) - \ln y_i! +
\frac{\mu_i^{2-p}}{\alpha} \ln \left(
\frac{\frac{\mu_i^{2-p}}{\alpha}}{\frac{\mu_i^{2-p}}{\alpha} + \mu_i}
\right) + y_i \ln \left( \frac{\mu_i}{\frac{\mu_i^{2-p}}{\alpha} +
\mu_i} \right) \right\]\$\$

The mean is: \$\$\mu = exp(X\beta)\$\$

The variance is: \$\$\text{Var}(Y) = \mu + \alpha\mu^p\$\$

### Poisson-Lognormal (PLN) Model

The compound Probability Mass Function(PMF) for the Poisson-Lognormal
distribution is: \$\$f(y\|\lambda,\sigma)=\int_0^\infty \frac{\lambda^y
x^y e^{-\lambda x}} {y!}\frac{exp\left(-\frac{ln^2(x)}{2\sigma^2}
\right)}{x\sigma\sqrt{2\pi}}dx\$\$

Where \\\sigma\\ is a parameter for the lognormal distribution with the
restriction \\\sigma\>0\\, and \\y\\ is a non-negative integer.

The expected value of the distribution is:
\$\$E\[y\]=e^{X\beta+\sigma^2/2} = \mu e^{\sigma^2/2}\$\$

The t-statistics and p-values for the coefficients related to ln(sigma)
are, by default, testing if the coefficients are different from a value
of 0. This has little practical meaning given that they are coefficients
for ln(sigma). They are not testing if the coefficients have statistical
significance in terms of improvement over a Poisson model. To test
whether the Poisson-Lognormal model fits significantly better than the
Poisson model, use a likelihood-ratio test with
[`regCompTest`](https://jwood-iastate.github.io/flexCountReg/reference/regCompTest.md).

### Poisson Generalized-Exponential (PGE) Model

The Generalized Exponential distribution can be written as a function
with a shape parameter \\\alpha\>0\\ and scale parameter \\\gamma\>0\\.
The distribution has strictly positive continuous values. The PDF of the
distribution is:
\$\$f(x\|\alpha,\gamma)=\frac{\alpha}{\gamma}\left(1-e^{-\frac{x}{\gamma}}
\right)^{\alpha-1}e^{-\frac{x}{\gamma}}\$\$

Thus, the compound Probability Mass Function(PMF) for the PGE
distribution is: \$\$f(y\|\lambda,\alpha,\beta)=\int_0^\infty
\frac{\lambda^y x^y e^{-\lambda
x}}{y!}\frac{\alpha}{\gamma}\left(1-e^{-\frac{x}{\gamma}}
\right)^{\alpha-1}e^{-\frac{x}{\gamma}} dx\$\$

The expected value of the distribution is: \$\$E\[y\]=\mu=\lambda
\left(\frac{\psi(\alpha+1)-\psi(1)}{\gamma}\right)\$\$

Where \\\psi(\cdot)\\ is the digamma function.

The variance is: \$\$\sigma^2=\lambda
\left(\frac{\psi(\alpha+1)-\psi(1)}{\gamma}\right) +
\left(\frac{-\psi'(\alpha+1)+\psi'(1)}{\gamma^2}\right)\lambda^2\$\$

Where \\\psi'(\cdot)\\ is the trigamma function.

To ensure that \\\mu=e^{X\beta}\\, \\\lambda\\ is replaced with:
\$\$\lambda=\frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)}\$\$

This results in: \$\$f(y\|\mu,\alpha,\beta)=\int_0^\infty
\frac{\left(\frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)}\right)^y
x^y e^{-\left(\frac{\gamma e^{X\beta}}{\psi(\alpha+1)-\psi(1)}\right)
x}}{y!}\frac{\alpha}{\gamma}
\left(1-e^{-\frac{x}{\gamma}}\right)^{\alpha-1}e^{-\frac{x}{\gamma}}
dx\$\$

Halton draws are used to perform simulation over the lognormal
distribution to solve the integral.

### Poisson-Inverse-Gaussian Type 1 (PIG1) and Type 2 (PIG2) Models

The Poisson-Inverse-Gaussian regression model is based on the
Poisson-Inverse-Gaussian Distribution.

The expected value of the distribution in the regression utilizes a
log-link function. Thus, the mean is: \$\$\mu=e^{X\beta}\$\$

The variance function for the Type 1 distribution (which is the default)
is: \$\$\sigma^2=\mu+\eta\mu\$\$

While the variance for the Type 2 distribution is:
\$\$\sigma^2=\mu+\eta\mu^2\$\$

The parameter \\\eta\\ is estimated as the natural logarithm transformed
value, \\\ln(\eta)\\, to ensure that \\\eta\>0\\.

### Poisson-Inverse-Gamma (PIG) Model

The PDF of the distribution is:
\$\$f(x\|\eta,\mu)=\frac{2\left(\mu\left(\frac{1}{\eta}+1\right)\right)
^{\frac{x+\frac{1}{\eta}+2}{2}}}{x!\Gamma\left(\frac{1}{\eta}+2\right)}
K\_{x-\frac{1}{\eta}-2}\left(2\sqrt{\mu\left(\frac{1}{\eta}+1\right)}\right)\$\$

Where \\\eta\\ is a shape parameter with the restriction that
\\\eta\>0\\, \\\mu\>0\\ is the mean value, \\y\\ is a non-negative
integer, and \\K_i(z)\\ is the modified Bessel function of the second
kind. This formulation uses the mean directly.

The variance of the distribution is: \$\$\sigma^2=\mu+\eta\mu^2\$\$

### Poisson-Lindley (PL) Model

The Poisson-Lindley regression is based on a compound Poisson-Lindley
distribution. It handles count outcomes with high levels of zero
observations (or other high densities at low outcome values) that
standard count regression methods, including the negative binomial, may
struggle to adequately capture or model.

The compound Probability Mass Function(PMF) for the Poisson-Lindley (PL)
distribution is:
\$\$f(y\|\theta,\lambda)=\frac{\theta^2\lambda^y(\theta+\lambda+
y+1)}{(\theta+1)(\theta+\lambda)^{y+2}}\$\$

Where \\\theta\\ and \\\lambda\\ are distribution parameters with the
restrictions that \\\theta\>0\\ and \\\lambda\>0\\, and \\y\\ is a
non-negative integer.

The expected value of the distribution is:
\$\$\mu=\frac{\lambda(\theta+2)}{\theta(\theta+1)}\$\$

If a log-link function is used, the mean is:
\$\$\mu=e^{X\beta}=\frac{\lambda(\theta+2)}{\theta(\theta+1)}\$\$

Thus, the parameter \\\lambda\\ in the PL distribution when applied to
regression analysis is:
\$\$\lambda=\frac{\mu\theta(\theta+1)}{\theta+2}\$\$

Using the replacement and simplifying results in: \$\$f(y \mid \theta,
\mu) = \\ \frac{\theta^2 (\mu \theta (\theta+1))^y (\theta^2 (1+\mu) +
\theta (2+\mu) + (\theta+2) (y+1))}{(\theta+1) (\theta+2)^{y+1}
(\theta^2 (1+\mu) + \theta (2+\mu))^{y+2}}\$\$ And \$\$LL=2
\log(\theta) + y (\log(\mu) + \log(\theta) + \log(\theta+1)) +
\log(\theta^2 (1+\mu) + \theta (2+\mu) + (\theta+2) (y+1)) -
\log(\theta+1) - (y+1) \log(\theta+2) - (y+2) \log(\theta^2 (1+\mu) +
\theta (2+\mu))\$\$

The variance function is defined as:
\$\$\sigma^2=\mu+\left(1-\frac{2}{(\theta+2)^2}\right)\mu^2\$\$

It should be noted that the p-value for the parameter `ln(theta)` in the
model summary is testing if the parameter `theta` is equal to a value
of 1. This has no practical meaning. The Likelihood-Ratio (LR) test
compares the Poisson-Lindley regression with a Poisson regression with
the same independent variables. Thus, the LR test result indicates the
statistical significance for the improvement in how well the model fits
the data over a Poisson regression. This indicates the statistical
significance of the `theta` parameter.

### Poisson-Lindley-Gamma (PLG) Model

The Poisson-Lindley-Gamma regression is based on a compound
Poisson-Lindley- Gamma distribution. Details of the distribution can be
seen at
[`dplindGamma`](https://jwood-iastate.github.io/flexCountReg/reference/NegativeBinomialLindley.md).

The mean for the regression model is: \$\$\mu=e^{X\beta}\$\$

The variance function is defined as:
\$\$\sigma^2=\mu+\left(2\alpha+1-\frac{2(1+\alpha)}
{(\theta+2)^2}\right)\mu^2\$\$

It should be noted that the p-value for the parameters `ln(theta)` and
`ln(alpha)` in the model summary are testing if the parameter `theta`
and `alpha` are equal to a value of 1.

### Poisson-Lindley-Lognormal (PLL) Model

The Poisson-Lindley-Lognormal regression is based on a compound Poisson-
Lindley-Lognormal distribution. Details of the distribution can be seen
at
[`dplindLnorm`](https://jwood-iastate.github.io/flexCountReg/reference/PoissonLindleyLognormal.md).

The mean for the regression model is: \$\$\mu=e^{X\beta}\$\$

The variance function is defined as:
\$\$\sigma^2=\mu+\left(\frac{1-\frac{2}{(\theta+2)^2}}{e^{\frac{\sigma^2}
{2}}}+e^{\sigma^2}-1\right)\mu^2\$\$

It should be noted that the p-value for the parameters `ln(theta)` and
`ln(sigma)` in the model summary are testing if the parameter `theta`
and `sigma` are equal to a value of 1.

### Poisson-Weibull (PW) Model

The Poisson-Weibull distribution uses the Weibull distribution as a
mixing distribution for a Poisson process. It is useful for modeling
overdispersed count data. The density function (probability mass
function) for the Poisson-Weibull distribution is given by:
\$\$P(y\|\lambda,\alpha,\sigma) = \int_0^\infty \frac{e^{-\lambda x}
\lambda^y x^y }{y!}
\left(\frac{\alpha}{\sigma}\right)\left(\frac{x}{\sigma}
\right)^{\alpha-1}e^{-\left(\frac{x}{\sigma}\right)^\alpha} dx\$\$ where
\\f(x\| \alpha, \sigma)\\ is the PDF of the Weibull distribution and
\\\lambda\\ is the mean of the Poisson distribution.

For the Poisson-Weibull Regression model, the expected values is:
\$\$E\[Y\]=\lambda\sigma\Gamma\left(1+\frac{1}{\alpha}\right)\$\$ Where
\\\lambda\\ is the mean of the Poisson distribution, \\\alpha\\ is the
shape parameter, and \\\sigma\\ is the scale parameter.

To ensure that the regression model predicts the mean value, the
regression utilizes:
\$\$\mu=\exp{X\gamma}=\lambda\sigma\Gamma\left(1+\frac{1}{\alpha}\right)\$\$
Where \\X\\ is a matrix of independent variables and \\\gamma\\ is a
vector of coefficients.

This leads to:
\$\$\lambda=\frac{\mu}{\sigma\Gamma\left(1+\frac{1}{\alpha}\right)}\$\$

The variance for the Poisson-Weibull regression is:
\$\$V\[Y\]=\mu+\left(\frac{\Gamma\left(1+\frac{2}{\alpha}\right)}
{\Gamma\left(1+\frac{1}{\alpha}\right)^2}-1\right)\mu^2\$\$

### Sichel (SI) Model

The compound Probability Mass Function (PMF) for the Sichel distribution
uses the formulation from Zou, Lord and Zhang (2012) and Rigby et al.
(2008): \$\$f(y\|\mu, \sigma, \gamma)=\frac{\left(\frac{\mu}{c}\right)^y
K\_{y+\gamma}(\alpha)}{K\_\gamma(1/\sigma)y!(\alpha\sigma)^{y+\gamma}}\$\$

Where \\\sigma\\ and \\\gamma\\ are distribution parameters with
\\\sigma\>0\\,
\$\$c=\frac{K\_{\gamma+1}(1/\sigma)}{K\_\gamma(1/\sigma)}\$\$,
\$\$\alpha^2=\sigma^{-2}+2\mu(c\sigma)^{-1}\$\$, a mean value of
\\\mu\\, \\y\\ is a non-negative integer, and \$\$K_j(x)\$\$ is a
modified Bessel function of the third kind with order \$\$j\$\$ and
argument \$\$x\$\$.

The variance of the distribution is:
\$\$\sigma^2=\mu+\left(\frac{2\sigma(\gamma+1)}{c}+\frac{1}{c^2}-1\right)
\mu^2\$\$

### Generalized Waring (GW) Model

The following are the PMF, mean, and variance used for the Generalized
Waring model. The regression formulation replaces parameter \\\alpha\\
with the mean \\\mu\\, while retaining parameters \\k\\ and \\\rho\\.

In the original parameterization, the PMF is: \$\$ P(Y=y)= \frac{
\Gamma(\alpha+y)\Gamma(k+y)\Gamma(\rho+k)\Gamma(\alpha+\rho) }{
y!\Gamma(\alpha)\Gamma(k)\Gamma(\rho) \Gamma(\alpha+k+\rho+y) }, \$\$
where \\y=0,1,2,\ldots\\ and \\\alpha,k,\rho\>0\\.

The mean and variance are: \$\$ E(Y)=\frac{\alpha k}{\rho-1}, \qquad
\rho\>1, \$\$ \$\$ \operatorname{Var}(Y)= \frac{\alpha
k(\alpha+\rho-1)(k+\rho-1)} {(\rho-1)^2(\rho-2)}, \qquad \rho\>2. \$\$

In accident applications, the variance can be decomposed into components
attributed to randomness, proneness, and liability. These are variance
components, rather than individual interpretations of \\\alpha\\, \\k\\,
and \\\rho\\.

To parameterize the model in terms of its mean, set: \$\$
\mu=e^{X\beta}, \qquad \alpha=\frac{\mu(\rho-1)}{k}. \$\$

Substituting into the PMF gives: \$\$ P(Y=y \mid X)= \frac{
\Gamma\left(\frac{\mu(\rho-1)}{k}+y\right) \Gamma(k+y)\Gamma(\rho+k)
\Gamma\left(\frac{\mu(\rho-1)}{k}+\rho\right) }{
y!\Gamma\left(\frac{\mu(\rho-1)}{k}\right) \Gamma(k)\Gamma(\rho)
\Gamma\left(\frac{\mu(\rho-1)}{k}+k+\rho+y\right) }. \$\$

This results in a regression model with: \$\$ E(Y \mid
X)=\mu=e^{X\beta}, \$\$ \$\$ \operatorname{Var}(Y \mid X)=
\frac{\mu(k+\mu)(k+\rho-1)}{k(\rho-2)} = \frac{k+\rho-1}{\rho-2}\mu+
\frac{k+\rho-1}{k(\rho-2)}\mu^2, \qquad \rho\>2. \$\$

The original distribution is valid for all \\\alpha,k,\rho\>0\\,
including \\\rho=1\\ and \\\rho=2\\. Its mean is infinite when
\\0\<\rho\leq1\\. When \\1\<\rho\leq2\\, the mean is finite but the
variance is infinite. A finite variance requires \\\rho\>2\\. The
mean-based regression parameterization requires \\\mu\>0\\, \\k\>0\\,
and \\\rho\>1\\.

### Conway-Maxwell-Poisson (COM) Model

The following is the the PMF for the COM model. \$\$f(x\|\lambda,
\nu)=\frac{\lambda^x}{(x!)^\nu Z(\lambda,\nu)}\$\$

Where \\\lambda\\ and \\\nu\\ are distribution parameters with
\\\lambda\>0\\ and \\\nu\>0\\, and \\Z(\lambda,\nu)\\ is the normalizing
constant.

The normalizing constant is given by:
\$\$Z(\lambda,\nu)=\sum\_{n=0}^{\infty}\frac{\lambda^n}{(n!)^\nu}\$\$

The mean and variance are: \$\$\mu=e^{X\beta}=\lambda
\frac{\delta}{\delta \lambda} \log(Z(\lambda, \nu))\$\$
\$\$\sigma^2=\lambda \frac{\delta}{\delta \lambda} \mu\$\$

Note that the COM distribution parameter \\\lambda\\ is solved for using
\\\mu\\ and \\\nu\\, so the regression model provides direct predictions
for the mean.

### Underreporting

Models for underreporting combine a binary probability model (logit or
probit) with a count model. This is accomplished using a model for the
probability of crashes being reported multiplied by the estimated mean
for the count model, based on the observed data. This is discussed in
Wood et. al. (2016), Pararai et. al., (2006), and Pararai et. al.,
(2010). The underreporting model is based on:
\$\$\mu\_{observed}=\mu\_{true}\cdot P(\text{event is reported})\$\$

This allows the inference of both the true event count and the
probability of the event being reported as a function of independent
variables.

## References

Greene, W. (2008). Functional forms for the negative binomial model for
count data. Economics Letters, 99(3), 585-590.

Pararai, M., Famoye, F., & Lee, C. (2006). Generalized Poisson
regression model for underreported counts. Advances and applications in
Statistics, 6(3), 305-322.

Pararai, M., Famoye, F., & Lee, C. (2010). Generalized Poisson-Poisson
mixture model for misreported counts with an application to smoking
data. Journal of Data Science, 8(4), 607-617.

Rigby, R. A., Stasinopoulos, D. M., & Akantziliotou, C. (2008). A
framework for modelling overdispersed count data, including the
Poisson-shifted generalized inverse Gaussian distribution. Computational
Statistics & Data Analysis, 53(2), 381-393.

Wood, J.S., Eric T. Donnell, & Christopher J. Fariss. "A method to
account for and estimate underreporting in crash frequency research."
Accident Analysis & Prevention 95 (2016): 57-66.

Zou, Y., Lord, D., & Zhang, Y. (2012). Analyzing highly dispersed crash
data using the Sichel generalized additive models for location, scale
and shape. Accident Analysis & Prevention, 45.

## Examples

``` r
# \donttest{
# Load the Washington data
data("washington_roads")
washington_roads$AADT10kplus <- ifelse(washington_roads$AADT > 10000, 1, 0)
# Estimate an NB2 model with a dispersion parameter as a function of the 
# variable `speed50` (i.e., generalized NB2), verbose output, and use the 
# BFGS optimization method
nb2 <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
               data = washington_roads, family = "NB2",
               dis_param_formula_1 = ~ speed50, verbose = TRUE, 
               method='BFGS')
#> initial  value 1087.972055 
#> iter   2 value 1086.761180
#> iter   3 value 1082.460817
#> iter   4 value 1081.386383
#> iter   5 value 1080.983877
#> iter   6 value 1079.708975
#> iter   7 value 1070.603399
#> iter   8 value 1068.239211
#> iter   9 value 1065.881517
#> iter  10 value 1065.612937
#> iter  11 value 1065.068995
#> iter  12 value 1064.893660
#> iter  13 value 1064.876776
#> iter  14 value 1064.876161
#> iter  15 value 1064.876105
#> iter  15 value 1064.876105
#> iter  15 value 1064.876105
#> final  value 1064.876105 
#> converged
summary(nb2)
#> Call:
#> countreg(formula = Total_crashes ~ lnaadt + lnlength + speed50 + 
#>     AADT10kplus, data = washington_roads, family = "NB2", verbose = TRUE, 
#>     dis_param_formula_1 = ~speed50, method = "BFGS")
#> 
#> Model: countreg | family: NB2 | optimizer: BFGS
#> Iterations: 44 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>                       Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)            -7.4015     0.4729  -15.65  < 2e-16 ***
#> lnaadt                  0.9116     0.0566   16.09  < 2e-16 ***
#> lnlength                0.8428     0.0700   12.04  < 2e-16 ***
#> speed50                -0.4699     0.1151   -4.08  4.5e-05 ***
#> AADT10kplus             0.7697     0.1299    5.92  3.1e-09 ***
#> ln(alpha):(Intercept)  -1.6192     0.3740   -4.33  1.5e-05 ***
#> ln(alpha):speed50       1.3059     0.5963    2.19    0.029 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1064.876 | AIC: 2143.752 | BIC: 2180.949 
#> N: 1501
#> LR test vs Poisson: chi-sq = 22.092 on 2 df, p-value = 1.59e-05 


# Estimate a Poisson-Lognormal model (integration with the poilog package)
if (requireNamespace("poilog", quietly = TRUE)) {
  pln <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                  data = washington_roads, family = "PLN", engine="poilog")
  summary(pln)
}
#> Call:
#> countreg(formula = Total_crashes ~ lnaadt + lnlength + speed50 + 
#>     AADT10kplus, data = washington_roads, family = "PLN", engine = "poilog")
#> 
#> Model: countreg | family: PLN | optimizer: NM
#> Iterations: 733 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   -7.507         NA      NA       NA    
#> lnaadt         0.910         NA      NA       NA    
#> lnlength       0.847      0.110    7.70  1.4e-14 ***
#> speed50       -0.459      0.105   -4.37  1.2e-05 ***
#> AADT10kplus    0.809      0.133    6.09  1.1e-09 ***
#> ln(sigma)     -0.734         NA      NA       NA    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1066.67 | AIC: 2145.341 | BIC: 2177.224 
#> N: 1501
#> LR test vs Poisson: chi-sq = 18.504 on 1 df, p-value = 8.48e-06 (boundary-corrected) 

# Estimate a Negative Binomial with underreporting (probit)
plogn_underreport <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + 
               AADT10kplus,
              data = washington_roads, family = "NB2",
              underreport_formula = ~ speed50 + AADT10kplus, 
              underreport_family = "probit")
#> Warning: The optimizer did not report convergence (iteration limit exceeded). Consider a different `method` or more `max.iters`.
summary(plogn_underreport)
#> Call:
#> countreg(formula = Total_crashes ~ lnaadt + lnlength + speed50 + 
#>     AADT10kplus, data = washington_roads, family = "NB2", underreport_formula = ~speed50 + 
#>     AADT10kplus, underreport_family = "probit")
#> 
#> Model: countreg | family: NB2 | optimizer: NM
#> Iterations: 1002 
#> Convergence: iteration limit exceeded 
#> Warning: the optimizer did not report convergence.
#> 
#> Coefficients (standard errors: Hessian-based):
#>                            Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)                 -6.6903     1.1072   -6.04  1.5e-09 ***
#> lnaadt                       0.9189     0.0551   16.68  < 2e-16 ***
#> lnlength                     0.8470     0.0700   12.11  < 2e-16 ***
#> speed50                      0.0534     0.4474    0.12     0.90    
#> AADT10kplus                 -0.0274     1.0722   -0.03     0.98    
#> ln(alpha)                   -1.4959     0.3417   -4.38  1.2e-05 ***
#> Underreporting:(Intercept)  -0.0898     1.1909   -0.08     0.94    
#> Underreporting:speed50      -0.5394     0.6136   -0.88     0.38    
#> Underreporting:AADT10kplus   2.5350     4.5067    0.56     0.57    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1065.202 | AIC: 2148.404 | BIC: 2196.229 
#> N: 1501

# Estimate a Conway-Maxwell-Poisson model
com_model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + 
                      AADT10kplus,
                      data = washington_roads, family = "COM", method="BHHH")
summary(com_model)
#> Call:
#> countreg(formula = Total_crashes ~ lnaadt + lnlength + speed50 + 
#>     AADT10kplus, data = washington_roads, family = "COM", method = "BHHH")
#> 
#> Model: countreg | family: COM | optimizer: BHHH
#> Iterations: 6 
#> Convergence: successive function values within relative tolerance limit (reltol) 
#> 
#> Coefficients (standard errors: Hessian-based):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -7.4214     0.4652  -15.95  < 2e-16 ***
#> lnaadt        0.9124     0.0562   16.23  < 2e-16 ***
#> lnlength      0.8399     0.0667   12.59  < 2e-16 ***
#> speed50      -0.4451     0.1031   -4.32  1.6e-05 ***
#> AADT10kplus   0.7875     0.1259    6.25  4.0e-10 ***
#> ln(nu)       -0.6194     0.1714   -3.61    3e-04 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1066.063 | AIC: 2144.126 | BIC: 2176.009 
#> N: 1501
#> LR test vs Poisson: chi-sq = 19.719 on 1 df, p-value = 8.97e-06 

# Estimate a Generalized-Poisson Version 1 model
gp1_model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + 
                      AADT10kplus,
                      data = washington_roads, family = "GP1", method="BHHH")
summary(gp1_model)
#> Call:
#> countreg(formula = Total_crashes ~ lnaadt + lnlength + speed50 + 
#>     AADT10kplus, data = washington_roads, family = "GP1", method = "BHHH")
#> 
#> Model: countreg | family: GP1 | optimizer: BHHH
#> Iterations: 6 
#> Convergence: successive function values within relative tolerance limit (reltol) 
#> 
#> Coefficients (standard errors: Hessian-based):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -7.2956     0.4747  -15.37  < 2e-16 ***
#> lnaadt        0.8963     0.0565   15.86  < 2e-16 ***
#> lnlength      0.8200     0.0616   13.30  < 2e-16 ***
#> speed50      -0.4448     0.0976   -4.56  5.2e-06 ***
#> AADT10kplus   0.7684     0.1093    7.03  2.0e-12 ***
#> ln(phi+1)     0.0875     0.0277    3.16   0.0016 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1069.222 | AIC: 2150.445 | BIC: 2182.328 
#> N: 1501
#> LR test vs Poisson: chi-sq = 13.4 on 1 df, p-value = 0.000252 

# Estimate a Generalized-Poisson Version 2 model
gp2_model <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + 
                      AADT10kplus,
                      data = washington_roads, family = "GP2", method="BHHH")
summary(gp2_model)
#> Call:
#> countreg(formula = Total_crashes ~ lnaadt + lnlength + speed50 + 
#>     AADT10kplus, data = washington_roads, family = "GP2", method = "BHHH")
#> 
#> Model: countreg | family: GP2 | optimizer: BHHH
#> Iterations: 7 
#> Convergence: successive function values within relative tolerance limit (reltol) 
#> 
#> Coefficients (standard errors: Hessian-based):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -7.3820     0.4606  -16.03  < 2e-16 ***
#> lnaadt        0.9088     0.0558   16.28  < 2e-16 ***
#> lnlength      0.8470     0.0692   12.23  < 2e-16 ***
#> speed50      -0.4599     0.1027   -4.48  7.5e-06 ***
#> AADT10kplus   0.8123     0.1350    6.02  1.8e-09 ***
#> ln(alpha+1)   0.1056     0.0284    3.72    2e-04 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1066.757 | AIC: 2145.513 | BIC: 2177.397 
#> N: 1501
#> LR test vs Poisson: chi-sq = 18.331 on 1 df, p-value = 1.86e-05 
# }
```
