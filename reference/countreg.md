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

- Discrete Weibull (DW)

- Double Poisson (DPO)

- Delaporte (DELAPORTE)

- Hermite of order 2, 3 or 4 (HERMITE, HERMITE3, HERMITE4)

- Hyper-Poisson (HYPERPOIS)

- Gamma renewal count (GAMMACOUNT)

- Polya-Aeppli, i.e. Poisson-Tweedie with power 1.5 (PA)

- Negative binomial-generalized exponential (NBGE)

See
[count_families](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md)
for a table of all families with their `family` codes, distribution
parameters, variances, distribution functions and available model
features.

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
  engine = "quadrature",
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
  is the standard negative binomial distribution with a log link. See
  [count_families](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md)
  for all families, their parameters and the features available for
  each.

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

  an optional formula for the probability that a true count is reported
  (binomial thinning of the true counts), as a logit or probit function
  of its predictors. Available for the families whose reported counts
  follow the same family (see
  [count_families](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md);
  for NB1 and PIG1 the dispersion becomes \\p\alpha\\, for NBP \\\alpha
  p^{2-P}\\, for HERMITE \\\rho\\ becomes \\p\rho\\ and for PA \\\phi\\
  becomes \\\sqrt{p}\phi\\). It is identified only by variables that
  affect reporting but not the true counts (and, weakly, by the shape of
  the link function): an intercept-only formula is an error, and a
  warning is given when the reporting and count coefficients are nearly
  collinear. Cannot be combined with `zero_formula`. The default is
  `NULL`.

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

  the numerical method for the mixed Poisson families PLN, PGE, PW and
  PLL: `"quadrature"` (default; adaptive quadrature on the log scale,
  accurate for every count), `"halton"` (the average over Halton draws
  used by flexCountReg 0.1.2 and earlier, kept to reproduce earlier
  fits; it is badly biased for counts far in the tail) or, for PLN only,
  `"poilog"` (the poilog package).

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

The `family` codes are listed in
[count_families](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md).
Family names are not case-sensitive.

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
\$\$Param = \exp((Intercept) + \sum \beta X)\$\$ except for the shares
\\\rho\\ of the Delaporte and Hermite families (logit link) and the
Sichel \\\gamma\\ (identity link). The parameter modeled by each formula
in each family is listed in
[count_families](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md).

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

Underdispersion: when the dispersion parameter is negative, the PMF is
truncated at this point and renormalized, so its mean and variance
differ from the formulas above and \\\mu\\ is a location parameter
rather than the mean (e.g. GP1 with \\\mu = 0.75\\ and \\\phi = -0.5\\
has mean 0.803; if \\-\mu/\phi \< 1\\ all probability is at 0).
Predictions, Pearson residuals, simulations and marginal effects use the
exact moments of the renormalized PMF.

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

Underdispersion: when the dispersion parameter is negative, the PMF is
truncated at this point and renormalized, so its mean and variance
differ from the formulas above and \\\mu\\ is a location parameter
rather than the mean (e.g. GP1 with \\\mu = 0.75\\ and \\\phi = -0.5\\
has mean 0.803; if \\-\mu/\phi \< 1\\ all probability is at 0).
Predictions, Pearson residuals, simulations and marginal effects use the
exact moments of the renormalized PMF.

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
\frac{\Gamma(y + 1/\alpha)}{y! \\ \Gamma(1/\alpha)} \left( \frac{1}{1 +
\alpha\mu} \right)^{1/\alpha} \left( \frac{\alpha\mu}{1 + \alpha\mu}
\right)^y\$\$ \$\$LL\_{\text{NB2}} = \sum\_{i=1}^n \left\[ \ln
\Gamma(y_i + 1/\alpha) - \ln \Gamma(1/\alpha) - \ln y_i! -
\frac{1}{\alpha} \ln(1 + \alpha\mu_i) + y_i \ln \left(
\frac{\alpha\mu_i}{1 + \alpha\mu_i} \right) \right\]\$\$

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

Halton draws are used to integrate over the generalized-exponential
distribution.

Substituting \\x = \gamma u\\ shows that the scale \\\gamma\\ cancels
out of this PMF, so it cannot be estimated: it is fixed at 1 and only
the shape \\\alpha\\ is estimated (`dis_param_formula_2` does not
apply).

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

The PMF above is evaluated at this \\\lambda\\.

The variance function is defined as:
\$\$\sigma^2=\mu+\left(1-\frac{2}{(\theta+2)^2}\right)\mu^2\$\$

It should be noted that the p-value for the parameter `ln(theta)` in the
model summary is testing if the parameter `theta` is equal to a value
of 1. This has no practical meaning. The PL model does not reduce to the
Poisson model for any \\\theta\\ (its variance is at least \\\mu +
\mu^2/2\\), so a likelihood-ratio test against Poisson is not valid;
compare the models with AIC and BIC.

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

The Poisson-Lindley mean parameter is \\\mu=e^{X\beta}\\, and the mean
of the regression model (used by
[`predict()`](https://rdrr.io/r/stats/predict.html)) is:
\$\$E\[Y\]=m=e^{X\beta+\sigma^2/2}\$\$

The variance function is:
\$\$Var(Y)=m+\left(\left(2-\frac{2}{(\theta+2)^2}\right)
e^{\sigma^2}-1\right)m^2\$\$

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

The scale \\\sigma\\ cancels out of the PMF after this substitution, so
it cannot be estimated: it is fixed at 1 and only the shape \\\alpha\\
is estimated (`dis_param_formula_2` does not apply).

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

### Families added in version 0.1.2

All of these families are parameterized by their mean, \\E\[Y\] = \mu =
\exp(X\beta)\\, so the coefficients are effects on the expected count
for every family, and an offset multiplies the mean. Their distribution
functions are documented separately (see the links).

- **DW**
  ([`ddweib`](https://jwood-iastate.github.io/flexCountReg/reference/DiscreteWeibull.md)):
  discrete Weibull with shape \\\beta\\; \\\theta\\ is found from
  \\\mu\\ and \\\beta\\. It can be under- or overdispersed; \\\beta =
  1\\ is geometric.

- **DPO**
  ([`ddpo`](https://jwood-iastate.github.io/flexCountReg/reference/DoublePoisson.md)):
  double Poisson, exactly normalized, with precision \\\phi\\; its
  location \\\lambda\\ is found from \\\mu\\ and \\\phi\\. \\\phi = 1\\
  is Poisson, \\\phi \> 1\\ underdispersed, \\\phi \< 1\\ overdispersed.

- **DELAPORTE**
  ([`ddelaporte`](https://jwood-iastate.github.io/flexCountReg/reference/Delaporte.md)):
  \\Var(Y) = \mu + (\rho\mu)^2 / k\\. \\\rho\\ and \\k\\ are weakly
  identified when \\\rho\\ is near 0, where Poisson is a limit.

- **HERMITE**, **HERMITE3**, **HERMITE4**
  ([`dghermite`](https://jwood-iastate.github.io/flexCountReg/reference/GeneralizedHermite.md)):
  events occurring singly or in batches of \\m = 2, 3, 4\\; \\Var(Y) =
  \mu(1 + (m - 1)\rho)\\, so the variance-to-mean ratio is between 1 and
  \\m\\.

- **HYPERPOIS**
  ([`dhpois`](https://jwood-iastate.github.io/flexCountReg/reference/HyperPoisson.md)):
  hyper-Poisson with dispersion \\\gamma\\; its rate \\\lambda\\ is
  found from \\\mu\\ and \\\gamma\\. \\\gamma = 1\\ is Poisson, \\\gamma
  \> 1\\ overdispersed, \\\gamma \< 1\\ underdispersed.

- **GAMMACOUNT**
  ([`dgcount`](https://jwood-iastate.github.io/flexCountReg/reference/GammaCount.md)):
  the number of events in a unit window of a renewal process with gamma
  inter-arrival times of shape \\a\\; the rate is found from \\\mu\\ and
  \\a\\. \\a = 1\\ is Poisson, \\a \> 1\\ underdispersed, \\a \< 1\\
  overdispersed. An offset scales the mean (a regression convention); it
  does not lengthen the window of one renewal process, whose expected
  count is not proportional to the window length.

- **PA**
  ([`dpaeppli`](https://jwood-iastate.github.io/flexCountReg/reference/PolyaAeppli.md)):
  the Poisson-Tweedie distribution with power 1.5, \\Var(Y) = \mu +
  \phi\mu^{1.5}\\.

- **NBGE**
  ([`dnbge`](https://jwood-iastate.github.io/flexCountReg/reference/NBGE.md)):
  NB2 with a mean-one generalized exponential random effect, \\Var(Y) =
  \mu + ((1 + \alpha)M_2 - 1)\mu^2\\. Its two dispersion parameters have
  confounded effects on the variance; expect weak identification (large
  or missing standard errors) in small samples.

For DW, DPO, HYPERPOIS and GAMMACOUNT the natural parameter of every
observation is found by solving the mean equation numerically (Newton
iterations in compiled code), so these fits are slower than closed-form
families, especially in
[`countreg.rp`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.rp.md).
Their shape parameter is restricted to \\\[10^{-6}, 10^6\]\\ and the
mean to at most \\10^9\\; outside this numerical domain the likelihood
is treated as 0.

Likelihood-ratio tests against the Poisson model are reported for DPO,
HYPERPOIS and GAMMACOUNT (Poisson is an interior special case). Binomial
thinning (`underreport_formula`) is available for DELAPORTE and NBGE
(same distribution parameters), HERMITE (\\\rho\\ becomes \\p\rho\\) and
PA (\\\phi\\ becomes \\\sqrt{p}\phi\\); the other new families are not
closed under thinning.

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
#> initial  value 1066.772641 
#> iter   2 value 1066.770434
#> iter   3 value 1066.651116
#> iter   4 value 1066.591134
#> iter   5 value 1066.227309
#> iter   6 value 1066.021057
#> iter   7 value 1065.598383
#> iter   8 value 1065.169686
#> iter   9 value 1064.882770
#> iter  10 value 1064.878491
#> iter  11 value 1064.877220
#> iter  12 value 1064.876265
#> iter  13 value 1064.876106
#> iter  13 value 1064.876105
#> iter  13 value 1064.876105
#> final  value 1064.876105 
#> converged
summary(nb2)
#> Call:
#> countreg(formula = Total_crashes ~ lnaadt + lnlength + speed50 + 
#>     AADT10kplus, data = washington_roads, family = "NB2", verbose = TRUE, 
#>     dis_param_formula_1 = ~speed50, method = "BFGS")
#> 
#> Model: countreg | family: NB2 | optimizer: BFGS
#> Iterations: 41 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>                       Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)            -7.4015     0.4729  -15.65  < 2e-16 ***
#> lnaadt                  0.9116     0.0566   16.09  < 2e-16 ***
#> lnlength                0.8428     0.0700   12.04  < 2e-16 ***
#> speed50                -0.4699     0.1151   -4.08  4.5e-05 ***
#> AADT10kplus             0.7697     0.1299    5.92  3.1e-09 ***
#> ln(alpha):(Intercept)  -1.6191     0.3740   -4.33  1.5e-05 ***
#> ln(alpha):speed50       1.3058     0.5963    2.19    0.029 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1064.876 | AIC: 2143.752 | BIC: 2180.949 
#> N: 1501


# Estimate a Poisson-Lognormal model (integration with the poilog package)
if (requireNamespace("poilog", quietly = TRUE)) {
  pln <- countreg(Total_crashes ~ lnaadt + lnlength + speed50 + AADT10kplus,
                  data = washington_roads, family = "PLN", engine="poilog")
  summary(pln)
}
#> Warning: The Hessian is not negative definite at the estimates, so the standard errors of (Intercept), lnaadt, ln(sigma) are NA. The estimate may be on the boundary of the parameter space or not at a maximum.
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
#> lnlength       0.847      0.112    7.58  3.5e-14 ***
#> speed50       -0.459      0.106   -4.34  1.4e-05 ***
#> AADT10kplus    0.809      0.133    6.07  1.3e-09 ***
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
#> (Intercept)                 -6.6923     1.2044   -5.56  2.7e-08 ***
#> lnaadt                       0.9358     0.0562   16.66  < 2e-16 ***
#> lnlength                     0.8517     0.0698   12.21  < 2e-16 ***
#> speed50                      0.5700     1.0308    0.55     0.58    
#> AADT10kplus                  0.1580     0.8404    0.19     0.85    
#> ln(alpha)                   -1.5854     0.3837   -4.13  3.6e-05 ***
#> Underreporting:(Intercept)  -0.2370     1.1619   -0.20     0.84    
#> Underreporting:speed50      -0.8505     0.9645   -0.88     0.38    
#> Underreporting:AADT10kplus   0.7909     0.7482    1.06     0.29    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1065.466 | AIC: 2148.933 | BIC: 2196.758 
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
#> (Intercept)  -7.4230     0.4677  -15.87  < 2e-16 ***
#> lnaadt        0.9127     0.0565   16.16  < 2e-16 ***
#> lnlength      0.8400     0.0668   12.58  < 2e-16 ***
#> speed50      -0.4451     0.1031   -4.32  1.6e-05 ***
#> AADT10kplus   0.7872     0.1259    6.25  4.1e-10 ***
#> ln(nu)       -0.6199     0.1715   -3.62    3e-04 ***
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
