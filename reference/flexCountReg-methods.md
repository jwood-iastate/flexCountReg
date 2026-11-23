# Standard methods for flexCountReg models

Methods that make flexCountReg models work like other fitted models in
R, including with [`AIC()`](https://rdrr.io/r/stats/AIC.html),
[`BIC()`](https://rdrr.io/r/stats/AIC.html),
[`confint()`](https://rdrr.io/r/stats/confint.html),
[`update()`](https://rdrr.io/r/stats/update.html) and packages that
build on these generics.

## Usage

``` r
# S3 method for class 'flexCountReg'
coef(object, ...)

# S3 method for class 'flexCountReg'
nobs(object, ...)

# S3 method for class 'flexCountReg'
logLik(object, ...)

# S3 method for class 'flexCountReg'
formula(x, ...)

# S3 method for class 'flexCountReg'
terms(x, ...)

# S3 method for class 'flexCountReg'
confint(object, parm, level = 0.95, ...)

# S3 method for class 'flexCountReg'
fitted(object, ...)

# S3 method for class 'flexCountReg'
residuals(object, type = c("response", "pearson", "quantile"), ...)

# S3 method for class 'flexCountReg'
simulate(object, nsim = 1, seed = NULL, ...)

# S3 method for class 'flexCountReg'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- object, x:

  a model estimated with this package.

- ...:

  further arguments passed to or from other methods. For
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`residuals()`](https://rdrr.io/r/stats/residuals.html) of
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  models, this can include `method` (see
  [`predict.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md)).

- parm:

  a specification of which parameters are to be given confidence
  intervals, either a vector of numbers or a vector of names. If
  missing, all parameters are considered.

- level:

  the confidence level required.

- type:

  the type of residuals: `"response"` (observed minus fitted),
  `"pearson"` (divided by the standard deviation implied by the fitted
  distribution) or `"quantile"` (randomized quantile residuals of Dunn
  and Smyth, 1996, which are standard normal when the model is correct).
  Pearson and quantile residuals are available for
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
  models.

- nsim:

  number of response vectors to simulate.

- seed:

  an optional random seed, as in
  [`simulate`](https://rdrr.io/r/stats/simulate.html).

- digits:

  the number of significant digits to print.

## Value

[`coef()`](https://rdrr.io/r/stats/coef.html): named numeric vector.
[`vcov()`](https://rdrr.io/r/stats/vcov.html): matrix.
[`logLik()`](https://rdrr.io/r/stats/logLik.html): an object of class
`logLik`. [`nobs()`](https://rdrr.io/r/stats/nobs.html): integer.
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html): numeric
vectors. [`confint()`](https://rdrr.io/r/stats/confint.html): matrix
with lower and upper limits.
[`formula()`](https://rdrr.io/r/stats/formula.html): the main model
formula. [`terms()`](https://rdrr.io/r/stats/terms.html): its terms
object. [`simulate()`](https://rdrr.io/r/stats/simulate.html): a data
frame with `nsim` columns.

## Details

[`coef()`](https://rdrr.io/r/stats/coef.html) returns all estimated
parameters, including distribution, underreporting, random-parameter and
heterogeneity parameters, named as in
[`summary()`](https://rdrr.io/r/base/summary.html).
[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns the covariance
matrix chosen at fit time: Hessian-based, robust (sandwich),
cluster-robust or bootstrapped (see the `stderr` argument of each
estimator); see
[`vcov.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/vcov.flexCountReg.md)
to recompute it after the fit.

[`nobs()`](https://rdrr.io/r/stats/nobs.html) is the number of rows used
in estimation. For panel models
([`renb()`](https://jwood-iastate.github.io/flexCountReg/reference/renb.md),
[`poisLind.re()`](https://jwood-iastate.github.io/flexCountReg/reference/poisLindRE.md),
[`repois()`](https://jwood-iastate.github.io/flexCountReg/reference/repois.md),
[`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)
and
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
with `panel_id`), the `nobs` attribute of
[`logLik()`](https://rdrr.io/r/stats/logLik.html), which
[`BIC()`](https://rdrr.io/r/stats/AIC.html) uses, is the number of
panels, because the panels are the independent units. With frequency
weights (`weights`), it is the sum of the weights of the rows (or
panels), the size of the expanded data.

[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) are for the rows
used in estimation (rows with missing values are not included).
[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws new
responses from the fitted distribution for the same rows; it is
available for
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
models.

## Examples

``` r
# \donttest{
data("washington_roads")
m <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
              family = "NB2", method = "BFGS")
coef(m)
#> (Intercept)      lnaadt    lnlength   ln(alpha) 
#>  -9.2124976   1.1159470   0.7440805  -0.9162359 
confint(m)
#>                   2.5 %     97.5 %
#> (Intercept) -10.0837227 -8.3412726
#> lnaadt        1.0122318  1.2196622
#> lnlength      0.6076599  0.8805011
#> ln(alpha)    -1.3742068 -0.4582650
AIC(m)
#> [1] 2203.92
BIC(m)
#> [1] 2225.176
head(residuals(m, type = "pearson"))
#> [1] -0.8946318  0.7475095  0.2732480 -0.6512653 -0.9049217  0.3294260
sims <- simulate(m, nsim = 2, seed = 1)
# }
```
