# Function for estimating a Random Effects Poisson-Lindley regression model

Function for estimating a Random Effects Poisson-Lindley regression
model

## Usage

``` r
poisLind.re(
  formula,
  group_var,
  data,
  method = "NM",
  max.iters = 1000,
  verbose = FALSE,
  stderr = "normal",
  bootstraps = NULL,
  offset = NULL,
  weights = NULL,
  keep_data = TRUE,
  print.level = NULL
)
```

## Arguments

- formula:

  an R formula.

- group_var:

  the grouping variable(s) indicating random effects (e.g., individual
  ID). Required. Several columns are combined into one panel identifier.

- data:

  a dataframe that has all of the variables in the `formula`.

- method:

  a method to use for optimization in the maximum likelihood estimation.
  For options, see
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html). Note that
  "BHHH" is not available for this function due to the implementation
  for the random effects.

- max.iters:

  the maximum number of iterations to allow the optimization method to
  perform.

- verbose:

  logical; if `TRUE`, print the progress of the optimization.

- stderr:

  the type of standard errors: `"normal"` (Hessian-based, the default),
  `"robust"` (sandwich; clustered by panel because the likelihood has
  one contribution per panel), or `"boot"` (bootstrap that resamples
  whole panels). Giving `bootstraps` without `stderr` also requests
  bootstrapped standard errors.

- bootstraps:

  Optional integer specifying the number of bootstrap samples to be used
  for estimating standard errors. Whole panels (groups) are resampled,
  so the standard errors account for the correlation within panels. If
  not specified, no bootstrapping is performed.

- offset:

  an optional character vector of offset column names (summed when
  several are given; on the log scale). An offset can instead be given
  with [`offset()`](https://rdrr.io/r/stats/offset.html) in `formula`;
  use only one of the two.

- weights:

  the name of a column of \`data\`, or a numeric vector with one value
  per row of \`data\`, of frequency weights of the panels: a panel with
  weight w counts as w independent panels in the log-likelihood, the
  robust covariance, BIC and the bootstrap. Weights must be constant
  within each panel. Bootstrap standard errors need whole-number
  weights.

- keep_data:

  logical; if `FALSE`, the training data are not stored in the fitted
  object (saving memory).
  [`predict()`](https://rdrr.io/r/stats/predict.html) then needs
  `newdata`, and
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html)/[`residuals()`](https://rdrr.io/r/stats/residuals.html)
  are unavailable.

- print.level:

  deprecated; use `verbose`.

## Value

An object of class `flexCountReg` which is a list with the following
components:

- model: the fitted model object.

- data: the data frame used to fit the model.

- call: the matched call.

- formula: the formula used to fit the model.

## Details

`poislind_re()` is the same function as `poisLind.re()`, under a name
that tools do not mistake for an S3 method.

The function `poisLindRE` is similar to the `poisLind` function, but it
includes an additional argument `group_var` that specifies the grouping
variable for the random effects. The function estimates a Random Effects
Poisson-Lindley regression model using maximum likelihood. It is similar
to `poisLind`, but includes additional terms to account for the random
effects.

The Random Effects Poisson-Lindley model is useful for panel data and
assumes that the random effects follow a gamma distribution. The PDF is
\$\$ f(y\_{it}\|\mu\_{it},\theta)=\frac{\theta^2}{\theta+1}
\prod\_{t=1}^{n_i}\frac{\left(\mu\_{it}\frac{\theta(\theta+1)}
{\theta+2}\right)^{y\_{it}}}{y\_{it}!} \cdot \frac{
\left(\sum\_{t=1}^{n_i}y\_{it}\right)!
\left(\sum\_{t=1}^{n_i}\mu\_{it}\frac{\theta(\theta+1)}{\theta+2} +
\theta + \sum\_{t=1}^{n_i}y\_{it} + 1\right) }{
\left(\sum\_{t=1}^{n_i}\mu\_{it}\frac{\theta(\theta+1)}{\theta+2} +
\theta\right)^{\sum\_{t=1}^{n_i}y\_{it}+2} } \$\$

The log-likelihood function is: \$\$ LL = 2\log(\theta) -
\log(\theta+1) + \sum\_{t=1}^{n_i} y\_{it}\log(\mu\_{it}) +
\sum\_{t=1}^{n_i} y\_{it}\log\\\left( \frac{\theta(\theta+1)}{\theta+2}
\right) - \sum\_{t=1}^{n_i}\log(y\_{it}!) + \log\\\left(
\left(\sum\_{t=1}^{n_i}y\_{it}\right)! \right) + \log\\\left(
\sum\_{t=1}^{n_i}\mu\_{it}\frac{\theta(\theta+1)}{\theta+2} + \theta +
\sum\_{t=1}^{n_i}y\_{it} + 1 \right) - \left(\sum\_{t=1}^{n_i}y\_{it} +
2\right) \log\\\left(
\sum\_{t=1}^{n_i}\mu\_{it}\frac{\theta(\theta+1)}{\theta+2} + \theta
\right) \$\$

The mean and variance are: \$\$\mu\_{it}=\exp(X\_{it} \beta)\$\$ \$\$
V(\mu\_{it})=\mu\_{it}+ \left(1-\frac{2}{(\theta+2)^2}\right)\mu\_{it}^2
\$\$

## Examples

``` r
# \donttest{
data("washington_roads")
washington_roads$AADTover10k <-
  ifelse(washington_roads$AADT > 10000, 1, 0)

poislind.mod <- poisLind.re(
  Animal ~ lnaadt + lnlength + speed50 +
    ShouldWidth04 + AADTover10k,
  data      = washington_roads,
  group_var = "ID",
  method    = "NM",
  max.iters = 1000
)
summary(poislind.mod)
#> Call:
#> poisLind.re(formula = Animal ~ lnaadt + lnlength + speed50 + 
#>     ShouldWidth04 + AADTover10k, group_var = "ID", data = washington_roads, 
#>     method = "NM", max.iters = 1000)
#> 
#> Model: poislind_re | optimizer: NM
#> Iterations: 828 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)     -9.623      1.426   -6.75  1.5e-11 ***
#> lnaadt           1.036      0.167    6.19  5.9e-10 ***
#> lnlength         1.456      0.236    6.18  6.4e-10 ***
#> speed50         -0.855      0.335   -2.55    0.011 *  
#> ShouldWidth04   -0.434      0.291   -1.49    0.136    
#> AADTover10k     -0.756      0.532   -1.42    0.155    
#> ln(theta)        6.068         NA      NA       NA    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -263.098 | AIC: 540.195 | BIC: 569.795 
#> N: 1501 | panels: 507
# }
```
