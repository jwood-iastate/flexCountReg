# Estimate a Random Effects Poisson Panel Model

Estimate a random effects Poisson panel model for clustered count data.
The model is also known as a negative multinomial model. A more flexible
version of this model can be estimated using the
[`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)
function.

## Usage

``` r
repois(
  formula,
  group_var,
  data,
  method = "NM",
  max.iters = 1000,
  verbose = FALSE,
  stderr = "normal",
  bootstraps = NULL,
  offset = NULL,
  keep_data = TRUE,
  print.level = NULL
)
```

## Arguments

- formula:

  an R formula.

- group_var:

  the grouping variable(s) for the random effects (e.g., individual ID
  or other panel ID variables).

- data:

  a dataframe that has all of the variables in the `formula`.

- method:

  a method to use for optimization in the maximum likelihood estimation.
  For options, see
  [`maxLik`](https://rdrr.io/pkg/maxLik/man/maxLik.html).

- max.iters:

  the maximum number of iterations to allow the optimization method to
  perform.

- verbose:

  logical; if `TRUE`, print the progress of the optimization.

- stderr:

  the type of standard errors: `"normal"` (Hessian-based, the default),
  `"robust"` (sandwich, clustered by group), or `"boot"` (bootstrap that
  resamples whole groups). Giving `bootstraps` without `stderr` also
  requests bootstrapped standard errors.

- bootstraps:

  optional integer specifying the number of bootstrap samples to be used
  for estimating standard errors. If not specified, no bootstrapping is
  performed.

- offset:

  an optional offset term provided as a string.

- keep_data:

  logical; if `FALSE`, the training data are not stored in the fitted
  object.

- print.level:

  deprecated; use `verbose`.

## Value

An object of class `countreg` which is a list with the following
components:

- model: the fitted model object.

- data: the data frame used to fit the model.

- call: the matched call.

- formula: the formula used to fit the model.

## Details

This function estimates a random effects Poisson panel model. The
likelihood is integrated analytically over a gamma-distributed cluster
effect with mean 1 and variance `alpha`. The model estimates the
regression coefficients and one heterogeneity parameter `ln(alpha)`.

The PMF of this distribution is: \$\$P(y\_{it} \mid \mu\_{it}) =
\frac{\Gamma\left(\sum_t y\_{it} + \frac{1}{\alpha}\right)}
{\Gamma\left(\frac{1}{\alpha}\right)\prod_t y\_{it}!}
\left(\frac{1}{\alpha\sum_t \mu\_{it} + 1}\right)^{\frac{1}{\alpha}}
\prod_t \left(\frac{\alpha\mu\_{it}}{\alpha\sum_t \mu\_{it} + 1}
\right)^{y\_{it}}\$\$

## Examples

``` r
# \donttest{
## RE Poisson Panel Model
data("washington_roads")
washington_roads$AADTover10k <-
  ifelse(washington_roads$AADT > 10000, 1, 0)

repois.mod <- repois(
  Animal ~ lnaadt + speed50 + ShouldWidth04 + AADTover10k,
  data = washington_roads,
  offset = "lnlength",
  group_var = "ID",
  method = "BFGS",
  max.iters = 1000
)

summary(repois.mod)
#> Call:
#> repois(formula = Animal ~ lnaadt + speed50 + ShouldWidth04 + 
#>     AADTover10k, group_var = "ID", data = washington_roads, method = "BFGS", 
#>     max.iters = 1000, offset = "lnlength")
#> 
#> Model: repois | optimizer: BFGS
#> Iterations: 7 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>               Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)     -9.818      1.270   -7.73  1.1e-14 ***
#> lnaadt           1.026      0.147    7.00  2.5e-12 ***
#> speed50         -0.919      0.290   -3.16   0.0016 ** 
#> ShouldWidth04   -0.486      0.243   -2.00   0.0456 *  
#> AADTover10k     -0.870      0.459   -1.89   0.0582 .  
#> lnalpha        -13.816         NA      NA       NA    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> Offset: lnlength (coefficient fixed at 1)
#> 
#> Log-likelihood: -271.719 | AIC: 555.438 | BIC: 580.809 
#> N: 1501 | panels: 507
# }
```
