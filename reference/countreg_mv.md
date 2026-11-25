# Multivariate count regression

Fits several count outcomes observed on the same units, such as crashes
by severity or by type at each site, each with its own count family,
mean formula and distribution-parameter formulas.

## Usage

``` r
countreg_mv(
  formula,
  data,
  family = "NB2",
  dis_param_formula_1 = NULL,
  dis_param_formula_2 = NULL,
  offset = NULL,
  weights = NULL,
  group_var = NULL,
  stderr = "normal",
  method = "BFGS",
  max.iters = 1000,
  start.vals = NULL,
  ndraws = 1500,
  engine = "quadrature",
  verbose = FALSE,
  .setup_only = FALSE
)
```

## Arguments

- formula:

  either a formula with the outcomes on the left as
  `cbind(y_1, ..., y_J)` and one right-hand side for every equation, or
  a list of two-sided formulas, one per outcome (named by the list
  names, or by their responses). An
  [`offset()`](https://rdrr.io/r/stats/offset.html) term (an exposure)
  belongs in the formulas.

- data:

  a data frame.

- family:

  the count family of each outcome (see
  [count_families](https://jwood-iastate.github.io/flexCountReg/reference/count_families.md)):
  one code for every outcome, or one per outcome (named by outcome or in
  their order). Families made with
  [`new_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/new_count_family.md)
  are accepted.

- dis_param_formula_1, dis_param_formula_2:

  formulas for the outcomes' distribution parameters, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md):
  `NULL` (constant), one one-sided formula for every outcome that has
  the parameter, or a list named by outcome (outcomes not in the list
  keep a constant parameter).

- offset:

  optional names of offset columns (log scale) used by every equation,
  instead of [`offset()`](https://rdrr.io/r/stats/offset.html) in the
  formulas.

- weights:

  optional frequency weights: a column name or a numeric vector.

- group_var:

  optional column naming the units for cluster-robust standard errors
  (such as a site observed over several years).

- stderr:

  `"normal"` (from the Hessian, the default), `"robust"` (sandwich) or
  `"cluster"` (sandwich by `group_var`).

- method, max.iters, start.vals, verbose:

  optimizer settings, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md);
  the default method is `"BFGS"`.

- ndraws, engine:

  integration settings of families without a closed form, as in
  [`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).

- .setup_only:

  internal use: return the model setup without fitting.

## Value

An object of class `"countreg_mv"`, with methods
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html).

## Details

**Independent equations.** The outcomes are independent given the
covariates, so the log-likelihood is the sum of the equations'
log-likelihoods and every equation's estimates equal those of a separate
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
fit. Fitting them together gives one object with joint predictions,
simulations and tests across equations (for example, a Wald test that a
coefficient is the same in two equations). The rows used are those
complete for every equation.

Different outcomes often need different families: frequent property
damage crashes may be close to NB2, while rare fatal crashes are more
skewed and zero-heavy. Each equation can use any family of the package.

## See also

[`predict.countreg_mv()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.countreg_mv.md),
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
[`countreg_alloc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_alloc.md)
for a total and its split across categories,
[`negmultinom()`](https://jwood-iastate.github.io/flexCountReg/reference/negmultinom.md)

## Examples

``` r
# \donttest{
data("washington_roads")
w <- washington_roads
w$pdo <- w$Total_crashes - w$Injury_crashes
m <- countreg_mv(cbind(pdo, Injury_crashes) ~ lnaadt + speed50 +
                   offset(lnlength), data = w,
                 family = c(pdo = "NB2", Injury_crashes = "PIG"))
summary(m)
#> Multivariate count model
#> 
#> Rows:        1501 
#> Std. errors: normal 
#> Equations:   independent given the covariates
#> 
#> pdo (NB2):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept) -9.42403    0.47838 -19.700  < 2e-16 ***
#> lnaadt       1.17438    0.05483  21.420  < 2e-16 ***
#> speed50     -0.51500    0.11414  -4.512 6.42e-06 ***
#> ln(alpha)   -0.86569    0.23122  -3.744 0.000181 ***
#> 
#> Injury_crashes (PIG):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -7.4115     1.1903  -6.227 4.77e-10 ***
#> lnaadt        0.6661     0.1418   4.698 2.63e-06 ***
#> speed50      -1.2201     0.4146  -2.943  0.00325 ** 
#> ln(eta)       0.7900     1.9900   0.397  0.69137    
#> 
#> Log-likelihood: -1247.0544  AIC: 2510.1088  BIC: 2552.6199 
head(predict(m))
#>            pdo Injury_crashes
#> [1,] 0.7745464    0.030057397
#> [2,] 0.6844829    0.026562351
#> [3,] 1.1348006    0.044037581
#> [4,] 0.2521779    0.009786129
#> [5,] 0.8105719    0.031455415
#> [6,] 0.3782669    0.014679194
# }
```
