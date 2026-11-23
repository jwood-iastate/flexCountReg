# Define a count distribution for the regression functions

Creates a count family from its log probability mass function, for use
as `family` in
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md),
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
and
[`countreg_mix()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mix.md).
Only `logpmf` is required: the mean, variance, CDF and random draws are
computed from it when they are not given.

## Usage

``` r
new_count_family(
  name,
  logpmf,
  params = list(),
  mean_is_mu = TRUE,
  mean = NULL,
  variance = NULL,
  cdf = NULL,
  rand = NULL,
  valid = NULL
)
```

## Arguments

- name:

  a short name for the family (used in printed output).

- logpmf:

  `function(y, mu, par)` returning \\\log P(Y = y)\\ for vectors `y` and
  `mu` and a named list `par` of parameter vectors of the same length.
  Invalid parameter values should give `-Inf` or `NaN` (errors are also
  caught).

- params:

  a named list of up to two
  [`count_param()`](https://jwood-iastate.github.io/flexCountReg/reference/count_param.md)
  objects. The first is modelled by `dis_param_formula_1`, the second by
  `dis_param_formula_2`. The coefficients are labelled `ln(name)`,
  `logit(name)` or `name`, by link.

- mean_is_mu:

  `TRUE` (default) if \\\mu = \exp(X\beta)\\ is the mean of the
  distribution. Otherwise give `mean`, or the mean is computed by
  summing the PMF.

- mean, variance:

  optional `function(mu, par)`: the mean and variance. Without them they
  are computed by summing the PMF.

- cdf:

  optional `function(q, mu, par)`: \\P(Y \le q)\\. If it has arguments
  `lower.tail` and `log.p` they are used, which is more accurate in the
  tails. Without it the PMF is summed.

- rand:

  optional `function(n, mu, par)`: random draws. Without it draws are
  made by inverting the CDF.

- valid:

  optional `function(mu, par)`: `TRUE` where the parameters define a
  distribution. Without it every finite value is accepted and `logpmf`
  must signal invalid values.

## Value

An object of class `"fcr_count_family"`.

## Details

The family supports everything that needs only probabilities: fixed and
random-parameter models, zero-inflated and hurdle models, finite
mixtures, predictions of probabilities, quantile residuals, rootograms
and simulation. Underreporting (binomial thinning) is not available,
because a user family is not known to be closed under thinning.

The fit stores the family, so a saved model works in a new R session.
Two families with the same `name` but different definitions are kept
apart.

Use
[`check_count_family()`](https://jwood-iastate.github.io/flexCountReg/reference/check_count_family.md)
to test a new family before fitting.

## Examples

``` r
# Negative binomial with size theta (variance mu + mu^2 / theta)
nb_theta <- new_count_family(
  name = "NB-theta",
  params = list(theta = count_param("log")),
  logpmf = function(y, mu, par)
    dnbinom(y, size = par$theta, mu = mu, log = TRUE),
  variance = function(mu, par) mu + mu^2 / par$theta)
# \donttest{
data("washington_roads")
m <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
              family = nb_theta, method = "BFGS")
summary(m)
#> Call:
#> countreg(formula = Total_crashes ~ lnaadt + lnlength, data = washington_roads, 
#>     family = nb_theta, method = "BFGS")
#> 
#> Model: countreg | family: NB-theta | optimizer: BFGS
#> Iterations: 33 
#> Convergence: successful convergence 
#> 
#> Coefficients (standard errors: Hessian-based):
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)  -9.2125     0.4442  -20.74  < 2e-16 ***
#> lnaadt        1.1159     0.0529   21.10  < 2e-16 ***
#> lnlength      0.7441     0.0696   10.69  < 2e-16 ***
#> ln(theta)     0.9162     0.2335    3.92  8.7e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Log-likelihood: -1097.96 | AIC: 2203.92 | BIC: 2225.176 
#> N: 1501
# }
```
