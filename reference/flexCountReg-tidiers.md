# Tidy, glance and augment flexCountReg models

Methods for the
[`tidy()`](https://generics.r-lib.org/reference/tidy.html),
[`glance()`](https://generics.r-lib.org/reference/glance.html) and
[`augment()`](https://generics.r-lib.org/reference/augment.html)
generics from the generics package (also used by broom), so flexCountReg
models work with broom-style workflows, tables such as modelsummary, and
tidymodels.

## Usage

``` r
# S3 method for class 'flexCountReg'
tidy(x, conf.int = FALSE, conf.level = 0.95, exponentiate = FALSE, ...)

# S3 method for class 'flexCountReg'
glance(x, ...)

# S3 method for class 'flexCountReg'
augment(
  x,
  data = NULL,
  new_data = NULL,
  newdata = NULL,
  type.residuals = c("response", "pearson"),
  ...
)
```

## Arguments

- x:

  a model estimated with this package.

- conf.int:

  logical; include confidence intervals?

- conf.level:

  the confidence level for `conf.int = TRUE`.

- exponentiate:

  logical; exponentiate the mean-model coefficients (giving incidence
  rate ratios)? Other components are not transformed.

- ...:

  further arguments passed to
  [`predict()`](https://rdrr.io/r/stats/predict.html) (for example
  `method` for
  [`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
  models).

- data:

  the data to augment. Defaults to the rows of the training data used in
  estimation.

- new_data, newdata:

  new data to augment instead of `data`; the outcome is optional.

- type.residuals:

  the residual type for `.resid`: `"response"` or `"pearson"` (see
  [`residuals.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/flexCountReg-methods.md)).

## Value

[`tidy()`](https://generics.r-lib.org/reference/tidy.html): a tibble
with columns `term`, `component`, `estimate`, `std.error`, `statistic`,
`p.value` (and `conf.low`, `conf.high`).
[`glance()`](https://generics.r-lib.org/reference/glance.html): a
one-row tibble of model statistics.
[`augment()`](https://generics.r-lib.org/reference/augment.html): the
data as a tibble with `.fitted` (and `.resid` when the outcome is
available).

## Details

[`tidy()`](https://generics.r-lib.org/reference/tidy.html) has one row
per parameter, with a `component` column because these models have
several equations:

- mean:

  coefficients of the mean (count) model;

- distribution:

  parameters of the count distribution (for example `ln(alpha)`),
  including their own regression coefficients;

- underreporting:

  coefficients of the reporting model;

- zero:

  coefficients of the zero model of zero-inflated and hurdle models;

- random_mean, random_scale, random_cholesky:

  means, scales and Cholesky elements of random parameters;

- heterogeneity:

  heterogeneity-in-means and -variances coefficients.

Standard errors are those chosen when the model was fit (see
[`vcov()`](https://rdrr.io/r/stats/vcov.html)), and `p.value` uses the
normal distribution.

## Examples

``` r
# \donttest{
data("washington_roads")
m <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
              family = "NB2", method = "BFGS")
tidy(m, conf.int = TRUE)
#> # A tibble: 4 × 8
#>   term        component estimate std.error statistic  p.value conf.low conf.high
#>   <chr>       <chr>        <dbl>     <dbl>     <dbl>    <dbl>    <dbl>     <dbl>
#> 1 (Intercept) mean        -9.21     0.445     -20.7  2.06e-95  -10.1      -8.34 
#> 2 lnaadt      mean         1.12     0.0529     21.1  1.01e-98    1.01      1.22 
#> 3 lnlength    mean         0.744    0.0696     10.7  1.13e-26    0.608     0.881
#> 4 ln(alpha)   distribu…   -0.916    0.234      -3.92 8.81e- 5   -1.37     -0.458
tidy(m, exponentiate = TRUE)
#> # A tibble: 4 × 6
#>   term        component      estimate std.error statistic  p.value
#>   <chr>       <chr>             <dbl>     <dbl>     <dbl>    <dbl>
#> 1 (Intercept) mean          0.0000998 0.0000444    -20.7  2.06e-95
#> 2 lnaadt      mean          3.05      0.162         21.1  1.01e-98
#> 3 lnlength    mean          2.10      0.146         10.7  1.13e-26
#> 4 ln(alpha)   distribution -0.916     0.234         -3.92 8.81e- 5
glance(m)
#> # A tibble: 1 × 10
#>   model    family logLik   AIC   BIC    df  nobs n_groups converged vcov_type
#>   <chr>    <chr>   <dbl> <dbl> <dbl> <int> <int>    <int> <lgl>     <chr>    
#> 1 countreg NB2    -1098. 2204. 2225.     4  1501       NA TRUE      hessian  
head(augment(m))
#> # A tibble: 6 × 15
#>   ID     Year  AADT Length Total_crashes lnaadt lnlength speed50 ShouldWidth04
#>   <fct> <int> <dbl>  <dbl>         <int>  <dbl>    <dbl>   <int>         <int>
#> 1 1      2016  7819  0.430             0   8.96   -0.844       1             0
#> 2 2      2016  7819  0.380             2   8.96   -0.968       1             0
#> 3 3      2016  7819  0.630             2   8.96   -0.462       1             0
#> 4 4      2016  7819  0.140             0   8.96   -1.97        1             0
#> 5 5      2016  7819  0.450             0   8.96   -0.799       1             0
#> 6 6      2016  7819  0.210             1   8.96   -1.56        1             0
#> # ℹ 6 more variables: Fatal_crashes <dbl>, Injury_crashes <dbl>, Animal <dbl>,
#> #   Rollover <dbl>, .fitted <dbl>, .resid <dbl>
# }
```
