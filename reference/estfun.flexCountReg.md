# Scores and bread for the sandwich package

Methods that make `flexCountReg` fits work with
[`sandwich::vcovCL()`](https://zeileis.codeberg.page/sandwich/reference/vcovCL.html),
[`sandwich::vcovHC()`](https://zeileis.codeberg.page/sandwich/reference/vcovHC.html),
[`sandwich::sandwich()`](https://zeileis.codeberg.page/sandwich/reference/sandwich.html)
and, through them,
[`lmtest::coeftest()`](https://rdrr.io/pkg/lmtest/man/coeftest.html) and
`modelsummary(vcov = ~ cluster)`.

## Usage

``` r
# S3 method for class 'flexCountReg'
estfun(x, ...)

# S3 method for class 'flexCountReg'
bread(x, ...)
```

## Arguments

- x:

  a fitted `flexCountReg` model.

- ...:

  not used.

## Value

`estfun()`: a matrix with one row per unit and one column per parameter.
`bread()`: a square matrix.

## Details

`estfun()` returns the scores of the model's independent units: one row
per row of the data used in estimation, or per panel for panel models,
multiplied by the weights. `bread()` is \\n (-H)^{-1}\\ with \\n\\ the
number of units, the convention of the sandwich package, so
`sandwich::sandwich(m)` is \\(-H)^{-1} \sum_i s_i s_i' (-H)^{-1}\\.

The sandwich package treats each row of `estfun()` as an independent
unit. That is right without weights and for sampling weights; with
frequency weights, use `vcov(m, type = "robust")` or `type = "cluster"`,
which count a row with weight \\w\\ as \\w\\ observations. For panel
models, a cluster variable for `vcovCL()` must have one value per panel
(in the order of the panel ids);
`vcov(m, type = "cluster", cluster = ~ v)` takes the variable from the
data.

## Examples

``` r
# \donttest{
data("washington_roads")
m <- countreg(Total_crashes ~ lnaadt + lnlength, data = washington_roads,
              family = "NB2", method = "BFGS")
sandwich::vcovCL(m, cluster = washington_roads$ID, type = "HC0")
#>              (Intercept)        lnaadt      lnlength    ln(alpha)
#> (Intercept)  0.343855392 -3.923457e-02  7.031254e-03  0.022182837
#> lnaadt      -0.039234573  4.610843e-03  3.194782e-05 -0.002843661
#> lnlength     0.007031254  3.194782e-05  7.591455e-03 -0.001220079
#> ln(alpha)    0.022182837 -2.843661e-03 -1.220079e-03  0.058516309
# }
```
