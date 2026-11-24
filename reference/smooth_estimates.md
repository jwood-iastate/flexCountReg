# Estimated smooth functions

Evaluates the smooth terms (`s()`, `te()`, ...) of a model fitted by
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
or
[`countreg_rp()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_rp.md)
on a grid over the observed range of their variables, with pointwise and
simultaneous confidence bands.

## Usage

``` r
smooth_estimates(
  object,
  term = NULL,
  n = 100,
  level = 0.95,
  nsim = 10000,
  seed = NULL
)
```

## Arguments

- object:

  a fitted model with smooth terms.

- term:

  optional term labels (as in
  [`summary()`](https://rdrr.io/r/base/summary.html), e.g.
  `"s(lnaadt)"`); by default every smooth.

- n:

  number of grid points for one-variable smooths (the grid of a
  two-variable smooth has `round(sqrt(n) * 3)` points per variable).

- level:

  confidence level of the bands.

- nsim:

  number of draws for the simultaneous bands.

- seed:

  optional seed for those draws (the caller's random-number stream is
  left unchanged).

## Value

A data frame with columns `parameter`, `term`, the smooth's variable(s),
`fit`, `se`, `lwr`, `upr` (pointwise) and `slwr`, `supr` (simultaneous).

## Details

Each smooth is shown on the scale of its linear predictor (log for the
mean, the parameter's link for a distribution parameter), centered as
fitted (the constant is in the intercept). Pointwise bands are \\\hat f
\pm z\\SE\\, with \\SE^2 = x' V x\\ from the fit's covariance matrix
(Bayesian by default). Simultaneous bands use the 95th (or `level`)
percentile of \\\max \|f^\* - \hat f\| / SE\\ over the grid, for
coefficient draws \\f^\*\\ from \\N(\hat b, V)\\ (Ruppert, Wand and
Carroll 2003), so they cover the whole curve with the stated
probability.

## References

Ruppert, D., Wand, M.P. and Carroll, R.J. (2003). *Semiparametric
Regression*. Cambridge University Press.

## See also

[`plot_smooths()`](https://jwood-iastate.github.io/flexCountReg/reference/plot_smooths.md)

## Examples

``` r
# \donttest{
data("washington_roads")
m <- countreg(Total_crashes ~ s(lnaadt) + offset(lnlength),
              data = washington_roads, family = "NB2")
head(smooth_estimates(m))
#>   parameter      term   lnaadt       fit        se       lwr        upr
#> 1      mean s(lnaadt) 5.796058 -2.108904 0.9008161 -3.874471 -0.3433364
#> 2      mean s(lnaadt) 5.837581 -2.003193 0.8112092 -3.593133 -0.4132518
#> 3      mean s(lnaadt) 5.879105 -1.897555 0.7235963 -3.315778 -0.4793326
#> 4      mean s(lnaadt) 5.920628 -1.792169 0.6391599 -3.044899 -0.5394385
#> 5      mean s(lnaadt) 5.962152 -1.687271 0.5592993 -2.783477 -0.5910644
#> 6      mean s(lnaadt) 6.003675 -1.583175 0.4855303 -2.534797 -0.6315527
#>        slwr         supr
#> 1 -4.826182  0.608375272
#> 2 -4.450175  0.443790163
#> 3 -4.080257  0.285146400
#> 4 -3.720171  0.135833542
#> 5 -3.374377 -0.000165082
#> 6 -3.047759 -0.118590152
# }
```
