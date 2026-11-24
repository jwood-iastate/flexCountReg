# Derivatives of smooth functions

The first derivative \\df/dx\\ of each one-variable smooth, with
pointwise and simultaneous confidence bands: where the band excludes 0,
the smooth is increasing (or decreasing).

## Usage

``` r
smooth_derivatives(
  object,
  term = NULL,
  n = 100,
  level = 0.95,
  nsim = 10000,
  seed = NULL,
  h = NULL
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

- h:

  step for the central differences of the basis (default `1e-4` times
  the range of the variable).

## Value

A data frame like that of
[`smooth_estimates()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_estimates.md),
with the derivative in `fit`.

## Details

The derivative is exact for the linear predictor up to the
central-difference error of the basis: with \\D = (X(x + h) - X(x - h))
/ 2h\\, \\df/dx = D b\\ and its standard error is \\\sqrt{D V D'}\\.
Simultaneous bands are computed as in
[`smooth_estimates()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_estimates.md).

## See also

[`plot_derivatives()`](https://jwood-iastate.github.io/flexCountReg/reference/plot_derivatives.md)
