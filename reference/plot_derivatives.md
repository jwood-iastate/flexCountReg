# Plot derivatives of smooth functions

Plot derivatives of smooth functions

## Usage

``` r
plot_derivatives(
  object,
  term = NULL,
  n = 100,
  level = 0.95,
  nsim = 10000,
  seed = 1,
  ...
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

- ...:

  further arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

The data of
[`smooth_derivatives()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_derivatives.md),
invisibly.
