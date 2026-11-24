# Plot estimated smooth functions

Plots the smooth terms of a model with
[`smooth_estimates()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_estimates.md):
one panel per smooth, with the pointwise band shaded, the simultaneous
band dashed and the observed values as a rug. Two-variable smooths are
drawn as contour plots of the estimate.

## Usage

``` r
plot_smooths(
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
[`smooth_estimates()`](https://jwood-iastate.github.io/flexCountReg/reference/smooth_estimates.md),
invisibly.

## Examples

``` r
# \donttest{
data("washington_roads")
m <- countreg(Total_crashes ~ s(lnaadt) + offset(lnlength),
              data = washington_roads, family = "NB2")
plot_smooths(m)

# }
```
