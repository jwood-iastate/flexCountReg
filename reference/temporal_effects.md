# Period effects of time-varying terms

Computes the effect of each
[`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
term in each period, with its standard error from the covariance of the
coefficients (\\\mathrm{Var}(\hat\beta_t) = A V A^\top\\, where \\A\\
maps the basis coefficients to the period effects), and pointwise
confidence intervals. The effects are the same for every basis; only the
coefficients differ.

## Usage

``` r
temporal_effects(object, term = NULL, level = 0.95, vcov = NULL, ...)
```

## Arguments

- object:

  a fitted model with
  [`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
  terms.

- term:

  optional names of the varying variables (`"(Intercept)"` for
  `tv(1, ...)`); all by default.

- level:

  confidence level.

- vcov:

  optional covariance matrix of the coefficients, such as a
  cluster-robust one; `vcov(object)` by default.

- ...:

  unused.

## Value

A data frame with columns `part` (`"mean"` for the mean, otherwise the
coefficient block, such as `"ln(alpha)"`), `term`, `time`, `basis`,
`period`, `estimate`, `std.error`, `conf.low` and `conf.high`.

## See also

[`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md),
[`temporal_tests()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_tests.md)
