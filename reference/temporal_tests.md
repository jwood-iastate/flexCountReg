# Tests of time-varying effects

Wald tests of whether the effect of each
[`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
term changes over time: a joint test that the effects are equal in all
periods, and a test of each change between adjacent periods. The tests
concern the period effects, so they are the same for every basis.

## Usage

``` r
temporal_tests(object, term = NULL, vcov = NULL, adjust = "none", ...)
```

## Arguments

- object:

  a fitted model with
  [`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md)
  terms.

- term:

  optional names of the varying variables (`"(Intercept)"` for
  `tv(1, ...)`); all by default.

- vcov:

  optional covariance matrix of the coefficients, such as a
  cluster-robust one; `vcov(object)` by default.

- adjust:

  p-value adjustment for the adjacent-period tests (see
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html)); none by
  default.

- ...:

  unused.

## Value

A data frame with columns `part`, `term`, `test` (`"all equal"` or the
pair of periods), `estimate` (the change, for adjacent periods),
`std.error`, `statistic` (\\\chi^2\\ for the joint test, \\z\\
otherwise), `df` and `p.value`.

## See also

[`tv()`](https://jwood-iastate.github.io/flexCountReg/reference/tv.md),
[`temporal_effects()`](https://jwood-iastate.github.io/flexCountReg/reference/temporal_effects.md)
