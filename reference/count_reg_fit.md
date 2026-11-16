# Fit wrapper used by the parsnip engine

Internal function called by parsnip; use
[`count_reg`](https://jwood-iastate.github.io/flexCountReg/reference/count_reg.md)
or
[`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
instead.

## Usage

``` r
count_reg_fit(formula, data, family = "NB2", weights = NULL, ...)
```

## Arguments

- formula, data, family, weights, ...:

  passed to
  [`countreg`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md).

## Value

A
[`countreg()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg.md)
model.
