# Predictions from total-plus-allocation models

Predictions from total-plus-allocation models

## Usage

``` r
# S3 method for class 'countreg_alloc'
predict(
  object,
  newdata = NULL,
  type = c("counts", "prob", "total", "total_positive"),
  conditional = FALSE,
  ...
)
```

## Arguments

- object:

  a
  [`countreg_alloc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_alloc.md)
  fit.

- newdata:

  optional new data; the estimation data by default.

- type:

  `"counts"` (default): the expected count in each category,
  \\E\[c\]\\\pi_j\\; `"prob"`: the allocation probabilities \\\pi_j\\;
  `"total"`: the expected total \\E\[c\]\\; `"total_positive"`: the
  expected total given at least one event, \\E\[c\] / (1 - P(0))\\.

- conditional:

  for `type = "counts"`: `TRUE` gives the expected counts given at least
  one event, \\E\[c \mid c \> 0\]\\\pi_j\\. The default is `FALSE`, the
  unconditional expected counts, also for a fit with `truncated = TRUE`
  (whose estimated distribution is that of the untruncated total).

- ...:

  unused.

## Value

A matrix with one column per category (`"counts"`, `"prob"`), or a
vector (`"total"`, `"total_positive"`).
