# Predictions from multivariate count models

Predictions from multivariate count models

## Usage

``` r
# S3 method for class 'countreg_mv'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "variance", "prob"),
  max_count = NULL,
  ...
)
```

## Arguments

- object:

  a
  [`countreg_mv()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_mv.md)
  fit.

- newdata:

  optional new data; the estimation data by default.

- type:

  `"response"` (default): a matrix of the expected counts, one column
  per outcome; `"link"`: the linear predictors (log of the mean
  parameter); `"variance"`: the variances; `"prob"`: a list with one
  matrix per outcome of the probabilities of the counts `0:max_count`.

- max_count:

  for `type = "prob"`: the largest count (default: the largest observed
  count of each outcome).

- ...:

  unused.

## Value

A matrix with one column per outcome, or a list of matrices (`"prob"`).
