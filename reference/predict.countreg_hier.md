# Predictions from a hierarchical count model

Predictions from a hierarchical count model

## Usage

``` r
# S3 method for class 'countreg_hier'
predict(
  object,
  newdata = NULL,
  type = c("response", "link", "prob"),
  re.form = NULL,
  allow.new.levels = TRUE,
  marginal = FALSE,
  history = NULL,
  new_data = NULL,
  counts = NULL,
  max_count = NULL,
  overflow = FALSE,
  underflow = FALSE,
  interval = c("none", "confidence", "prediction"),
  level = 0.95,
  ...
)
```

## Arguments

- object:

  a
  [`countreg_hier()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_hier.md)
  fit.

- newdata, new_data:

  optional data frame, under either name (default: the estimation rows).

- type:

  `"response"` (expected count), `"link"` (log of the mean parameter) or
  `"prob"` (a matrix of count probabilities, given the estimated
  effects; with `marginal` the other effects are integrated out).

- re.form:

  which random effects to condition on: `NULL` (all, the default), `NA`
  (none, population level) or a formula of terms, such as
  `~ (1 | county)`.

- allow.new.levels:

  allow groups not seen in estimation; their effects are set to zero (or
  integrated out with `marginal`).

- marginal:

  integrate the random effects that are not conditioned on (those left
  out by `re.form`, and new groups) out of the expected count, instead
  of setting them to zero. For the log link the marginal mean is
  \\\exp(\eta + z'\Sigma z / 2)\\, larger than \\\exp(\eta)\\.

- history:

  optional data frame of past counts (with the outcome) for the groups
  to predict. The effects of the groups are then estimated from
  `history` alone, at the fitted parameters, instead of taken from the
  estimation data, so a forecast never uses its own outcomes. Groups
  without history are treated as new groups.

- counts, max_count:

  for `"prob"`: the counts, such as `0:10` or `5:8`, or the largest one
  of `0:max_count` (default: the largest training count). For intervals,
  `max_count` is the largest count evaluated.

- overflow, underflow:

  for `"prob"`: add columns with the probabilities above and below the
  counts.

- interval, level:

  confidence interval for the expected count or prediction interval for
  a new count (a data frame with `fit`, `lwr` and `upr`). They include
  the uncertainty of the fixed coefficients, the distribution parameters
  and the estimated effects of known groups, but not of the
  random-effect variances, which are held at their estimates. Prediction
  intervals use `marginal = TRUE` unless it is given. See the section
  *Intervals of other estimators* of
  [`predict.flexCountReg`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md).

- ...:

  `ndraws` and `seed` for intervals.

## Value

A numeric vector with one value per row of `newdata`, a matrix
(`"prob"`) or a data frame (intervals).

## Details

Conditional predictions use the conditional modes of the random effects
(as `lme4` does; posterior means with `engine = "msl"`). With
`marginal = FALSE` (the default, as in `lme4`) the effects that are not
conditioned on are zero, so a population-level prediction is the count
of a median group, not the mean over groups.
