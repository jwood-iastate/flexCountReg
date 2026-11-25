# Predictions from total-plus-allocation models

Predictions from total-plus-allocation models

## Usage

``` r
# S3 method for class 'countreg_alloc'
predict(
  object,
  newdata = NULL,
  type = c("counts", "prob", "total", "total_positive", "total_prob", "count_prob"),
  conditional = FALSE,
  random = c("marginal", "posterior"),
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
  [`countreg_alloc()`](https://jwood-iastate.github.io/flexCountReg/reference/countreg_alloc.md)
  fit.

- newdata, new_data:

  optional new data (either name); the estimation data by default.

- type:

  `"counts"` (default): the expected count in each category,
  \\E\[c\]\\\pi_j\\; `"prob"`: the allocation probabilities \\\pi_j\\;
  `"total"`: the expected total \\E\[c\]\\; `"total_positive"`: the
  expected total given at least one event, \\E\[c\] / (1 - P(0))\\;
  `"total_prob"`: a matrix of the probabilities of the total's counts;
  `"count_prob"`: a list with one such matrix per category, the total
  thinned by the category's probability (not for the fractional
  quasi-likelihood, which defines no distribution of the category
  counts).

- conditional:

  for the probabilities, `TRUE` conditions on at least one event; for
  `type = "counts"`: `TRUE` gives the expected counts given at least one
  event, \\E\[c \mid c \> 0\]\\\pi_j\\. The default is `FALSE`, the
  unconditional expected counts, also for a fit with `truncated = TRUE`
  (whose estimated distribution is that of the untruncated total).

- random:

  for fits with random parameters: `"marginal"` (default) averages over
  the distribution of the random parameters (the prediction for a new
  unit); `"posterior"` averages over each unit's posterior distribution
  given its counts in the estimation data (Empirical Bayes; units not in
  the estimation data get the marginal prediction). `newdata` then needs
  the `group_var` column(s). The shared site-year term
  (`shared_formula`) is always averaged over its distribution, as it is
  new for every row.

- counts, max_count:

  for the probabilities: the counts, such as `0:10` or `5:8`, or the
  largest one of `0:max_count` (default: the largest observed total).
  For intervals, `max_count` is the largest count evaluated.

- overflow, underflow:

  add columns with the probabilities above and below the counts.

- interval, level:

  confidence intervals for the expected count or prediction intervals
  for a new count: for `type = "counts"` a list of data frames (`fit`,
  `lwr`, `upr`) named by category, for `type = "total"` one data frame.
  See the section *Intervals of other estimators* of
  [`predict.flexCountReg()`](https://jwood-iastate.github.io/flexCountReg/reference/predict.flexCountReg.md).
  With random parameters, prediction intervals use 50 parameter draws by
  default and are slow.

- ...:

  `ndraws` and `seed` for intervals.

## Value

A matrix with one column per category (`"counts"`, `"prob"`), or a
vector (`"total"`, `"total_positive"`). With random parameters, `"prob"`
is the average allocation probability and `"counts"` the average
expected count over the draws. For a fit with `truncated = TRUE` and
random parameters, the likelihood truncates the total given each unit's
coefficients \\b\\, so `"total_positive"` (and `conditional = TRUE`)
average the ratio over the draws, \\E_b\[E(c \mid b) / (1 - P_0(b))\]\\;
without truncation they give \\E_b\[E(c \mid b)\] / (1 -
E_b\[P_0(b)\])\\. The shared site-year term is integrated out before the
ratio is taken, as in the likelihood. Posterior predictions
(`random = "posterior"`) need random parameters.
