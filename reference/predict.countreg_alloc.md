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
  random = c("marginal", "posterior"),
  new_data = NULL,
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
  expected total given at least one event, \\E\[c\] / (1 - P(0))\\.

- conditional:

  for `type = "counts"`: `TRUE` gives the expected counts given at least
  one event, \\E\[c \mid c \> 0\]\\\pi_j\\. The default is `FALSE`, the
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

- ...:

  unused.

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
