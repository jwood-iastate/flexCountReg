# Predictions from a hierarchical count model

Predictions from a hierarchical count model

## Usage

``` r
# S3 method for class 'countreg_hier'
predict(
  object,
  newdata = NULL,
  type = c("response", "link"),
  re.form = NULL,
  allow.new.levels = TRUE,
  marginal = FALSE,
  history = NULL,
  new_data = NULL,
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

  `"response"` (expected count) or `"link"` (log of the mean parameter).

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

- ...:

  unused.

## Value

A numeric vector with one value per row of `newdata`.

## Details

Conditional predictions use the conditional modes of the random effects
(as `lme4` does; posterior means with `engine = "msl"`). With
`marginal = FALSE` (the default, as in `lme4`) the effects that are not
conditioned on are zero, so a population-level prediction is the count
of a median group, not the mean over groups.
